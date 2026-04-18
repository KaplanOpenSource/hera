# Hera Makefile
# Usage: make <target>

MONGO_CONTAINER = hera-mongo
MONGO_IMAGE     = mongo:5.0
MONGO_DATA      ?= $(HOME)/mongo-db-datadir
MONGO_INIT      = $(CURDIR)/mongo-init.d
MONGO_PORT      = 27017

SERVER_CONTAINER = hera-server-instance
SERVER_IMAGE     = hera-server

TEST_HERA       ?= $(HOME)/hera_unittest_data

VENV_DIR ?= $(HOME)/.pyhera/environment

.PHONY: install install-env setup install-docs install-rag populate populate-project \
        help mongo-up mongo-down mongo-status mongo-logs mongo-clean \
        build run stop test test-setup \
        install-deps install-deps-all install-paraview install-freecad install-openfoam \
        mermaid-pull render-diagrams render-diagrams-force \
        docs-deps docs-build docs-build-strict docs-serve docs-deploy docs-clean \
        rag-index rag-reindex rag-index-docs rag-search rag-search-raw \
        rag-serve rag-serve-watch rag-watch rag-docs-serve rag-docs-build

# --- Top-level install targets ---

install: mongo-up
	@echo ""
	@echo "=== Hera installed ==="
	@echo "  MongoDB running on port $(MONGO_PORT)"
	@echo "  Run 'make populate' to load repositories into all projects."
	@echo "  Run 'make help' for all available targets."

populate:
	@echo "Populating all projects with registered repositories..."
	hera-project project populate --overwrite
	@echo "Done."

populate-project:
	@echo "Populating project $(PROJECT)..."
	hera-project project populate --projectName $(PROJECT) --overwrite
	@echo "Done."

install-env:
	@echo "Creating Python virtual environment at $(VENV_DIR)..."
	@mkdir -p $(HOME)/.pyhera
	python3 -m venv $(VENV_DIR)
	@echo "Pinning build tools inside the venv (needed for --no-build-isolation)..."
	$(VENV_DIR)/bin/pip install --upgrade pip setuptools wheel
	@echo "Installing requirements (this may take a while)..."
	$(VENV_DIR)/bin/pip install --no-build-isolation -r $(CURDIR)/requirements.txt
	$(VENV_DIR)/bin/pip install --no-build-isolation -e $(CURDIR)
	@echo ""
	@echo "=== Python environment installed at $(VENV_DIR) ==="
	@echo "  Activate with: make setup"
	@echo "  Or manually:   source $(VENV_DIR)/bin/activate"
	@echo ""
	@echo "  GDAL is not installed automatically (its version is bound to"
	@echo "  the system libgdal). After activating the venv, run:"
	@echo "    make install-deps"

setup:
	@if [ -n "$$VIRTUAL_ENV" ]; then \
		echo "Deactivating current environment: $$VIRTUAL_ENV"; \
		echo "deactivate" ; \
	fi
	@if [ ! -f "$(VENV_DIR)/bin/activate" ]; then \
		echo "ERROR: Hera environment not found at $(VENV_DIR)"; \
		echo "Run 'make install-env' first to create it."; \
		exit 1; \
	fi
	@echo ""
	@echo "To activate the Hera environment, run:"
	@echo ""
	@echo "  source $(VENV_DIR)/bin/activate"
	@echo ""
	@echo "(make cannot activate a venv in your current shell — use source directly,"
	@echo " or run: source activate_hera.sh)"

install-docs: docs-build
	@echo ""
	@echo "=== Documentation built ==="
	@echo "  Site in site/ — run 'make docs-serve' to preview."

install-rag: rag-setup
	@echo ""
	@echo "=== RAG installed ==="
	@echo "  Run 'make rag-search RAG_QUERY=\"your question\"'"
	@echo "  Run 'make rag-serve' for REST API."

help:
	@echo "Hera targets:"
	@echo ""
	@echo "  Quick start:"
	@echo "    make install-env         Create Python venv and install all dependencies"
	@echo "    make setup               Show how to activate the Hera environment"
	@echo "    make install             Set up MongoDB (Docker)"
	@echo "    make populate            Load all repositories into all projects"
	@echo "    make populate-project PROJECT=name  Populate a single project"
	@echo "    make install-docs        Build documentation (pulls mermaid, renders diagrams)"
	@echo "    make install-rag         Full RAG setup (services + model + index)"
	@echo ""
	@echo "  MongoDB:"
	@echo "    make mongo-up      Start MongoDB container (data at $(MONGO_DATA))"
	@echo "    make mongo-down    Stop MongoDB container"
	@echo "    make mongo-status  Show MongoDB container status"
	@echo "    make mongo-logs    Tail MongoDB container logs"
	@echo "    make mongo-clean   Stop container and delete data volume"
	@echo ""
	@echo "  Hera Server:"
	@echo "    make build         Build the hera-server Docker image"
	@echo "    make run           Run the hera-server container"
	@echo "    make stop          Stop the hera-server container"
	@echo ""
	@echo "  Testing:"
	@echo "    make test-setup    Create the test data directory structure"
	@echo "    make test          Run all tests (requires MongoDB + test data)"
	@echo ""
	@echo "  Third-party dependencies:"
	@echo "    make install-deps        Install system packages and GDAL Python binding"
	@echo "    make install-paraview    Download and set up ParaView"
	@echo "    make install-freecad     Install FreeCad Python3 bindings"
	@echo "    make install-openfoam    Install OpenFOAM 10"
	@echo "    make install-deps-all    Install all of the above"
	@echo ""
	@echo "  Documentation:"
	@echo "    make docs-build          Build full site (render diagrams + mkdocs build)"
	@echo "    make docs-serve          Start local MkDocs preview server"
	@echo "    make docs-deploy         Deploy to GitHub Pages"
	@echo "    make docs-clean          Remove built site/"
	@echo "    make render-diagrams     Render new/changed Mermaid diagrams to SVG"
	@echo "    make render-diagrams-force  Re-render ALL diagrams"
	@echo "    make mermaid-pull        Pull the mermaid-cli Docker image"
	@echo "    make docs-deps           Install documentation dependencies"
	@echo ""
	@echo "  RAG Search:"
	@echo "    make rag-setup           Full setup: install + services + model + index"
	@echo "    make rag-services-up     Start Qdrant + Cassandra + Ollama"
	@echo "    make rag-services-down   Stop all RAG services"
	@echo "    make rag-services-status Show RAG service status"
	@echo "    make rag-index           Build RAG index from docs/ + hera/"
	@echo "    make rag-reindex         Wipe and rebuild RAG index"
	@echo "    make rag-search          Search with default query"
	@echo "    make rag-serve           Start RAG REST API server"
	@echo "    make rag-serve-watch     Serve + auto re-index on changes"
	@echo "    make rag-clean           Stop services + remove all data"

# --- MongoDB ---

mongo-up:
	@mkdir -p $(MONGO_DATA)
	@if docker ps --format '{{.Names}}' | grep -q '^$(MONGO_CONTAINER)$$'; then \
		echo "$(MONGO_CONTAINER) is already running."; \
	elif docker ps -a --format '{{.Names}}' | grep -q '^$(MONGO_CONTAINER)$$'; then \
		echo "Starting existing $(MONGO_CONTAINER) container..."; \
		docker start $(MONGO_CONTAINER); \
	else \
		echo "Creating and starting $(MONGO_CONTAINER)..."; \
		docker run --name $(MONGO_CONTAINER) \
			-v $(MONGO_DATA):/data/db \
			-v $(MONGO_INIT):/docker-entrypoint-initdb.d \
			-p 127.0.0.1:$(MONGO_PORT):27017 \
			-d $(MONGO_IMAGE); \
	fi
	@echo "Waiting for MongoDB to be ready..."
	@for i in $$(seq 1 30); do \
		if docker exec $(MONGO_CONTAINER) mongosh --quiet --eval 'db.runCommand({ping:1})' >/dev/null 2>&1; then \
			echo "MongoDB is ready on port $(MONGO_PORT)."; \
			break; \
		fi; \
		if [ "$$i" -eq 30 ]; then \
			echo "ERROR: MongoDB did not become ready in time."; \
			exit 1; \
		fi; \
		sleep 1; \
	done

mongo-down:
	@if docker ps --format '{{.Names}}' | grep -q '^$(MONGO_CONTAINER)$$'; then \
		docker stop $(MONGO_CONTAINER); \
		echo "$(MONGO_CONTAINER) stopped."; \
	else \
		echo "$(MONGO_CONTAINER) is not running."; \
	fi

mongo-status:
	@if docker ps --format '{{.Names}}' | grep -q '^$(MONGO_CONTAINER)$$'; then \
		echo "$(MONGO_CONTAINER) is running."; \
		docker ps --filter name=$(MONGO_CONTAINER) --format 'table {{.Status}}\t{{.Ports}}'; \
	elif docker ps -a --format '{{.Names}}' | grep -q '^$(MONGO_CONTAINER)$$'; then \
		echo "$(MONGO_CONTAINER) exists but is stopped."; \
	else \
		echo "$(MONGO_CONTAINER) does not exist."; \
	fi

mongo-logs:
	docker logs -f $(MONGO_CONTAINER)

mongo-clean:
	@echo "This will stop the container and DELETE all MongoDB data at $(MONGO_DATA)."
	@read -p "Are you sure? [y/N] " confirm && [ "$$confirm" = "y" ] || (echo "Aborted."; exit 1)
	-docker rm -f $(MONGO_CONTAINER) 2>/dev/null
	rm -rf $(MONGO_DATA)
	@echo "Container removed and data deleted."

# --- Hera Server ---

build:
	docker build -t $(SERVER_IMAGE) -f dockerfile .

run:
	docker run -it --network host -v .:/app --rm --name $(SERVER_CONTAINER) $(SERVER_IMAGE)

stop:
	-docker stop $(SERVER_CONTAINER)

# --- Tests ---

test-setup:
	@echo "Setting up test directory at $(TEST_HERA)..."
	mkdir -p $(TEST_HERA)/measurements/GIS/raster
	mkdir -p $(TEST_HERA)/measurements/GIS/vector
	mkdir -p $(TEST_HERA)/measurements/meteorology/lowfreqdata
	mkdir -p $(TEST_HERA)/measurements/meteorology/highfreqdata
	mkdir -p $(TEST_HERA)/expected/BASELINE
	@if [ ! -f $(TEST_HERA)/data_config.json ]; then \
		echo '{"default_result_set": "BASELINE"}' > $(TEST_HERA)/data_config.json; \
		echo "Created data_config.json"; \
	else \
		echo "data_config.json already exists, skipping."; \
	fi
	@if [ ! -f $(TEST_HERA)/test_repository.json ]; then \
		echo '{}' > $(TEST_HERA)/test_repository.json; \
		echo "Created test_repository.json"; \
	else \
		echo "test_repository.json already exists, skipping."; \
	fi
	@echo ""
	@echo "Test directory ready at $(TEST_HERA)"
	@echo "Add test data files under measurements/ and expected outputs under expected/BASELINE/"

test:
	TEST_HERA=$(TEST_HERA) pytest hera/tests/ -v

# --- Third-party Dependencies ---

install-deps:
	@echo "Installing system packages..."
	sudo apt-get update
	sudo apt-get install -y libcairo2-dev pkg-config python3-dev \
		libgirepository1.0-dev libgdal-dev gdal-bin python3-gdal
	@echo "Installing GDAL Python binding matching system version..."
	pip install GDAL==$$(gdal-config --version)
	@echo "System dependencies installed."

install-paraview:
	@echo "Installing ParaView..."
	@if [ -d "$(HOME)/paraview" ]; then \
		echo "ParaView directory already exists at $(HOME)/paraview, skipping download."; \
	else \
		echo "Downloading ParaView 5.11.0..."; \
		wget -q --show-progress -O /tmp/paraview.tar.gz \
			"https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.11&type=binary&os=Linux&downloadFile=ParaView-5.11.0-MPI-Linux-Python3.9-x86_64.tar.gz"; \
		mkdir -p $(HOME)/paraview; \
		tar -xzf /tmp/paraview.tar.gz -C $(HOME)/paraview --strip-components=1; \
		rm -f /tmp/paraview.tar.gz; \
	fi
	@echo ""
	@echo "ParaView installed at $(HOME)/paraview"
	@echo "Add to your environment:"
	@echo "  export PYTHONPATH=$(HOME)/paraview/lib/python3.9/site-packages:\$$PYTHONPATH"

install-freecad:
	@echo "Installing FreeCad Python3 bindings..."
	sudo apt-get update
	sudo apt-get install -y libfreecad-python3-0.19
	@echo ""
	@echo "FreeCad installed. Add to your code or environment:"
	@echo "  export PYTHONPATH=/usr/lib/freecad-python3/lib/:\$$PYTHONPATH"

install-openfoam:
	@echo "Installing OpenFOAM 10..."
	sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
	sudo add-apt-repository -y http://dl.openfoam.org/ubuntu
	sudo apt-get update
	sudo apt-get install -y openfoam10
	@echo ". /opt/openfoam10/etc/bashrc" > $(CURDIR)/of10
	@echo ""
	@echo "OpenFOAM 10 installed."
	@echo "Run 'source of10' to set up the OpenFOAM environment."

install-deps-all: install-deps install-paraview install-freecad install-openfoam
	@echo ""
	@echo "All third-party dependencies installed."

# --- Documentation ---

mermaid-pull:
	@echo "Pulling mermaid-cli Docker image..."
	docker pull minlag/mermaid-cli
	@echo "Done. Run 'make render-diagrams' to render all Mermaid diagrams."

render-diagrams: mermaid-pull
	@echo "Rendering Mermaid diagrams to SVG (skips unchanged)..."
	python render_diagrams.py
	@echo "Done. SVGs saved to docs/images/diagrams/"

render-diagrams-force: mermaid-pull
	@echo "Re-rendering ALL Mermaid diagrams to SVG..."
	python render_diagrams.py --force
	@echo "Done. All SVGs regenerated."

docs-deps:
	@echo "Installing documentation dependencies..."
	pip install -r docs/requirements-docs.txt
	pip install -e . --no-deps
	@echo "Done."

docs-build: render-diagrams docs-deps
	@echo "Building MkDocs site..."
	mkdocs build
	@echo "Done. Static site built in site/"

docs-build-strict: render-diagrams docs-deps
	@echo "Building MkDocs site (strict mode)..."
	mkdocs build --strict
	@echo "Done. Static site built in site/"

docs-serve: render-diagrams
	@echo "Starting MkDocs development server at http://127.0.0.1:8000"
	./serve_docs.sh

docs-deploy: render-diagrams docs-deps
	@echo "Deploying documentation to GitHub Pages..."
	mkdocs gh-deploy --force
	@echo "Done. Site deployed."

docs-clean:
	@echo "Removing built site..."
	rm -rf site/
	@echo "Done."

# --- RAG ---

RAG_DOCS  ?= docs
RAG_CODE  ?= hera
RAG_QUERY ?= "How do I get started?"

rag-index:
	hera-rag-search index --docs $(RAG_DOCS) --code $(RAG_CODE)

rag-reindex:
	hera-rag-search index --docs $(RAG_DOCS) --code $(RAG_CODE) --clean

rag-index-docs:
	hera-rag-search index --docs $(RAG_DOCS) --docs-only

rag-search:
	hera-rag-search search $(RAG_QUERY)

rag-search-raw:
	hera-rag-search search $(RAG_QUERY) --raw

rag-serve:
	hera-rag-search serve

rag-serve-watch:
	hera-rag-search serve --with-watcher --watch-root $(RAG_DOCS)

rag-watch:
	hera-rag-search watch --root $(RAG_DOCS)

rag-docs-serve:
	RAG_ENABLED=true hera-rag-search serve &
	sleep 2
	RAG_ENABLED=true mkdocs serve

rag-docs-build:
	RAG_ENABLED=true mkdocs build

# ── RAG Infrastructure ────────────────────────────────────────────────────────

QDRANT_CONTAINER = hera-qdrant
QDRANT_IMAGE     = qdrant/qdrant
QDRANT_PORT      = 6333
QDRANT_DATA      ?= $(HOME)/qdrant-data

CASS_CONTAINER   = hera-cassandra
CASS_IMAGE       = cassandra:4
CASS_PORT        = 9042
CASS_DATA        ?= $(HOME)/cassandra-data

OLLAMA_CONTAINER = hera-ollama
OLLAMA_IMAGE     = ollama/ollama
OLLAMA_PORT      = 11434
OLLAMA_DATA      ?= $(HOME)/ollama-data
OLLAMA_MODEL     ?= llama3

rag-install:
	@echo "Installing RAG dependencies..."
	pip install -e .[rag]
	@echo "Done."

rag-qdrant-up:
	@mkdir -p $(QDRANT_DATA)
	@if docker ps --format '{{.Names}}' | grep -q '^$(QDRANT_CONTAINER)$$'; then \
		echo "$(QDRANT_CONTAINER) is already running."; \
	elif docker ps -a --format '{{.Names}}' | grep -q '^$(QDRANT_CONTAINER)$$'; then \
		echo "Starting existing $(QDRANT_CONTAINER)..."; \
		docker start $(QDRANT_CONTAINER); \
	else \
		echo "Creating $(QDRANT_CONTAINER)..."; \
		docker run --name $(QDRANT_CONTAINER) \
			-v $(QDRANT_DATA):/qdrant/storage \
			-p 127.0.0.1:$(QDRANT_PORT):6333 \
			-d $(QDRANT_IMAGE); \
	fi
	@echo "Qdrant ready on port $(QDRANT_PORT)."

rag-qdrant-down:
	@docker stop $(QDRANT_CONTAINER) 2>/dev/null || true
	@echo "$(QDRANT_CONTAINER) stopped."

rag-cassandra-up:
	@mkdir -p $(CASS_DATA)
	@if docker ps --format '{{.Names}}' | grep -q '^$(CASS_CONTAINER)$$'; then \
		echo "$(CASS_CONTAINER) is already running."; \
	elif docker ps -a --format '{{.Names}}' | grep -q '^$(CASS_CONTAINER)$$'; then \
		echo "Starting existing $(CASS_CONTAINER)..."; \
		docker start $(CASS_CONTAINER); \
	else \
		echo "Creating $(CASS_CONTAINER)..."; \
		docker run --name $(CASS_CONTAINER) \
			-v $(CASS_DATA):/var/lib/cassandra \
			-p 127.0.0.1:$(CASS_PORT):9042 \
			-d $(CASS_IMAGE); \
	fi
	@echo "Waiting for Cassandra to be ready..."
	@for i in $$(seq 1 60); do \
		if docker exec $(CASS_CONTAINER) cqlsh -e "DESCRIBE KEYSPACES" >/dev/null 2>&1; then \
			echo "Cassandra ready on port $(CASS_PORT)."; \
			break; \
		fi; \
		if [ "$$i" -eq 60 ]; then \
			echo "WARNING: Cassandra may still be starting up."; \
		fi; \
		sleep 2; \
	done

rag-cassandra-down:
	@docker stop $(CASS_CONTAINER) 2>/dev/null || true
	@echo "$(CASS_CONTAINER) stopped."

rag-ollama-up:
	@mkdir -p $(OLLAMA_DATA)
	@if docker ps --format '{{.Names}}' | grep -q '^$(OLLAMA_CONTAINER)$$'; then \
		echo "$(OLLAMA_CONTAINER) is already running."; \
	elif docker ps -a --format '{{.Names}}' | grep -q '^$(OLLAMA_CONTAINER)$$'; then \
		echo "Starting existing $(OLLAMA_CONTAINER)..."; \
		docker start $(OLLAMA_CONTAINER); \
	else \
		echo "Creating $(OLLAMA_CONTAINER)..."; \
		docker run --name $(OLLAMA_CONTAINER) \
			-v $(OLLAMA_DATA):/root/.ollama \
			-p 127.0.0.1:$(OLLAMA_PORT):11434 \
			-d $(OLLAMA_IMAGE); \
	fi
	@echo "Ollama ready on port $(OLLAMA_PORT)."

rag-ollama-down:
	@docker stop $(OLLAMA_CONTAINER) 2>/dev/null || true
	@echo "$(OLLAMA_CONTAINER) stopped."

rag-ollama-pull:
	@echo "Pulling model $(OLLAMA_MODEL)..."
	@docker exec $(OLLAMA_CONTAINER) ollama pull $(OLLAMA_MODEL)
	@echo "Model $(OLLAMA_MODEL) ready."

rag-services-up: rag-qdrant-up rag-cassandra-up rag-ollama-up
	@echo "All RAG services running."

rag-services-down: rag-qdrant-down rag-cassandra-down rag-ollama-down
	@echo "All RAG services stopped."

rag-services-status:
	@echo "=== RAG Services ==="
	@docker ps --filter name=hera-qdrant --filter name=hera-cassandra --filter name=hera-ollama \
		--format 'table {{.Names}}\t{{.Status}}\t{{.Ports}}' 2>/dev/null || echo "No RAG containers found."

rag-setup: rag-install rag-services-up
	@echo "Waiting for services to stabilize..."
	@sleep 5
	@$(MAKE) rag-ollama-pull
	@$(MAKE) rag-reindex
	@echo ""
	@echo "=== RAG Setup Complete ==="
	@echo "  make rag-search RAG_QUERY=\"your question\""
	@echo "  make rag-serve"

rag-clean: rag-services-down
	@echo "Removing RAG containers and data..."
	-docker rm -f $(QDRANT_CONTAINER) $(CASS_CONTAINER) $(OLLAMA_CONTAINER) 2>/dev/null
	rm -rf $(QDRANT_DATA) $(CASS_DATA) $(OLLAMA_DATA)
	@echo "RAG data cleaned."

.PHONY: rag-install rag-qdrant-up rag-qdrant-down rag-cassandra-up rag-cassandra-down \
        rag-ollama-up rag-ollama-down rag-ollama-pull \
        rag-services-up rag-services-down rag-services-status rag-setup rag-clean
