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

.PHONY: help mongo-up mongo-down mongo-status mongo-logs mongo-clean \
        build run stop test test-setup

help:
	@echo "Hera targets:"
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
