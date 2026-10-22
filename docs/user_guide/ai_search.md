# AI-Powered Documentation Search

Hera includes an optional AI search feature that lets you ask natural language questions about the documentation and source code. It uses a self-hosted LLM (via Ollama) to generate answers based on the indexed content.

The entire pipeline runs locally — no data is sent to external services.

---

## How it works

When enabled, a floating "Ask AI" button appears on every documentation page. Click it to open a search panel where you can type questions in plain English. The system:

1. Searches the indexed docs and docstrings for relevant content (via Qdrant vector search)
2. Fetches full text from Cassandra (source of truth)
3. Sends the top matches + your question to a local LLM (Ollama)
4. Streams the answer back in real-time
5. Shows source references with relevance scores

---

## Quick start (all-in-one)

If you just want everything running with one command:

```bash
make rag-setup
make rag-docs-serve
```

`make rag-setup` installs the `rag` extra, starts all services, downloads the LLM model and builds the index (clean rebuild). `make rag-docs-serve` then starts the API in the background and serves the docs with `RAG_ENABLED=true` — which the settings ignore (see the note below); export `RAG_RAG_ENABLED=true` yourself.

!!! warning "The widget needs the MkDocs plugin registered"
    The widget is injected by the `hera_rag_search` MkDocs plugin (`hera.utils.rag.serve:HeraMkDocsPlugin`, registered as an entry point in `setup.py`). The repository's `mkdocs.yml` does not list it under `plugins:` at the moment, so to see the widget add it yourself:

    ```yaml
    plugins:
      - search
      - hera_rag_search
    ```

For step-by-step control, follow the three phases below.

---

## Phase 1: Set up the RAG environment

### Prerequisites

- **Docker** — all three backend services run as containers
- **~6 GB disk** — for the LLM model (~4 GB) + databases (~2 GB)
- **~4 GB RAM** — for embedding model + Ollama inference

### Start the services

The RAG system uses three backend services:

| Service | Purpose | Default port |
|---------|---------|-------------|
| **Qdrant** | Vector store for fast similarity search | `6333` |
| **Cassandra** | Text store (source of truth for chunk content) | `9042` |
| **Ollama** | Self-hosted LLM inference | `11434` |

Start all three at once:

```bash
make rag-services-up
```

Or start individually:

```bash
make rag-qdrant-up       # Vector store
make rag-cassandra-up    # Text store (takes ~30s to become ready)
make rag-ollama-up       # LLM server
```

### Download the LLM model

```bash
make rag-ollama-pull
```

This downloads the default model (`llama3`, ~4 GB) into the `hera-ollama` container. To use a different model, pull it and point Hera at it with `RAG_OLLAMA_MODEL` (the Makefile variable `OLLAMA_MODEL` only controls what is pulled):

```bash
OLLAMA_MODEL=mistral make rag-ollama-pull
export RAG_OLLAMA_MODEL=mistral
```

If you already run Ollama natively on the host, use `make rag-services-up-native` (Qdrant + Cassandra only) and set `RAG_OLLAMA_BASE_URL` to your Ollama instance.

### Verify services are running

```bash
make rag-services-status
```

You should see three running containers: `hera-qdrant`, `hera-cassandra`, `hera-ollama`.

### Data locations

| Service | Data stored at |
|---------|---------------|
| Qdrant | `~/qdrant-data/` |
| Cassandra | `~/cassandra-data/` |
| Ollama | `~/ollama-data/` |

Override with environment variables: `QDRANT_DATA`, `CASS_DATA`, `OLLAMA_DATA`.

### Configuration

All RAG settings (`hera/utils/rag/config.py`) use the `RAG_` prefix and can be set as environment variables or in a `.env` file in the directory you run from (see `.env.example`):

```bash
# .env (in project root)
RAG_OLLAMA_MODEL=llama3
RAG_OLLAMA_BASE_URL=http://localhost:11434
RAG_QDRANT_HOST=localhost
RAG_QDRANT_PORT=6333
RAG_CASSANDRA_HOSTS=["localhost"]
RAG_CASSANDRA_PORT=9042
RAG_EMBED_MODEL=BAAI/bge-small-en-v1.5
RAG_TOP_K=5
RAG_CHUNK_SIZE=1000
```

Other settings and their defaults:

| Variable | Default | Purpose |
|----------|---------|---------|
| `RAG_RAG_ENABLED` | `false` | Inject the docs widget (MkDocs plugin) |
| `RAG_RAG_API_URL` | `http://localhost:8765` | API URL the widget calls |
| `RAG_API_HOST` / `RAG_API_PORT` | `0.0.0.0` / `8765` | Where `hera-rag-search serve` listens (all interfaces by default) |
| `RAG_OLLAMA_TIMEOUT` | `120` | Ollama request timeout (seconds) |
| `RAG_CHUNK_OVERLAP` | `150` | Overlap of the sliding window for long sections |
| `RAG_MIN_CHUNK_CHARS` | `80` | Chunks shorter than this are skipped |
| `RAG_CASSANDRA_KEYSPACE` | `hera_rag` | Cassandra keyspace (created automatically) |
| `RAG_QDRANT_COLLECTION` | `hera_docs` | Qdrant collection (created automatically) |
| `RAG_EMBED_DIM` | `384` | Embedding size; must match `RAG_EMBED_MODEL` |
| `RAG_MODELS_CACHE_DIR` | `~/.pyhera/rag_models` | Local cache for the embedding model |
| `RAG_OFFLINE_MODE` | `false` | Never download; the embedding model must already be cached |

!!! note "Doubled prefix"
    pydantic-settings prepends `RAG_` to each field name, and these two fields are named `rag_enabled` / `rag_api_url` (`hera/utils/rag/config.py:55-58`), so the variables are `RAG_RAG_ENABLED` and `RAG_RAG_API_URL`. Plain `RAG_ENABLED` (used by the `rag-docs-*` Makefile targets) has no effect.

---

## Phase 2: Set up the RAG index

### Install Python dependencies

```bash
pip install -e .[rag]
```

This installs `sentence-transformers`, `qdrant-client`, `cassandra-driver`, `fastapi`, `typer`, and related packages.

### Build the index

```bash
make rag-index
```

Or via the CLI directly:

```bash
hera-rag-search index --docs docs --code hera
```

### What gets indexed

| Content type | Source | Chunking strategy |
|-------------|--------|------------------|
| **Markdown docs** | `docs/**/*.md`, `docs/**/*.mdx` | Split by headings (H1/H2/H3), sliding window (`RAG_CHUNK_SIZE`, `RAG_CHUNK_OVERLAP`) for long sections |
| **Python docstrings** | `hera/**/*.py` | One chunk per module/class/function docstring (AST-parsed) |
| **Jupyter notebooks** | `docs/**/*.ipynb` | One chunk per cell |

Small chunks (< `RAG_MIN_CHUNK_CHARS`, default 80 characters) are automatically skipped.

### Index output

After indexing, you'll see a summary like:

```
Done.  Markdown: 450  Docstrings: 380  Notebooks: 25  Total: 855
```

### Rebuild options

```bash
# Full rebuild: wipe everything and re-index (hera-rag-search index --clean)
make rag-reindex

# Index without wiping (use on an empty index)
make rag-index

# Docs only (skip Python docstrings)
make rag-index-docs
```

!!! warning "`rag-index` is not incremental"
    `hera-rag-search index` without `--clean` re-indexes every file and stores each chunk under a new random ID, so running it on an existing index duplicates all chunks. Use `make rag-reindex` (or `--clean`) to refresh the index; per-file updates are done by the watcher below.

### Auto re-index on file changes

```bash
make rag-watch
```

This watches `docs/` (`hera-rag-search watch --root docs`) for changes to `.md`, `.mdx`, `.py` and `.ipynb` files and re-indexes each modified file (debounced at 2 seconds); deleted files have their chunks removed. Changes under `hera/` are not watched.

---

## Phase 3: Activate the RAG

### Start the REST API server

```bash
make rag-serve
```

The API runs at `http://localhost:8765`. Verify it's working:

```bash
curl http://localhost:8765/health
# {"status": "ok", "model": "llama3"}
```

### Enable in MkDocs (documentation widget)

```bash
make rag-docs-serve
```

This starts both the RAG API and MkDocs (export `RAG_RAG_ENABLED=true` first; the target's own `RAG_ENABLED=true` is ignored). Once the `hera_rag_search` plugin is listed in `mkdocs.yml` (see the warning above), the floating "Ask AI" button appears in the bottom-right corner of every page.

To enable manually:

```bash
RAG_RAG_ENABLED=true mkdocs serve
```

### Use from the command line

```bash
# Ask a question (LLM generates the answer)
hera-rag-search search "How do I create a project and load a repository?"

# Search without LLM (show matching chunks only)
hera-rag-search search "demography API" --raw

# Filter by content type
hera-rag-search search "addDataSource" --type docstring
hera-rag-search search "project lifecycle" --type markdown

# Filter by source path prefix (applied to the listed chunks, not to the LLM answer)
hera-rag-search search "LSM" --source toolkits/

# Number of chunks (default 5)
hera-rag-search search "LSM" --top-k 10

# Pre-download the embedding model for offline use (RAG_OFFLINE_MODE=true)
hera-rag-search download-models
```

### Use from Python

```python
from hera.utils.rag import RAGSearch

rag = RAGSearch()

# Full answer
answer = rag.ask("How does the risk assessment toolkit work?")
print(answer)

# Streaming answer
for token in rag.stream("Show me the demography API"):
    print(token, end="", flush=True)

# Raw chunks (no LLM)
chunks = rag.retrieve("data source versioning", top_k=10)
for c in chunks:
    print(f"{c['score']:.3f}  {c['source']} § {c['section']}")
```

### Use the REST API directly

```bash
# Search with LLM answer
curl -X POST http://localhost:8765/search \
  -H "Content-Type: application/json" \
  -d '{"query": "How do I create a project?", "top_k": 5}'

# Stream answer (SSE)
curl -N "http://localhost:8765/stream?q=How+does+the+LSM+toolkit+work"

# Chunks only (no LLM)
curl -X POST http://localhost:8765/search \
  -H "Content-Type: application/json" \
  -d '{"query": "How do I create a project?", "no_llm": true}'

# Same via GET (k = number of chunks)
curl "http://localhost:8765/search?q=create+project&k=5&no_llm=true"
```

All search endpoints accept an optional `chunk_type` (`markdown`, `docstring`, `notebook`). The API also has `POST /ingest` (body `{"docs_root": ..., "code_root": ..., "clean": false}`, builds the index) and `DELETE /source` (body `{"source": "<relative path>"}`, removes a file's chunks).

API docs are available at `http://localhost:8765/docs` (OpenAPI/Swagger).

### Combined: serve + auto re-index

```bash
make rag-serve-watch
```

This starts the API server and watches `docs/` for file changes simultaneously (`hera-rag-search serve --with-watcher --watch-root docs`).

---

## Troubleshooting

### Cassandra takes a long time to start

Cassandra can take 30–60 seconds to become ready after the container starts. The `make rag-cassandra-up` target waits automatically. If you see connection errors, wait a minute and retry.

### Ollama model download is slow

The default `llama3` model is ~4 GB. Use a smaller model for faster setup:

```bash
OLLAMA_MODEL=phi3 make rag-ollama-pull
export RAG_OLLAMA_MODEL=phi3
```

### Port conflicts

If ports 6333, 9042, or 11434 are already in use, either stop the conflicting service or move the containers and tell Hera the new ports. The host ports of the containers are Makefile variables, and the client side is configured in `.env`:

```bash
make rag-qdrant-up QDRANT_PORT=16333
make rag-cassandra-up CASS_PORT=19042
make rag-ollama-up OLLAMA_PORT=21434
```

```bash
# .env
RAG_QDRANT_PORT=16333
RAG_CASSANDRA_PORT=19042
RAG_OLLAMA_BASE_URL=http://localhost:21434
```

If the Ollama clash is with a native Ollama install, use `make rag-services-up-native` instead.

### Embedding model download

The sentence-transformers model (`BAAI/bge-small-en-v1.5`, ~130 MB) downloads automatically on first use and is saved under `~/.pyhera/rag_models/` (`RAG_MODELS_CACHE_DIR`). Subsequent runs use the cached copy. To prepare an offline machine, run `make rag-download-models` (or `hera-rag-search download-models`) while online, then set `RAG_OFFLINE_MODE=true`.

### Re-indexing after documentation changes

If search results seem stale, rebuild the index:

```bash
make rag-reindex
```

---

## Cleanup

```bash
# Stop all services (keep data)
make rag-services-down

# Stop services AND delete all data (qdrant, cassandra, ollama)
make rag-clean
```

---

## Disabling AI search

Leave `RAG_RAG_ENABLED` unset or `false` (the default). The docs site works normally without it. With `RAG_RAG_ENABLED=true` the widget is injected regardless of whether the RAG API is running; if the API is down, questions fail with "Error connecting to RAG server."

---

## All Makefile targets

| Target | Purpose |
|--------|---------|
| `make rag-setup` | Full setup: install + services + model + index |
| `make install-rag` | Same as `rag-setup` |
| `make rag-install` | `pip install -e .[rag]` |
| `make rag-services-up` | Start Qdrant + Cassandra + Ollama |
| `make rag-services-up-native` | Start Qdrant + Cassandra only (use host Ollama) |
| `make rag-qdrant-up` / `rag-cassandra-up` / `rag-ollama-up` | Start a single service (matching `-down` targets stop it) |
| `make rag-services-down` | Stop all services |
| `make rag-services-status` | Show service status |
| `make rag-ollama-pull` | Download LLM model (`OLLAMA_MODEL`, default `llama3`) |
| `make rag-download-models` | Cache the embedding model for offline use |
| `make rag-index` | Build search index |
| `make rag-reindex` | Wipe and rebuild index |
| `make rag-index-docs` | Index docs only (skip docstrings) |
| `make rag-search` | Search with default query |
| `make rag-search-raw` | Search without LLM |
| `make rag-serve` | Start REST API server |
| `make rag-serve-watch` | Serve + auto re-index on changes |
| `make rag-watch` | Watch files for auto re-indexing |
| `make rag-docs-serve` | MkDocs + RAG widget enabled |
| `make rag-docs-build` | Build static site + RAG widget |
| `make rag-clean` | Stop services + delete all data |
