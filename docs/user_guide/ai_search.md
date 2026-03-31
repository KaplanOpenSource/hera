# AI-Powered Documentation Search

Hera includes an optional AI search feature that lets you ask natural language questions about the documentation and source code. It uses a self-hosted LLM (via Ollama) to generate answers based on the indexed content.

The entire pipeline runs locally — no data is sent to external services.

---

## How it works

When enabled, a floating search button appears on every documentation page. Click it to open a search panel where you can type questions in plain English. The system:

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

This starts all services, downloads the LLM model, builds the index, starts the API, and serves the docs with the search widget enabled.

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

This downloads the default model (`llama3`, ~4 GB). To use a different model:

```bash
OLLAMA_MODEL=mistral make rag-ollama-pull
```

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

All RAG settings use the `RAG_` prefix and can be set in a `.env` file (see `.env.example`):

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
| **Markdown docs** | `docs/**/*.md` | Split by headings (H1/H2/H3), sliding window for long sections |
| **Python docstrings** | `hera/**/*.py` | One chunk per module/class/function docstring (AST-parsed) |
| **Jupyter notebooks** | `docs/**/*.ipynb` | One chunk per cell |

Small chunks (< 80 characters) are automatically skipped.

### Index output

After indexing, you'll see a summary like:

```
Indexed: 450 doc chunks, 380 docstrings, 25 notebook cells
Total: 855 chunks in Qdrant + Cassandra
```

### Rebuild options

```bash
# Incremental: only index new/changed files
make rag-index

# Full rebuild: wipe everything and re-index
make rag-reindex

# Docs only (skip Python docstrings)
make rag-index-docs
```

### Auto re-index on file changes

```bash
make rag-watch
```

This watches `docs/` and `hera/` for changes and automatically re-indexes modified files (debounced at 2 seconds).

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

This starts both the RAG API and MkDocs with the search widget enabled. The floating "Ask AI" button appears in the bottom-right corner of every page.

To enable manually:

```bash
RAG_ENABLED=true mkdocs serve
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

# Filter by source path
hera-rag-search search "LSM" --source toolkits/
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
```

API docs are available at `http://localhost:8765/docs` (OpenAPI/Swagger).

### Combined: serve + auto re-index

```bash
make rag-serve-watch
```

This starts the API server and watches for file changes simultaneously.

---

## Troubleshooting

### Cassandra takes a long time to start

Cassandra can take 30–60 seconds to become ready after the container starts. The `make rag-cassandra-up` target waits automatically. If you see connection errors, wait a minute and retry.

### Ollama model download is slow

The default `llama3` model is ~4 GB. Use a smaller model for faster setup:

```bash
export RAG_OLLAMA_MODEL=phi3
make rag-ollama-pull
```

### Port conflicts

If ports 6333, 9042, or 11434 are already in use, either stop the conflicting service or change the RAG ports in `.env`:

```bash
RAG_QDRANT_PORT=16333
RAG_CASSANDRA_PORT=19042
RAG_OLLAMA_BASE_URL=http://localhost:21434
```

### Embedding model download

The sentence-transformers model (`BAAI/bge-small-en-v1.5`, ~130 MB) downloads automatically on first use. Subsequent runs use the cached version.

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

Set `RAG_ENABLED=false` (the default) or simply don't start the RAG server. The docs site works normally without it — the widget only appears when both `RAG_ENABLED=true` and the RAG API is accessible.

---

## All Makefile targets

| Target | Purpose |
|--------|---------|
| `make rag-setup` | Full setup: install + services + model + index |
| `make rag-services-up` | Start Qdrant + Cassandra + Ollama |
| `make rag-services-down` | Stop all services |
| `make rag-services-status` | Show service status |
| `make rag-ollama-pull` | Download LLM model |
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
