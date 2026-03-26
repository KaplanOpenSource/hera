# RAG Search Toolkit

The RAG (Retrieval-Augmented Generation) toolkit provides AI-powered search over the Hera documentation and source code. It indexes both markdown docs and Python docstrings into a dual-store system and answers natural language questions using a self-hosted LLM.

---

## Architecture

```
User question
    ↓
Embedding model (sentence-transformers)
    ↓
Qdrant (vector ANN search) → top-K chunk IDs
    ↓
Cassandra (full text lookup by ID) → chunk text + metadata
    ↓
Prompt: context chunks + question
    ↓
Ollama (self-hosted LLM) → answer
```

| Component | Role | Default |
|-----------|------|---------|
| **Qdrant** | Vector similarity search (ANN) | `localhost:6333` |
| **Cassandra** | Full text storage + metadata | `localhost:9042` |
| **Ollama** | Self-hosted LLM inference | `localhost:11434`, model `llama3` |
| **sentence-transformers** | Text → vector embedding | `BAAI/bge-small-en-v1.5` (384 dim) |

### Why dual stores?

- **Qdrant** stores only vectors + lightweight payload (source, section, type) — optimized for fast ANN search
- **Cassandra** stores the full chunk text + rich metadata — source of truth, no vector size limits
- Linked by a shared UUID (`chunk_id`) per chunk

---

## Package structure

```
hera/utils/rag/
    __init__.py          # public API: RAGSearch, retrieve, build_index
    config.py            # pydantic-settings config (all RAG_* env vars)
    indexer.py           # chunkers + dual-store write/delete
    search.py            # RAGSearch class (ask, stream, retrieve)
    serve.py             # FastAPI server + file watcher + MkDocs plugin
```

---

## Configuration (`config.py`)

All settings are overridable via environment variables with the `RAG_` prefix, or via a `.env` file:

| Setting | Env var | Default | Description |
|---------|---------|---------|-------------|
| `cassandra_hosts` | `RAG_CASSANDRA_HOSTS` | `["localhost"]` | Cassandra cluster hosts |
| `cassandra_port` | `RAG_CASSANDRA_PORT` | `9042` | Cassandra native port |
| `cassandra_keyspace` | `RAG_CASSANDRA_KEYSPACE` | `hera_rag` | Keyspace name |
| `qdrant_host` | `RAG_QDRANT_HOST` | `localhost` | Qdrant server host |
| `qdrant_port` | `RAG_QDRANT_PORT` | `6333` | Qdrant gRPC port |
| `qdrant_collection` | `RAG_QDRANT_COLLECTION` | `hera_docs` | Collection name |
| `ollama_base_url` | `RAG_OLLAMA_BASE_URL` | `http://localhost:11434` | Ollama API URL |
| `ollama_model` | `RAG_OLLAMA_MODEL` | `llama3` | LLM model name |
| `embed_model` | `RAG_EMBED_MODEL` | `BAAI/bge-small-en-v1.5` | Embedding model |
| `embed_dim` | `RAG_EMBED_DIM` | `384` | Embedding vector dimension |
| `chunk_size` | `RAG_CHUNK_SIZE` | `1000` | Max characters per chunk |
| `chunk_overlap` | `RAG_CHUNK_OVERLAP` | `150` | Overlap between sliding window chunks |
| `min_chunk_chars` | `RAG_MIN_CHUNK_CHARS` | `80` | Skip chunks shorter than this |
| `top_k` | `RAG_TOP_K` | `5` | Default number of chunks to retrieve |
| `api_port` | `RAG_API_PORT` | `8765` | FastAPI server port |
| `rag_enabled` | `RAG_ENABLED` | `false` | Enable MkDocs widget injection |

See `.env.example` for a complete template.

---

## Indexer (`indexer.py`)

### Chunking strategies

The indexer processes three types of content with different chunking strategies:

| Content type | Source | Chunking method |
|-------------|--------|----------------|
| **Markdown** | `docs/**/*.md` | Split by `#`/`##`/`###` headings, then sliding window if section exceeds `chunk_size` |
| **Python docstrings** | `hera/**/*.py` | AST parsing → extract module docstring + each class/function docstring |
| **Jupyter notebooks** | `docs/**/*.ipynb` | One chunk per cell (code or markdown) |

### Chunk schema

Each chunk stored in Cassandra has:

| Field | Type | Description |
|-------|------|-------------|
| `id` | UUID | Shared with Qdrant point ID |
| `text` | text | Full chunk content |
| `source` | text | Relative file path (e.g., `user_guide/concepts.md`) |
| `section` | text | Heading name or function name |
| `chunk_type` | text | `"markdown"`, `"docstring"`, or `"notebook"` |
| `metadata` | map | Additional key-value pairs (heading, cell_index, etc.) |
| `indexed_at` | timestamp | When the chunk was indexed |

Qdrant stores the same `id` as point ID, the embedding vector, and a lightweight payload (`source`, `section`, `chunk_type`) for pre-filtering.

### Public API

```python
from hera.utils.rag.indexer import build_index, update_file, delete_file
from pathlib import Path

# Full build
result = build_index(
    docs_root=Path("docs"),
    code_root=Path("hera"),
    clean=True,          # wipe and rebuild
    batch_size=64        # chunks per embedding batch
)
# {'docs': 350, 'docstrings': 800, 'notebooks': 50, 'total': 1200}

# Incremental update (single file)
update_file(Path("docs/user_guide/concepts.md"), root=Path("docs"))

# Delete chunks for a removed file
delete_file(Path("docs/old_page.md"), root=Path("docs"))
```

---

## Search (`search.py`)

### Retrieval pipeline

1. Embed the query using the same model as indexing
2. Search Qdrant for nearest vectors (with optional chunk_type filter)
3. Hydrate full text from Cassandra using the matched UUIDs
4. Optionally filter by `source_prefix`

### RAGSearch class

```python
from hera.utils.rag import RAGSearch

rag = RAGSearch()

# Full answer (retrieve + LLM)
answer = rag.ask("How does the risk assessment toolkit work?")

# Streaming answer (token by token)
for token in rag.stream("Show me the demography API"):
    print(token, end="", flush=True)

# Raw chunks only (no LLM)
chunks = rag.retrieve("data source versioning", top_k=10)
for c in chunks:
    print(f"{c['score']:.3f}  {c['source']} § {c['section']}")
```

### Filtering

```python
# Only search markdown docs (skip docstrings)
chunks = rag.retrieve("project lifecycle", chunk_type="markdown")

# Only search Python docstrings
chunks = rag.retrieve("addDataSource parameters", chunk_type="docstring")

# Only search files under a specific path
chunks = rag.retrieve("OpenFOAM templates", source_prefix="toolkits/simulations")
```

---

## REST API (`serve.py`)

### Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/health` | Health check (returns model name) |
| `POST` | `/ingest` | Build/rebuild index |
| `POST` | `/search` | Search + optional LLM answer |
| `GET` | `/search?q=...` | Same via query params |
| `GET` | `/stream?q=...` | SSE streaming answer |
| `DELETE` | `/source` | Delete chunks for a file |

### Starting the server

```bash
# Via CLI
hera-rag-search serve
hera-rag-search serve --port 9000 --reload

# Via Makefile
make rag-serve
make rag-serve-watch    # serve + auto re-index on file changes
```

### Example API calls

```bash
# Health check
curl http://localhost:8765/health

# Search
curl -X POST http://localhost:8765/search \
  -H "Content-Type: application/json" \
  -d '{"query": "How do I create a project?", "top_k": 5}'

# Stream answer (SSE)
curl -N "http://localhost:8765/stream?q=How+does+the+LSM+toolkit+work"

# Build index
curl -X POST http://localhost:8765/ingest \
  -H "Content-Type: application/json" \
  -d '{"docs_root": "docs", "code_root": "hera", "clean": true}'
```

---

## File watcher

Auto re-indexes files when they change on disk:

```bash
hera-rag-search watch --root docs
make rag-watch
```

- Watches for `.md`, `.py`, `.ipynb` file changes
- Debounced (default 2 seconds) to avoid re-indexing during rapid edits
- Deletes old chunks and inserts new ones for each changed file
- Can run alongside the API server: `make rag-serve-watch`

---

## MkDocs plugin

A floating "Ask AI" button can be injected into the MkDocs site:

1. Set `RAG_ENABLED=true` in `.env` or environment
2. Start the RAG API server: `make rag-serve`
3. Build or serve docs: `make rag-docs-serve`

The plugin injects a floating action button (bottom-right corner) that opens a search panel. Questions are answered via SSE streaming from the Ollama backend. Sources are displayed with relevance scores.

To register manually in `mkdocs.yml`:
```yaml
plugins:
  - hera_rag_search
```

---

## CLI reference

```bash
# Build index
hera-rag-search index                          # index docs/ + hera/
hera-rag-search index --clean                  # wipe and rebuild
hera-rag-search index --docs-only              # skip Python docstrings
hera-rag-search index --docs /other/docs       # custom docs path

# Search
hera-rag-search search "How do I use toolkits?"
hera-rag-search search "demography API" --raw  # chunks only, no LLM
hera-rag-search search "OpenFOAM" --type docstring  # filter by type
hera-rag-search search "LSM" --source toolkits/     # filter by path

# Watch
hera-rag-search watch                          # watch docs/ for changes
hera-rag-search watch --root hera              # watch source code

# Serve
hera-rag-search serve                          # start REST API
hera-rag-search serve --with-watcher           # serve + watch
hera-rag-search serve --port 9000 --reload     # dev mode
```

---

## Makefile targets

| Target | Description |
|--------|-------------|
| `make rag-index` | Build index from docs/ + hera/ |
| `make rag-reindex` | Wipe and rebuild |
| `make rag-index-docs` | Index docs only (skip docstrings) |
| `make rag-search` | Search with default query |
| `make rag-search-raw` | Search without LLM |
| `make rag-serve` | Start REST API server |
| `make rag-serve-watch` | Serve + auto re-index |
| `make rag-watch` | Watch docs/ for changes |
| `make rag-docs-serve` | Serve docs with RAG widget enabled |
| `make rag-docs-build` | Build docs with RAG widget enabled |

Override defaults: `make rag-search RAG_QUERY="my question"`

---

## Prerequisites

Install the RAG dependencies:
```bash
pip install -e .[rag]
```

Required services (can be Docker containers):
- **Qdrant**: `docker run -p 6333:6333 qdrant/qdrant`
- **Cassandra**: `docker run -p 9042:9042 cassandra:4`
- **Ollama**: [Install Ollama](https://ollama.ai/) then `ollama pull llama3`

---

## Extending

### Custom embedding model

```bash
RAG_EMBED_MODEL=intfloat/e5-large-v2 RAG_EMBED_DIM=1024 hera-rag-search index --clean
```

### Custom LLM

```bash
RAG_OLLAMA_MODEL=mistral hera-rag-search search "my question"
```

### Custom chunking

Modify `chunk_size`, `chunk_overlap`, and `min_chunk_chars` in `.env`:
```bash
RAG_CHUNK_SIZE=1500
RAG_CHUNK_OVERLAP=200
RAG_MIN_CHUNK_CHARS=50
```

### Adding new file types

Add a new chunker function in `indexer.py` following the pattern of `_chunk_markdown`, and register the extension in `SUPPORTED_EXTENSIONS` and the `_chunk_file` dispatch dict.
