# AI-Powered Documentation Search

Hera includes an optional AI search feature that lets you ask natural language questions about the documentation and source code. It uses a self-hosted LLM (via Ollama) to generate answers based on the indexed content.

---

## How it works

When enabled, a floating search button appears on every documentation page. Click it to open a search panel where you can type questions in plain English. The system:

1. Searches the indexed docs and docstrings for relevant content
2. Sends the top matches + your question to a local LLM
3. Streams the answer back in real-time
4. Shows source references with relevance scores

The entire pipeline runs locally — no data is sent to external services.

---

## Setup

### 1. Install RAG dependencies and start services

```bash
make install-rag
```

This runs the full setup:
- Installs Python dependencies (`pip install -e .[rag]`)
- Starts Qdrant (vector store), Cassandra (text store), and Ollama (LLM) via Docker
- Downloads the LLM model (default: `llama3`)
- Builds the search index from `docs/` and `hera/`

### 2. Start the RAG API server

```bash
make rag-serve
```

The API runs at `http://localhost:8765` by default.

### 3. Serve docs with AI search enabled

```bash
RAG_ENABLED=true mkdocs serve
```

Or use the Makefile shortcut:

```bash
make rag-docs-serve
```

This starts both the RAG API and the MkDocs dev server with the search widget injected.

---

## Using the search widget

1. Look for the search button in the bottom-right corner of any docs page
2. Click it to open the search panel
3. Type your question (e.g., "How do I create a project and load a repository?")
4. The answer streams in real-time from the local LLM
5. Source references are shown below the answer with relevance scores

---

## Configuration

All settings are controlled via environment variables (prefix `RAG_`):

| Variable | Default | Description |
|----------|---------|-------------|
| `RAG_ENABLED` | `false` | Set to `true` to enable the search widget |
| `RAG_API_URL` | `http://localhost:8765` | URL of the RAG API server |
| `RAG_OLLAMA_MODEL` | `llama3` | LLM model to use |
| `RAG_TOP_K` | `5` | Number of context chunks per query |

You can set these in a `.env` file (see `.env.example`) or export them:

```bash
export RAG_ENABLED=true
export RAG_OLLAMA_MODEL=mistral
```

---

## CLI search (without the web widget)

You can also search from the terminal:

```bash
# Search with LLM answer
hera-rag-search search "How does the risk assessment toolkit work?"

# Search without LLM (just show matching chunks)
hera-rag-search search "demography API" --raw

# Filter by content type
hera-rag-search search "addDataSource" --type docstring
hera-rag-search search "project lifecycle" --type markdown
```

---

## Rebuilding the index

After updating docs or source code, rebuild the index:

```bash
# Incremental rebuild
make rag-index

# Full rebuild (wipe and re-index everything)
make rag-reindex
```

Or use the file watcher to auto-rebuild on changes:

```bash
make rag-serve-watch
```

---

## Disabling AI search

Set `RAG_ENABLED=false` (the default) or simply don't start the RAG server. The docs site works normally without it — the widget only appears when both `RAG_ENABLED=true` and the RAG API is accessible.

---

## Requirements

- **Docker** — for Qdrant, Cassandra, and Ollama containers
- **~4GB disk** — for the LLM model (llama3)
- **~2GB RAM** — for the embedding model and services
- All traffic stays local — no external API calls
