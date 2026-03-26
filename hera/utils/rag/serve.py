"""
hera/utils/rag/serve.py

Optional serving layer — FastAPI REST API, SSE streaming, file watcher,
and MkDocs plugin registration.

Start with:
    hera-rag-search serve
    python -m hera.utils.rag.serve
"""
from __future__ import annotations
import threading
import time
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from .config import settings
from .search import RAGSearch, retrieve
from .indexer import build_index, update_file, delete_file, SUPPORTED_EXTENSIONS

app = FastAPI(
    title="Hera RAG API",
    description="RAG search over hera docs and source — Ollama + Cassandra + Qdrant.",
    version="1.0.0",
)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

_rag: Optional[RAGSearch] = None

def _get_rag() -> RAGSearch:
    """Get or create the RAGSearch singleton."""
    global _rag
    if _rag is None:
        _rag = RAGSearch()
    return _rag


class IngestRequest(BaseModel):
    """Request body for the /ingest endpoint."""
    docs_root: str
    code_root: Optional[str] = None
    clean: bool = False

class SearchRequest(BaseModel):
    """Request body for the /search POST endpoint."""
    query: str
    top_k: int = settings.top_k
    no_llm: bool = False
    chunk_type: Optional[str] = None

class DeleteRequest(BaseModel):
    """Request body for the /source DELETE endpoint."""
    source: str


@app.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "model": settings.ollama_model}


@app.post("/ingest")
def ingest(req: IngestRequest):
    """Build or rebuild the RAG index."""
    docs = Path(req.docs_root).resolve()
    code = Path(req.code_root).resolve() if req.code_root else None
    if not docs.exists():
        raise HTTPException(400, f"docs_root does not exist: {docs}")
    return build_index(docs, code, clean=req.clean)


@app.post("/search")
def search_post(req: SearchRequest):
    """Search the index and optionally generate an LLM answer."""
    chunks = retrieve(req.query, top_k=req.top_k, chunk_type=req.chunk_type)
    if not chunks:
        return {"answer": None, "chunks": [], "query": req.query}
    answer = None if req.no_llm else _get_rag().ask(req.query, req.top_k, req.chunk_type)
    return {
        "query": req.query,
        "answer": answer,
        "chunks": [
            {k: v for k, v in c.items() if k != "text"} | {"excerpt": c["text"][:300]}
            for c in chunks
        ],
    }


@app.get("/search")
def search_get(
    q: str = Query(...),
    k: int = Query(settings.top_k),
    no_llm: bool = Query(False),
    chunk_type: Optional[str] = Query(None),
):
    """Search via GET parameters."""
    return search_post(SearchRequest(query=q, top_k=k, no_llm=no_llm, chunk_type=chunk_type))


@app.get("/stream")
def search_stream(q: str = Query(...), k: int = Query(settings.top_k), chunk_type: Optional[str] = Query(None)):
    """SSE endpoint — streams Ollama tokens as `data: <token>\\n\\n`."""
    def _gen():
        try:
            for token in _get_rag().stream(q, top_k=k, chunk_type=chunk_type):
                yield f"data: {token}\n\n"
            yield "data: [DONE]\n\n"
        except Exception as e:
            yield f"data: [ERROR: {e}]\n\n"
    return StreamingResponse(_gen(), media_type="text/event-stream")


@app.delete("/source")
def delete_source(req: DeleteRequest):
    """Delete all chunks for a given source file."""
    from .indexer import _delete_source
    return {"deleted_chunks": _delete_source(req.source), "source": req.source}


def run_server(host=None, port=None, reload=False):
    """Start the FastAPI server with uvicorn."""
    import uvicorn
    uvicorn.run("hera.utils.rag.serve:app", host=host or settings.api_host,
                port=port or settings.api_port, reload=reload)


# ── File watcher ──────────────────────────────────────────────────────────────

class _RagWatchHandler:
    """Debounced file change handler for auto re-indexing."""

    def __init__(self, root: Path, debounce: float = 2.0):
        self._root = root
        self._debounce = debounce
        self._timers: dict[str, threading.Timer] = {}
        self._lock = threading.Lock()

    def dispatch(self, event):
        """Handle file system events."""
        from watchdog.events import FileDeletedEvent
        path = Path(getattr(event, "dest_path", event.src_path))
        if event.is_directory or path.suffix.lower() not in SUPPORTED_EXTENSIONS:
            return
        is_delete = isinstance(event, FileDeletedEvent)
        key = str(path)
        with self._lock:
            if key in self._timers:
                self._timers[key].cancel()
            fn = self._on_delete if is_delete else self._on_upsert
            t = threading.Timer(self._debounce, fn, args=[path])
            self._timers[key] = t
            t.start()

    def _on_upsert(self, path: Path):
        """Re-index a changed file."""
        try:
            r = update_file(path, self._root)
            print(f"[rag:watcher] re-indexed {r['source']} -> {r['chunks']} chunks")
        except Exception as e:
            print(f"[rag:watcher] error {path}: {e}")
        finally:
            self._timers.pop(str(path), None)

    def _on_delete(self, path: Path):
        """Remove chunks for a deleted file."""
        try:
            r = delete_file(path, self._root)
            print(f"[rag:watcher] deleted {r['deleted']} chunks for {r['source']}")
        except Exception as e:
            print(f"[rag:watcher] error {path}: {e}")
        finally:
            self._timers.pop(str(path), None)


def run_watcher(root: Path, debounce: float = 2.0):
    """Start a file watcher that auto re-indexes on changes."""
    from watchdog.observers import Observer
    handler = _RagWatchHandler(root, debounce=debounce)
    observer = Observer()
    observer.schedule(handler, str(root), recursive=True)
    observer.start()
    print(f"[rag:watcher] Watching {root} (debounce={debounce}s) — Ctrl-C to stop")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()
    print("[rag:watcher] Stopped.")


# ── MkDocs plugin ─────────────────────────────────────────────────────────────

_WIDGET_HTML = """
<style>
  #rag-fab {{
    position:fixed;bottom:24px;right:24px;z-index:9999;
    width:52px;height:52px;border-radius:50%;
    background:var(--md-primary-fg-color,#1976d2);color:#fff;
    border:none;cursor:pointer;font-size:22px;
    box-shadow:0 4px 12px rgba(0,0,0,.3);
    display:flex;align-items:center;justify-content:center;
    transition:transform .15s;
  }}
  #rag-fab:hover{{transform:scale(1.1)}}
  #rag-panel {{
    display:none;position:fixed;bottom:84px;right:24px;
    width:380px;max-height:520px;overflow-y:auto;
    background:var(--md-default-bg-color,#fff);
    border:1px solid var(--md-default-fg-color--lighter,#ddd);
    border-radius:12px;padding:16px;
    box-shadow:0 8px 32px rgba(0,0,0,.18);z-index:9998;
    font-family:var(--md-text-font,sans-serif);
  }}
  #rag-panel.open{{display:block}}
  #rag-input{{width:100%;padding:8px 10px;border-radius:8px;
    border:1px solid #ccc;font-size:14px;
    background:var(--md-default-bg-color,#fff);
    color:var(--md-default-fg-color,#000);}}
  #rag-btn{{margin-top:8px;width:100%;padding:8px;
    background:var(--md-primary-fg-color,#1976d2);
    color:#fff;border:none;border-radius:8px;cursor:pointer;font-size:14px;}}
  #rag-answer{{margin-top:12px;font-size:13px;line-height:1.6;white-space:pre-wrap;}}
  #rag-sources{{margin-top:10px;font-size:11px;color:#666;}}
  #rag-sources a{{color:var(--md-primary-fg-color,#1976d2);}}
</style>
<button id="rag-fab" title="Ask AI" aria-label="Ask the Hera docs">&#128269;</button>
<div id="rag-panel">
  <strong style="font-size:14px">Ask the Hera docs</strong>
  <input id="rag-input" type="text" placeholder="Ask a question..."/>
  <button id="rag-btn">Ask</button>
  <div id="rag-answer"></div>
  <div id="rag-sources"></div>
</div>
<script>
(function(){{
  const API="{api_url}";
  const fab=document.getElementById("rag-fab");
  const panel=document.getElementById("rag-panel");
  const input=document.getElementById("rag-input");
  const btn=document.getElementById("rag-btn");
  const answerEl=document.getElementById("rag-answer");
  const sourcesEl=document.getElementById("rag-sources");
  fab.addEventListener("click",()=>panel.classList.toggle("open"));
  async function ask(){{
    const q=input.value.trim();if(!q)return;
    answerEl.textContent="Thinking...";sourcesEl.innerHTML="";btn.disabled=true;
    let text="";
    const es=new EventSource(API+"/stream?q="+encodeURIComponent(q));
    es.onmessage=(e)=>{{
      if(e.data==="[DONE]"){{es.close();btn.disabled=false;return;}}
      if(e.data.startsWith("[")){{answerEl.textContent=e.data;es.close();btn.disabled=false;return;}}
      text+=e.data;answerEl.textContent=text;
    }};
    es.onerror=()=>{{answerEl.textContent="Error connecting to RAG server.";es.close();btn.disabled=false;}};
    try{{
      const r=await fetch(API+"/search?q="+encodeURIComponent(q)+"&no_llm=true");
      const d=await r.json();
      if(d.chunks&&d.chunks.length){{
        sourcesEl.innerHTML="<b>Sources:</b><br>"+
          d.chunks.map(c=>"<a href='#'>"+c.source+" - "+c.section+"</a> ("+(c.score*100).toFixed(1)+"%)").join("<br>");
      }}
    }}catch(e){{}}
  }}
  btn.addEventListener("click",ask);
  input.addEventListener("keydown",(e)=>{{if(e.key==="Enter")ask();}});
}})();
</script>
"""

try:
    from mkdocs.plugins import BasePlugin
    from mkdocs.config.defaults import MkDocsConfig

    class HeraMkDocsPlugin(BasePlugin):
        """
        MkDocs plugin — injects floating RAG search widget.
        Activated only when RAG_ENABLED=true.
        Register in mkdocs.yml as:  plugins: [hera_rag_search]
        """
        config_scheme = ()

        def on_page_content(self, html: str, *, page, config: MkDocsConfig, files) -> str:
            """Inject the RAG widget HTML before </body>."""
            if not settings.rag_enabled:
                return html
            widget = _WIDGET_HTML.format(api_url=settings.rag_api_url.rstrip("/"))
            return html.replace("</body>", widget + "\n</body>") if "</body>" in html else html + widget

except ImportError:
    pass


if __name__ == "__main__":
    run_server()
