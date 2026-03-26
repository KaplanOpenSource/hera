"""
hera/utils/rag/search.py

Query interface for the RAG index.

    from hera.utils.rag import RAGSearch

    rag = RAGSearch()
    answer = rag.ask("How does the risk assessment toolkit calculate casualties?")
    for token in rag.stream("What's the demography layer API?"):
        print(token, end="", flush=True)

Also exposes lower-level retrieve() if you just want the chunks.
"""
from __future__ import annotations
import uuid
from typing import Generator, Optional

import httpx
from qdrant_client.models import FieldCondition, Filter, MatchValue

from .config import settings
from .indexer import _cass, _qdrant, _embed


def retrieve(
    query: str,
    top_k: int | None = None,
    chunk_type: str | None = None,
    source_prefix: str | None = None,
) -> list[dict]:
    """
    Embed query, search Qdrant ANN, hydrate from Cassandra.

    Parameters
    ----------
    query : str
        Natural language question.
    top_k : int, optional
        Number of chunks (default from settings).
    chunk_type : str, optional
        Filter to "markdown", "docstring", or "notebook".
    source_prefix : str, optional
        Filter to sources starting with this path fragment.

    Returns
    -------
    list of dict
        Each dict: {text, source, section, chunk_type, metadata, score}
    """
    k = top_k or settings.top_k
    q_vec = _embed([query])[0]

    conditions = []
    if chunk_type:
        conditions.append(FieldCondition(key="chunk_type", match=MatchValue(value=chunk_type)))
    f = Filter(must=conditions) if conditions else None

    hits = _qdrant().search(
        settings.qdrant_collection,
        query_vector=q_vec,
        limit=k * 3 if source_prefix else k,
        query_filter=f,
    )

    results = []
    select = _cass().prepare(
        f"SELECT * FROM {settings.cassandra_keyspace}.chunks WHERE id = ?"
    )
    for hit in hits:
        if source_prefix and not hit.payload.get("source", "").startswith(source_prefix):
            continue
        row = _cass().execute(select, (uuid.UUID(str(hit.id)),)).one()
        if row:
            results.append({
                "text":       row.text,
                "source":     row.source,
                "section":    row.section,
                "chunk_type": row.chunk_type,
                "metadata":   dict(row.metadata or {}),
                "score":      round(hit.score, 4),
            })
        if len(results) >= k:
            break

    return results


def _build_prompt(query: str, chunks: list[dict]) -> str:
    """Build the LLM prompt from query and retrieved chunks."""
    context = "\n\n---\n\n".join(
        f"[{c['source']} § {c['section']}]\n{c['text']}"
        for c in chunks
    )
    return settings.rag_prompt_template.format(context=context, question=query)


def _ollama_generate(prompt: str, stream: bool = False):
    """Call Ollama API to generate a response."""
    with httpx.Client(timeout=settings.ollama_timeout) as client:
        if not stream:
            resp = client.post(
                f"{settings.ollama_base_url}/api/generate",
                json={"model": settings.ollama_model, "prompt": prompt, "stream": False},
            )
            resp.raise_for_status()
            return resp.json()["response"].strip()
        else:
            import json as _json
            with client.stream(
                "POST",
                f"{settings.ollama_base_url}/api/generate",
                json={"model": settings.ollama_model, "prompt": prompt, "stream": True},
            ) as resp:
                resp.raise_for_status()
                for line in resp.iter_lines():
                    if line:
                        data = _json.loads(line)
                        yield data.get("response", "")
                        if data.get("done"):
                            break


class RAGSearch:
    """
    High-level search interface.

    Example
    -------
    >>> from hera.utils.rag import RAGSearch
    >>> rag = RAGSearch()
    >>> answer = rag.ask("How does X work?")
    >>> for tok in rag.stream("Show me the Y API"):
    ...     print(tok, end="")
    >>> chunks = rag.retrieve("raw chunks only")
    """

    def retrieve(self, query, top_k=None, chunk_type=None, source_prefix=None) -> list[dict]:
        """Return ranked chunks without calling the LLM."""
        return retrieve(query, top_k=top_k, chunk_type=chunk_type, source_prefix=source_prefix)

    def ask(self, query, top_k=None, chunk_type=None) -> str:
        """Retrieve chunks and return a complete answer string from Ollama."""
        chunks = retrieve(query, top_k=top_k, chunk_type=chunk_type)
        if not chunks:
            return "No relevant documentation found."
        return _ollama_generate(_build_prompt(query, chunks), stream=False)

    def stream(self, query, top_k=None, chunk_type=None) -> Generator[str, None, None]:
        """Retrieve chunks and stream answer tokens from Ollama."""
        chunks = retrieve(query, top_k=top_k, chunk_type=chunk_type)
        if not chunks:
            yield "No relevant documentation found."
            return
        yield from _ollama_generate(_build_prompt(query, chunks), stream=True)
