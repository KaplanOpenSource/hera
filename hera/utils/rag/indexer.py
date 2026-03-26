"""
hera/utils/rag/indexer.py

Builds (or rebuilds) the RAG index from:
  - docs/        Markdown files, split by heading
  - hera/        Python docstrings, extracted via AST
  - docs/**      Jupyter notebooks, extracted cell-by-cell

Each chunk is dual-written:
  Cassandra  →  full text + metadata  (source of truth)
  Qdrant     →  embedding vector      (ANN retrieval index)

The two stores are linked by a shared UUID (chunk_id).

Public API:
    build_index(docs_root, code_root)   full build from scratch
    update_file(path, docs_root)        re-index a single file (used by watcher)
    delete_file(path, docs_root)        remove chunks for a single file
"""
from __future__ import annotations

import ast
import json
import re
import textwrap
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

from cassandra.cluster import Cluster, Session
from cassandra.policies import DCAwareRoundRobinPolicy
from qdrant_client import QdrantClient
from qdrant_client.models import (
    Distance, FieldCondition, Filter,
    MatchValue, PointIdsList, PointStruct, VectorParams,
)
from sentence_transformers import SentenceTransformer

from .config import settings

# ── Supported extensions ──────────────────────────────────────────────────────

SUPPORTED_EXTENSIONS = {".md", ".mdx", ".py", ".ipynb"}

# ── Cassandra DDL ─────────────────────────────────────────────────────────────

_KEYSPACE_DDL = (
    f"CREATE KEYSPACE IF NOT EXISTS {settings.cassandra_keyspace} "
    "WITH replication = {'class': 'SimpleStrategy', 'replication_factor': 1}"
)
_TABLE_DDL = f"""
CREATE TABLE IF NOT EXISTS {settings.cassandra_keyspace}.chunks (
    id         uuid PRIMARY KEY,
    text       text,
    source     text,
    section    text,
    chunk_type text,
    metadata   map<text, text>,
    indexed_at timestamp
)
"""
_SOURCE_INDEX_DDL = (
    f"CREATE INDEX IF NOT EXISTS ON {settings.cassandra_keyspace}.chunks (source)"
)


# ═════════════════════════════════════════════════════════════════════════════
# Backend connections (module-level singletons, lazy)
# ═════════════════════════════════════════════════════════════════════════════

_cass_session: Session | None = None
_qdrant_client: QdrantClient | None = None
_embedder: SentenceTransformer | None = None


def _cass() -> Session:
    """Get or create the Cassandra session singleton."""
    global _cass_session
    if _cass_session is None:
        cluster = Cluster(
            settings.cassandra_hosts,
            port=settings.cassandra_port,
            load_balancing_policy=DCAwareRoundRobinPolicy(),
            protocol_version=4,
        )
        s = cluster.connect()
        s.execute(_KEYSPACE_DDL)
        s.execute(_TABLE_DDL)
        s.execute(_SOURCE_INDEX_DDL)
        s.set_keyspace(settings.cassandra_keyspace)
        _cass_session = s
    return _cass_session


def _qdrant() -> QdrantClient:
    """Get or create the Qdrant client singleton."""
    global _qdrant_client
    if _qdrant_client is None:
        client = QdrantClient(host=settings.qdrant_host, port=settings.qdrant_port)
        existing = {c.name for c in client.get_collections().collections}
        if settings.qdrant_collection not in existing:
            client.create_collection(
                settings.qdrant_collection,
                vectors_config=VectorParams(
                    size=settings.embed_dim, distance=Distance.COSINE
                ),
            )
        _qdrant_client = client
    return _qdrant_client


def _embed_model() -> SentenceTransformer:
    """Get or create the embedding model singleton."""
    global _embedder
    if _embedder is None:
        _embedder = SentenceTransformer(settings.embed_model)
    return _embedder


def _embed(texts: list[str]) -> list[list[float]]:
    """Embed a list of texts and return normalized vectors."""
    return _embed_model().encode(
        texts, show_progress_bar=False, normalize_embeddings=True
    ).tolist()


# ═════════════════════════════════════════════════════════════════════════════
# Chunkers
# ═════════════════════════════════════════════════════════════════════════════

_HEADING_RE = re.compile(r"^(#{1,3})\s+(.+)", re.MULTILINE)


def _sliding_split(text: str) -> list[str]:
    """Sliding window split for bodies that exceed chunk_size."""
    size, overlap = settings.chunk_size, settings.chunk_overlap
    if len(text) <= size:
        return [text]
    parts, start = [], 0
    while start < len(text):
        parts.append(text[start : start + size])
        start += size - overlap
    return parts


def _chunk_markdown(path: Path, root: Path) -> list[dict]:
    """Split a markdown file into chunks by heading."""
    text = path.read_text(errors="replace")
    rel = str(path.relative_to(root))
    matches = list(_HEADING_RE.finditer(text))
    sections: list[tuple[str, str]] = []

    if not matches:
        sections = [("(root)", text)]
    else:
        for i, m in enumerate(matches):
            heading = m.group(2).strip()
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
            sections.append((heading, text[start:end].strip()))

    chunks = []
    for heading, body in sections:
        for idx, part in enumerate(_sliding_split(body)):
            if len(part) < settings.min_chunk_chars:
                continue
            label = heading if idx == 0 else f"{heading} (cont. {idx})"
            chunks.append({
                "text": part, "source": rel, "section": label,
                "chunk_type": "markdown",
                "metadata": {"heading": heading, "part": str(idx)},
            })
    return chunks


def _extract_docstrings(tree: ast.AST) -> Iterator[tuple[str, str]]:
    """Walk an AST and yield (name, docstring) pairs."""
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        doc = ast.get_docstring(node)
        if doc:
            yield node.name, doc


def _chunk_python(path: Path, root: Path) -> list[dict]:
    """Extract docstrings from a Python file as chunks."""
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return []

    rel = str(path.relative_to(root))
    chunks = []

    module_doc = ast.get_docstring(tree)
    if module_doc and len(module_doc) >= settings.min_chunk_chars:
        chunks.append({
            "text": module_doc, "source": rel, "section": "(module)",
            "chunk_type": "docstring", "metadata": {"name": "(module)"},
        })

    for name, doc in _extract_docstrings(tree):
        body = f"{name}:\n{textwrap.dedent(doc)}"
        if len(body) < settings.min_chunk_chars:
            continue
        chunks.append({
            "text": body, "source": rel, "section": name,
            "chunk_type": "docstring", "metadata": {"name": name},
        })
    return chunks


def _chunk_notebook(path: Path, root: Path) -> list[dict]:
    """Extract cells from a Jupyter notebook as chunks."""
    try:
        nb = json.loads(path.read_text(errors="replace"))
    except json.JSONDecodeError:
        return []

    rel = str(path.relative_to(root))
    chunks = []
    for idx, cell in enumerate(nb.get("cells", [])):
        src = "".join(cell.get("source", []))
        cell_type = cell.get("cell_type", "")
        if len(src) < settings.min_chunk_chars:
            continue
        chunks.append({
            "text": src, "source": rel,
            "section": f"cell_{idx}_{cell_type}",
            "chunk_type": "notebook",
            "metadata": {"cell_index": str(idx), "cell_type": cell_type},
        })
    return chunks


def _chunk_file(path: Path, root: Path) -> list[dict]:
    """Dispatch to the appropriate chunker based on file extension."""
    dispatch = {
        ".md": _chunk_markdown, ".mdx": _chunk_markdown,
        ".py": _chunk_python,
        ".ipynb": _chunk_notebook,
    }
    fn = dispatch.get(path.suffix.lower())
    return fn(path, root) if fn else []


# ═════════════════════════════════════════════════════════════════════════════
# Store helpers
# ═════════════════════════════════════════════════════════════════════════════

_INSERT = (
    f"INSERT INTO {settings.cassandra_keyspace}.chunks "
    "(id, text, source, section, chunk_type, metadata, indexed_at) "
    "VALUES (%s, %s, %s, %s, %s, %s, %s)"
)

_DELETE = (
    f"DELETE FROM {settings.cassandra_keyspace}.chunks WHERE id = %s"
)


def _store_chunks(chunks: list[dict]) -> int:
    """Embed chunks and write to both Cassandra and Qdrant."""
    if not chunks:
        return 0
    vectors = _embed([c["text"] for c in chunks])
    now = datetime.now(timezone.utc)
    points: list[PointStruct] = []

    for chunk, vec in zip(chunks, vectors):
        cid = uuid.uuid4()
        _cass().execute(_INSERT, (
            cid, chunk["text"], chunk["source"], chunk["section"],
            chunk["chunk_type"], chunk.get("metadata", {}), now,
        ))
        points.append(PointStruct(
            id=str(cid), vector=vec,
            payload={
                "source": chunk["source"],
                "section": chunk["section"],
                "chunk_type": chunk["chunk_type"],
            },
        ))

    _qdrant().upsert(settings.qdrant_collection, points=points)
    return len(chunks)


def _delete_source(rel_source: str) -> int:
    """Remove all chunks for a given relative source path from both stores."""
    scroll, offset = _qdrant().scroll(
        settings.qdrant_collection,
        scroll_filter=Filter(
            must=[FieldCondition(key="source", match=MatchValue(value=rel_source))]
        ),
        limit=10_000, with_payload=False,
    )
    ids = [str(p.id) for p in scroll]
    if not ids:
        return 0
    _qdrant().delete(
        settings.qdrant_collection,
        points_selector=PointIdsList(points=ids),
    )
    for cid in ids:
        _cass().execute(_DELETE, (uuid.UUID(cid),))
    return len(ids)


# ═════════════════════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════════════════════

def build_index(
    docs_root: Path,
    code_root: Path | None = None,
    *,
    clean: bool = False,
    batch_size: int = 64,
    progress_cb=None,
) -> dict:
    """
    Build (or rebuild) the full index.

    Parameters
    ----------
    docs_root : Path
        Path to docs/ directory (markdown, notebooks).
    code_root : Path or None
        Path to hera/ source directory (Python docstrings).
        Pass None to skip docstring indexing.
    clean : bool
        If True, wipe the entire collection before indexing.
    batch_size : int
        Chunks per embedding batch.
    progress_cb : callable, optional
        Optional callable(current, total, label) for progress reporting.

    Returns
    -------
    dict
        {"docs": int, "docstrings": int, "notebooks": int, "total": int}
    """
    if clean:
        _qdrant().delete_collection(settings.qdrant_collection)
        _cass().execute(
            f"TRUNCATE {settings.cassandra_keyspace}.chunks"
        )
        _qdrant().create_collection(
            settings.qdrant_collection,
            vectors_config=VectorParams(size=settings.embed_dim, distance=Distance.COSINE),
        )

    all_chunks: list[dict] = []

    for ext in (".md", ".mdx", ".ipynb"):
        for path in sorted(docs_root.rglob(f"*{ext}")):
            all_chunks.extend(_chunk_file(path, docs_root))

    if code_root and code_root.exists():
        for path in sorted(code_root.rglob("*.py")):
            all_chunks.extend(_chunk_python(path, code_root))

    counts = {"markdown": 0, "docstring": 0, "notebook": 0}
    total = len(all_chunks)

    for start in range(0, total, batch_size):
        batch = all_chunks[start : start + batch_size]
        _store_chunks(batch)
        for c in batch:
            counts[c["chunk_type"]] = counts.get(c["chunk_type"], 0) + 1
        if progress_cb:
            progress_cb(min(start + batch_size, total), total, "Indexing")

    return {
        "docs": counts.get("markdown", 0),
        "docstrings": counts.get("docstring", 0),
        "notebooks": counts.get("notebook", 0),
        "total": total,
    }


def update_file(path: Path, root: Path) -> dict:
    """
    Re-index a single file (delete old chunks, insert new ones).

    Returns
    -------
    dict
        {"source": str, "chunks": int}
    """
    rel = str(path.relative_to(root))
    _delete_source(rel)
    chunks = _chunk_file(path, root)
    n = _store_chunks(chunks)
    return {"source": rel, "chunks": n}


def delete_file(path: Path, root: Path) -> dict:
    """
    Remove all chunks for a file.

    Returns
    -------
    dict
        {"source": str, "deleted": int}
    """
    rel = str(path.relative_to(root))
    n = _delete_source(rel)
    return {"source": rel, "deleted": n}
