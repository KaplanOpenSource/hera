"""
RAG configuration — all values overridable via environment variables or .env.
"""
from __future__ import annotations
from typing import List
from pydantic import Field
from pydantic_settings import BaseSettings


class RagSettings(BaseSettings):
    # ── Cassandra ────────────────────────────────────────────────────────────
    cassandra_hosts: List[str] = Field(default=["localhost"])
    cassandra_port: int = 9042
    cassandra_keyspace: str = "hera_rag"

    # ── Qdrant ───────────────────────────────────────────────────────────────
    qdrant_host: str = "localhost"
    qdrant_port: int = 6333
    qdrant_collection: str = "hera_docs"

    # ── Ollama ───────────────────────────────────────────────────────────────
    ollama_base_url: str = "http://localhost:11434"
    ollama_model: str = "llama3"
    ollama_timeout: int = 120

    # ── Embedding ────────────────────────────────────────────────────────────
    embed_model: str = "BAAI/bge-small-en-v1.5"
    embed_dim: int = 384

    # ── Offline / local model cache ──────────────────────────────────────────
    # When offline_mode=True, no network fetches are attempted: embedding
    # models must already exist under models_cache_dir, and the Ollama server
    # must already be reachable at ollama_base_url (typically a native install).
    offline_mode: bool = False
    models_cache_dir: str = "~/.pyhera/rag_models"

    # ── Chunking ─────────────────────────────────────────────────────────────
    chunk_size: int = 1_000
    chunk_overlap: int = 150
    min_chunk_chars: int = 80

    # ── Search ───────────────────────────────────────────────────────────────
    top_k: int = 5
    rag_prompt_template: str = (
        "You are a helpful assistant for the Hera codebase. "
        "Answer the question based ONLY on the documentation and code excerpts below. "
        "If the answer is not in the excerpts, say so.\n\n"
        "Context:\n{context}\n\n"
        "Question: {question}\n\nAnswer:"
    )

    # ── API / MkDocs ─────────────────────────────────────────────────────────
    api_host: str = "0.0.0.0"
    api_port: int = 8765
    rag_enabled: bool = False
    rag_api_url: str = "http://localhost:8765"

    model_config = {"env_file": ".env", "env_prefix": "RAG_", "env_file_encoding": "utf-8"}


settings = RagSettings()
