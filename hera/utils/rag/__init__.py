"""
hera.utils.rag — RAG search over hera docs and source code.

Quick start::

    from hera.utils.rag import RAGSearch
    rag = RAGSearch()
    print(rag.ask("How do I create a project?"))
"""
from .search import RAGSearch, retrieve
from .indexer import build_index, update_file, delete_file
from .config import settings as rag_settings

__all__ = [
    "RAGSearch",
    "retrieve",
    "build_index",
    "update_file",
    "delete_file",
    "rag_settings",
]
