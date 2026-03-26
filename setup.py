from setuptools import setup, find_packages

setup(
    name="pyhera",
    url="http://mathsrv2:8081/edenn/pyhera.git",
    packages=find_packages(),
    author="Eden Nitzan",
    extras_require={
        "rag": [
            "sentence-transformers>=2.7",
            "cassandra-driver>=3.29",
            "qdrant-client>=1.9",
            "httpx>=0.27",
            "fastapi>=0.111",
            "uvicorn[standard]>=0.29",
            "typer[all]>=0.12",
            "rich>=13",
            "watchdog>=4.0",
            "pydantic-settings>=2.3",
        ],
    },
    entry_points={
        "mkdocs.plugins": [
            "hera_rag_search = hera.utils.rag.serve:HeraMkDocsPlugin",
        ],
    },
)