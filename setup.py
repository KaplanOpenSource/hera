import glob
from setuptools import setup, find_packages

setup(
    name="pyhera",
    url="https://github.com/KaplanOpenSource/hera",
    packages=find_packages(),
    author="Yehuda Arav",
    maintainer="Yehuda Arav",
    maintainer_email="",
    description="Scientific data management and analysis platform",
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.9",
    scripts=[s for s in glob.glob("hera/bin/hera-*") if not s.endswith(".old")],
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