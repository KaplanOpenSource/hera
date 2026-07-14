import glob
import re
from pathlib import Path
from setuptools import setup, find_packages

# Read version from the package without importing it
_version_match = re.search(
    r"^__version__\s*=\s*['\"]([^'\"]+)['\"]",
    Path("hera/__init__.py").read_text(encoding="utf-8"),
    re.M,
)
_VERSION = _version_match.group(1) if _version_match else "0.0.0"

setup(
    name="pyhera",
    version=_VERSION,
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
    install_requires=[
        "pandas>=1.3",
        "numpy>=1.21",
        "mongoengine>=0.24",
        "pymongo>=3.12",
        "pint>=0.19",
        "deprecated>=1.2",
        "scipy>=1.7",
        "xarray>=0.20",
        "dask>=2021.10",
        "geopandas>=0.10",
        "shapely>=1.8",
    ],
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