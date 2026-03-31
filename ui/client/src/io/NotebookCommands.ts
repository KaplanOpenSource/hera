import type { PythonCommand } from './fetchPython';

export const NotebookCommands = {
  list: (rootDir: string): PythonCommand => ({
    results: ['notebooks'],
    code: `
from pathlib import Path
notebooks_dir = Path("${rootDir}") / "notebooks"
if notebooks_dir.exists():
    notebooks = sorted(f.stem for f in notebooks_dir.iterdir() if f.suffix == ".ipynb" and f.is_file())
else:
    notebooks = []
`,
  }),

  create: (rootDir: string, name: string): PythonCommand => ({
    results: ['name'],
    code: `
import json
from pathlib import Path
name = "${name}"
notebooks_dir = Path("${rootDir}") / "notebooks"
notebooks_dir.mkdir(parents=True, exist_ok=True)
empty_notebook = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}},
    "cells": [],
}
(notebooks_dir / f"{name}.ipynb").write_text(json.dumps(empty_notebook, indent=2))
`,
  }),

  delete: (rootDir: string, name: string): PythonCommand => ({
    results: [],
    code: `
from pathlib import Path
(Path("${rootDir}") / "notebooks" / "${name}.ipynb").unlink()
`,
  }),
};
