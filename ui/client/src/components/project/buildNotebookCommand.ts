import { DocumentDesc } from "../../shared/types";

export const buildNotebookCommand = ({
  notebookResource,
  projectName,
  toolkitNames,
  desc,
}: {
  notebookResource: string,
  projectName: string,
  toolkitNames: string[],
  desc: DocumentDesc,
}) => {
  const toolkitSourceLines = toolkitNames.map(name => {
    const varName = name.replace(/[^a-zA-Z0-9_]/g, '_');
    return `                "${varName}_tk = hera.toolkit.ToolkitHome().getToolkit(\\"${name}\\", PROJECT_NAME)\\n",`;
  }).join('\n');
  const toolkitSection = toolkitNames.length > 0
    ? `\n                "\\n",\n${toolkitSourceLines}`
    : '';
  return `
import json
from pathlib import Path
notebook_path = Path("${notebookResource}")
if not notebook_path.exists():
    notebook_path.parent.mkdir(parents=True, exist_ok=True)
    empty_notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {"kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        }},
        "cells": [{
            "cell_type": "code",
            "metadata": {},
            "source": [
                "import hera\\n",
                "PROJECT_NAME = \\"${projectName}\\"\\n",
                "\\n",
                "p = hera.datalayer.Project(projectName=PROJECT_NAME)\\n",
                ${toolkitSection}
                "p.getMetadata()\\n",
            ],
            "outputs": [],
            "execution_count": None,
        }],
    }
    notebook_path.write_text(json.dumps(empty_notebook, indent=2))
Cache_Collection().addDocument(
    projectName="${projectName}",
    resource="${notebookResource}",
    dataFormat="JSON_dict",
    type="notebook",
    desc=${JSON.stringify(desc)},
)`;
};
