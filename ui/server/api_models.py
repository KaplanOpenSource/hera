"""Request/response models for the Hera UI API endpoints."""

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel


class JupyterStartPayload(BaseModel):
    root_dir: str
    dark: bool = False


class ExecPayload(BaseModel):
    code: str


class Problem(BaseModel):
    error: str
    traceback: str


class ExecResponse(BaseModel):
    data: Any = None
    problem: Optional[Problem] = None


class RunWorkflowPayload(BaseModel):
    projectName: str
    # The whole workflow document (desc.workflow, desc.workflowName, resource, ...).
    # The client always sends it, so the run builds straight from it, no DB lookup.
    doc: Dict[str, Any]


class WorkflowChunk(BaseModel):
    # One output segment: a task's name (or "__preamble__" / "__between__") and the
    # console output captured while that segment was current.
    name: str
    text: str


class RunWorkflowResponse(BaseModel):
    # start returns token (or status "busy"); poll returns status + chunks/error.
    token: Optional[str] = None
    status: Optional[str] = None
    error: str = ""
    # Per-task output segments, in run order. Grows live while the run is going;
    # None only before any output (or when the token is unknown).
    chunks: Optional[List[WorkflowChunk]] = None


class ValidateNodeParamsPayload(BaseModel):
    # The node's dotted Hermes type, e.g. "RiskAssessment.calculateThresholds".
    type: str
    # The node's input_parameters, exactly as the editor holds them.
    params: Dict[str, Any] = {}


class ValidateNodeParamsResponse(BaseModel):
    # ok with an empty message also means "nothing to say" (unknown type, or a
    # node with no checks of its own).
    ok: bool = True
    message: str = ""


class CreateProjectPayload(BaseModel):
    projectName: str
    # Empty means the default: <cwd>/projects/<projectName>.
    filesDirectory: str = ""
    loadRepositories: bool = True


class DeleteProjectPayload(BaseModel):
    projectName: str
    # Also remove the files under the project's own files directory.
    deleteFiles: bool = False


class DetectNotebooksPayload(BaseModel):
    projectName: str
    filesDirectory: str


class Collection(str, Enum):
    """The three collections a document can live in."""
    Measurements = "Measurements"
    Simulations = "Simulations"
    Cache = "Cache"


class NotebookSpec(BaseModel):
    # The toolkit lines the notebook's first cell opens with.
    toolkitNames: List[str] = []


class CreateDocumentPayload(BaseModel):
    projectName: str
    collection: Collection
    # The defaults are hera's own addDocument defaults.
    type: str = ""
    # A path, or the data itself (an agent's resource is a dict).
    resource: Any = ""
    dataFormat: str = "string"
    desc: Dict[str, Any] = {}
    # Only for a notebook document: write the .ipynb file if it is not there yet.
    notebook: Optional[NotebookSpec] = None


class UpdateDocumentPayload(BaseModel):
    docId: str
    # The document as the editor holds it, minus _id, _cls and projectName.
    fields: Dict[str, Any] = {}


class DocumentIdsPayload(BaseModel):
    docIds: List[str]
