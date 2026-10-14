"""The project entrypoints: the list, one project, toolkits, create, delete, notebooks.

Named routes replacing the project `/exec` calls the UI used to assemble in the
browser. See ui/client/ENTRYPOINTS_PROJECTS_STEPS.md.
"""
from fastapi import APIRouter, HTTPException
from fastapi.encoders import jsonable_encoder

from api_models import CreateProjectPayload, DeleteProjectPayload, DetectNotebooksPayload
from dependencies import warmup
from project_api.project_api import ProjectApi
from project_api.toolkit_documents import toolkit_documents

router = APIRouter(prefix="/project")


def _require_ready() -> None:
    """Every route here reaches hera, so none of them can answer while warming."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")


def _open_project(projectName: str) -> ProjectApi:
    """The project the request names, refusing one the database does not have."""
    _require_ready()
    if not projectName:
        raise HTTPException(status_code=400, detail="no project name given")
    project = ProjectApi(projectName)
    if not project.exists():
        raise HTTPException(status_code=404, detail=f"no such project: {projectName}")
    return project


@router.get("/names")
def names() -> dict:
    """The project list, the default project first. See ``ProjectApi.names``."""
    _require_ready()
    return {"projects": ProjectApi.names()}


@router.get("/documents")
def documents(projectName: str) -> dict:
    """One project and its documents. See ``ProjectApi.documents``."""
    project = _open_project(projectName)
    return jsonable_encoder({"project": project.documents()})


@router.get("/toolkits")
def toolkits() -> dict:
    """The toolkit documents. See ``project_api.toolkit_documents``."""
    _require_ready()
    return jsonable_encoder({"toolkitDocs": toolkit_documents()})


@router.post("/create")
def create(payload: CreateProjectPayload) -> dict:
    """Create the project, register the notebooks already in its files directory."""
    _require_ready()
    projectName = payload.projectName.strip()
    if not projectName:
        raise HTTPException(status_code=400, detail="no project name given")
    project = ProjectApi(projectName)
    filesDirectory = project.create(payload.filesDirectory, payload.loadRepositories)
    addedNotebooks = project.detect_notebooks(filesDirectory)
    # The list and the project itself, so creating costs one call: the client
    # shows the new project without loading it again.
    return jsonable_encoder({
        "projects": ProjectApi.names(),
        "project": project.documents(),
        "addedNotebooks": addedNotebooks,
    })


@router.post("/delete")
def delete(payload: DeleteProjectPayload) -> dict:
    """Delete the project, and answer with the project list that is left."""
    project = _open_project(payload.projectName)
    project.delete(payload.deleteFiles)
    remaining = ProjectApi.names()
    nextProject = None
    if remaining:
        nextProject = ProjectApi(remaining[0]["name"]).documents()
    return jsonable_encoder({"projects": remaining, "project": nextProject})


@router.post("/detectNotebooks")
def detect(payload: DetectNotebooksPayload) -> dict:
    """Register the notebooks under the project's files directory that are untracked."""
    project = _open_project(payload.projectName)
    if not payload.filesDirectory:
        raise HTTPException(status_code=400, detail="the project has no files directory")
    addedNotebooks = project.detect_notebooks(payload.filesDirectory)
    return jsonable_encoder({
        "project": project.documents(),
        "addedNotebooks": addedNotebooks,
    })
