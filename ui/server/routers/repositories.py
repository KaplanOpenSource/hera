"""The repository entrypoints: the registered list, one file, a folder, loading.

Named routes replacing the last repository `/exec` calls the UI used to assemble
in the browser. See ui/client/ENTRYPOINTS_REPOSITORIES_STEPS.md.

A repository name is only a nickname for a file, so no name crosses the wire:
every route that takes a repository takes its path.
"""
from typing import Any, Dict, Optional

from fastapi import APIRouter, HTTPException
from fastapi.encoders import jsonable_encoder

from api_models import LoadRepositoryPayload
from dependencies import warmup
from project_api.project_api import ProjectApi
from repository_api.errors import NoSuchRepository, Refused
from repository_api.load_log import LoadLog
from repository_api.repo_folder import RepoFolder
from repository_api.repo_json import RepositoryJson
from repository_api.repo_json_file import RepoJsonFile
from repository_api.repository_table import RepositoryTable
from repository_api.toolkit_check import ToolkitCheck

router = APIRouter(prefix="/repository")


def _require_ready() -> None:
    """Every route here reaches hera, so none of them can answer while warming."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")


def _load(payload: LoadRepositoryPayload, projectName: str) -> Optional[dict]:
    """Load what the payload points at, and answer the tree it loaded.

    See ``LoadRepositoryPayload``. The answer is None for the registered
    repositories, which are many trees and not one.
    """
    if payload.path and payload.repositoryJSON is not None:
        raise Refused("give either path or repositoryJSON, not both")
    if payload.path:
        repository = RepoJsonFile(payload.path).repository()
        repository.load_into(projectName, payload.overwrite)
        return repository.tree
    if payload.repositoryJSON is not None:
        repository = RepositoryJson(payload.repositoryJSON, payload.baseDir)
        repository.load_into(projectName, payload.overwrite)
        return repository.tree
    RepositoryTable().load_all_into(projectName, payload.overwrite)
    return None


@router.get("/registered")
def registered() -> dict:
    """The registered repositories. See ``RepositoryTable.rows``."""
    _require_ready()
    return jsonable_encoder({"repositories": RepositoryTable().rows()})


@router.get("/get")
def get(path: str) -> dict:
    """One repository file, as it is on disk and with its paths resolved.

    A path that cannot be read answers 404 with the reason the file gave.
    """
    _require_ready()
    try:
        repository = RepoJsonFile(path).repository()
    except NoSuchRepository as missing:
        raise HTTPException(status_code=404, detail=str(missing)) from missing
    return jsonable_encoder({"json": repository.tree, "resolved": repository.resolved()})


@router.get("/jsonsInFolder")
def jsonsInFolder(folder: str, withContents: bool = False) -> dict:
    """The repository files under the folder, with their contents when asked."""
    _require_ready()
    repositories = RepoFolder(folder).repositories()
    answer: Dict[str, Any] = {"files": list(repositories.keys())}
    if withContents:
        answer["repoJsons"] = repositories
    return jsonable_encoder(answer)


@router.post("/loadIntoProject")
def loadIntoProject(payload: LoadRepositoryPayload) -> dict:
    """Load datasources into the project.

    Answers the project they joined, how many documents were added, and the
    problems on the way: a load that added nothing must be able to say why.
    """
    _require_ready()
    projectName = payload.projectName.strip()
    if not projectName:
        raise HTTPException(status_code=400, detail="no project name given")
    before = len(ProjectApi(projectName).documents()["documents"])
    # The tree that was loaded, which the registered repositories do not have.
    tree: Optional[dict] = None
    log = LoadLog()
    try:
        # hera logs its failures instead of raising, so the log is part of the answer.
        with log:
            tree = _load(payload, projectName)
    except NoSuchRepository as missing:
        raise HTTPException(status_code=404, detail=str(missing)) from missing
    except Refused as refused:
        raise HTTPException(status_code=400, detail=str(refused)) from refused
    project = ProjectApi(projectName).documents()
    added = len(project["documents"]) - before
    problems = log.problems
    # Nothing added is the case that needs a reason, and a toolkit that cannot
    # be built is the usual one: hera skips it without saying so.
    if added == 0 and tree is not None:
        problems = problems + ToolkitCheck(tree, projectName).problems()
    return jsonable_encoder({"project": project, "added": added, "problems": problems})
