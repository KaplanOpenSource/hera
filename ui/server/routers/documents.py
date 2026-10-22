"""The document entrypoints: create, update, delete.

Named routes replacing the document `/exec` calls the UI used to assemble in the
browser. See ui/client/ENTRYPOINTS_DOCUMENTS_STEPS.md.
"""
from fastapi import APIRouter, HTTPException
from fastapi.encoders import jsonable_encoder

from api_models import CreateDocumentPayload, DocumentIdsPayload, UpdateDocumentPayload
from dependencies import warmup
from document_api.document_api import DocumentApi, NoSuchDocument, Refused
from document_api.notebook_file import NotebookFile
from project_api.project_api import ProjectApi

router = APIRouter(prefix="/document")


def _require_ready() -> None:
    """Every route here reaches hera, so none of them can answer while warming."""
    if not warmup.ready:
        raise HTTPException(status_code=503, detail="WARMING_UP")


def _project(projectName: str) -> dict:
    """The whole project, the same shape /project/documents answers."""
    return ProjectApi(projectName).documents()


@router.post("/create")
def create(payload: CreateDocumentPayload) -> dict:
    """Add one document, and answer it together with the project it joined."""
    _require_ready()
    projectName = payload.projectName.strip()
    if not projectName:
        raise HTTPException(status_code=400, detail="no project name given")
    # A notebook document owns a file; an existing one is always left alone.
    if payload.notebook is not None:
        if not isinstance(payload.resource, str) or not payload.resource:
            raise HTTPException(status_code=400, detail="a notebook needs a resource path")
        NotebookFile(payload.resource).write_if_missing(
            projectName, payload.notebook.toolkitNames)
    try:
        document = DocumentApi.create(
            projectName=projectName,
            collection=payload.collection,
            type=payload.type,
            resource=payload.resource,
            dataFormat=payload.dataFormat,
            desc=payload.desc)
    except Refused as refused:
        raise HTTPException(status_code=400, detail=str(refused))
    # The id at the top level, so the caller need not dig for it.
    return jsonable_encoder({
        "document": document,
        "project": _project(projectName),
        "addedDocId": document.get("_id", {}).get("$oid"),
    })


@router.post("/update")
def update(payload: UpdateDocumentPayload) -> dict:
    """Write the given fields, and answer the document and its project."""
    _require_ready()
    if not payload.docId:
        raise HTTPException(status_code=400, detail="no document id given")
    document = DocumentApi(payload.docId)
    try:
        projectName = document.project_name()
        saved = document.update(payload.fields)
    except NoSuchDocument as missing:
        raise HTTPException(status_code=404, detail=str(missing))
    except Refused as refused:
        raise HTTPException(status_code=400, detail=str(refused))
    return jsonable_encoder({"document": saved, "project": _project(projectName)})


@router.post("/delete")
def delete(payload: DocumentIdsPayload) -> dict:
    """Remove the given documents, and answer the project that is left."""
    _require_ready()
    if not payload.docIds:
        raise HTTPException(status_code=400, detail="no document ids given")
    try:
        # Read the project before deleting: afterwards the documents are gone.
        projectName = DocumentApi.one_project(payload.docIds)
        for docId in payload.docIds:
            DocumentApi(docId).delete()
    except NoSuchDocument as missing:
        raise HTTPException(status_code=404, detail=str(missing))
    except Refused as refused:
        raise HTTPException(status_code=400, detail=str(refused))
    return jsonable_encoder({"project": _project(projectName)})
