"""Serving a file from disk by its path."""
import pytest
from fastapi import HTTPException
from fastapi.responses import FileResponse

from routers import files


def test_serves_a_file_with_its_type(tmp_path):
    path = tmp_path / "a.json"
    path.write_text("{}")

    answer = files.serve_file(str(path).lstrip("/"))

    assert isinstance(answer, FileResponse)
    assert answer.path == str(path)
    assert answer.media_type == "application/json"


def test_missing_file_is_404(tmp_path):
    with pytest.raises(HTTPException) as raised:
        files.serve_file(str(tmp_path / "nope.json").lstrip("/"))

    assert raised.value.status_code == 404


def test_a_folder_is_404(tmp_path):
    with pytest.raises(HTTPException) as raised:
        files.serve_file(str(tmp_path).lstrip("/"))

    assert raised.value.status_code == 404
