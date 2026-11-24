"""Finding notebooks in a project's folder."""
import types

import pytest
from fastapi import HTTPException

from load_helpers import _notebook
from project_api.project_api import ProjectApi
from routers import projects as routes


class TestDetectNotebooks:
    def test_registers_untracked_notebooks(self, fake_hera, tmp_path):
        db = fake_hera(["Alpha"], {"Alpha": []})
        (tmp_path / "a.ipynb").write_text("{}")
        (tmp_path / "sub").mkdir()
        (tmp_path / "sub" / "b.ipynb").write_text("{}")

        added = ProjectApi("Alpha").detect_notebooks(str(tmp_path))

        assert added == 2
        names = sorted(d["desc"]["datasourceName"] for d in db.documents["Alpha"])
        assert names == ["a", "b"]

    def test_skips_notebooks_already_tracked(self, fake_hera, tmp_path):
        notebook = tmp_path / "a.ipynb"
        notebook.write_text("{}")
        fake_hera(["Alpha"], {"Alpha": [_notebook(str(notebook))]})

        assert ProjectApi("Alpha").detect_notebooks(str(tmp_path)) == 0

    def test_skips_checkpoint_copies(self, fake_hera, tmp_path):
        checkpoints = tmp_path / ".ipynb_checkpoints"
        checkpoints.mkdir()
        (checkpoints / "a-checkpoint.ipynb").write_text("{}")
        fake_hera(["Alpha"], {"Alpha": []})

        assert ProjectApi("Alpha").detect_notebooks(str(tmp_path)) == 0

    def test_route_answers_the_project_and_the_count(self, fake_hera, warmed, tmp_path):
        fake_hera(["Alpha"], {"Alpha": []})
        (tmp_path / "a.ipynb").write_text("{}")

        answer = routes.detect(types.SimpleNamespace(
            projectName="Alpha", filesDirectory=str(tmp_path)))

        assert answer["addedNotebooks"] == 1
        assert answer["project"]["name"] == "Alpha"
        assert len(answer["project"]["documents"]) == 1

    def test_route_rejects_a_project_with_no_files_directory(self, fake_hera, warmed):
        fake_hera(["Alpha"], {"Alpha": []})

        with pytest.raises(HTTPException) as raised:
            routes.detect(types.SimpleNamespace(projectName="Alpha", filesDirectory=""))

        assert raised.value.status_code == 400
