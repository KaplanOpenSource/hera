"""The project list and whether a project exists."""
from project_api.project_api import DEFAULT_PROJECT, ProjectApi
from routers import projects as routes


class TestProjectList:
    def test_default_project_comes_first(self, fake_hera):
        fake_hera(["zeta", DEFAULT_PROJECT, "Alpha"])

        assert [p["name"] for p in ProjectApi.names()] == [DEFAULT_PROJECT, "Alpha", "zeta"]

    def test_no_projects_is_an_empty_list(self, fake_hera):
        fake_hera([])

        assert ProjectApi.names() == []

    def test_route_answers_the_list(self, fake_hera, warmed):
        fake_hera(["Alpha", DEFAULT_PROJECT])

        assert routes.names() == {"projects": [{"name": DEFAULT_PROJECT}, {"name": "Alpha"}]}

    def test_exists_only_for_a_project_in_the_list(self, fake_hera):
        fake_hera(["Alpha"])

        assert ProjectApi("Alpha").exists()
        assert not ProjectApi("Beta").exists()
