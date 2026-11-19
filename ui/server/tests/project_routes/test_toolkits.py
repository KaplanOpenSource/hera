"""The toolkit documents."""
from project_api.toolkit_documents import toolkit_documents
from routers import projects as routes


class TestToolkits:
    def test_answers_the_toolkit_documents(self, fake_hera, warmed):
        fake_hera(["Alpha"])

        assert toolkit_documents()[0]["toolkit"] == "LSM"
        assert routes.toolkits()["toolkitDocs"][0]["toolkit"] == "LSM"
