"""The node catalog route, read from the real Hermes submodule."""
import node_catalog
from routers import nodes


def test_route_answers_the_catalog(monkeypatch):
    catalog = [{"type": "general.Parameters", "parameters": []}]
    monkeypatch.setattr(nodes, "get_node_catalog", lambda: catalog)

    assert nodes.node_catalog() == catalog


def test_catalog_lists_hermes_nodes_sorted_by_type():
    catalog = node_catalog.get_node_catalog()

    types = [entry["type"] for entry in catalog]
    assert types
    assert types == sorted(types)


def test_each_entry_has_its_parameters():
    catalog = node_catalog.get_node_catalog()

    for entry in catalog:
        assert isinstance(entry["parameters"], list)
        for parameter in entry["parameters"]:
            assert "name" in parameter
            assert "is_required" in parameter
