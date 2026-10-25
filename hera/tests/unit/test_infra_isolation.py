"""Each test must start from an empty database."""
import pytest


@pytest.mark.unit
def test_writes_a_document(unit_project):
    unit_project.addMeasurementsDocument(
        resource="/leak-check", dataFormat="string", type="LeakCheck", desc={}
    )
    assert len(unit_project.getMeasurementsDocuments(type="LeakCheck")) == 1


@pytest.mark.unit
def test_previous_documents_are_gone(unit_project):
    """Must ALWAYS see zero, whichever order the two tests run in.

    Length, not equality against []: the datalayer returns a mongoengine
    QuerySet, which is never == to a list even when both are empty.
    """
    assert len(unit_project.getMeasurementsDocuments(type="LeakCheck")) == 0


@pytest.mark.unit
def test_timeseries_factory_is_deterministic():
    from hera.tests.unit import _factories

    first = _factories.timeseries(n=5, seed=7)
    second = _factories.timeseries(n=5, seed=7)
    assert first.equals(second)
    assert list(first.columns) == ["u", "v", "w"]
    assert len(first) == 5


@pytest.mark.unit
def test_elevation_grid_shape():
    from hera.tests.unit import _factories

    grid = _factories.elevation_grid(nx=3, ny=2)
    assert len(grid) == 6
    assert set(grid.columns) == {"X", "Y", "Elevation"}
