# Basic tests for the shapes toolkit

import pathlib

import pandas as pd
import pytest

from hera.measurements.GIS.shapes import ShapesToolKit
from hera import toolkit


@pytest.fixture
def shapes_toolkit() -> ShapesToolKit:
    return toolkit.ToolkitHome().getToolkit("GIS_Shapes", projectName="test_hera")


@pytest.fixture
def shapes_file():
    doc_folder = pathlib.Path(__file__).parent.parent.parent / "doc" / "source"
    abs_path = (doc_folder / "measurements" / "GIS" / "locations" / "examples" / "Katsrin.shp").absolute()
    return str(abs_path)


def test_just_load_file_does_not_save(shapes_toolkit, shapes_file):
    shapes_toolkit.loadData(shapes_file, regionName="Narnia")
    # Now try to get the region from the toolkit
    shape = shapes_toolkit.getShape("Narnia")
    assert shape is None


def test_just_load_file(data_regression, shapes_toolkit, shapes_file):
    shapes = shapes_toolkit.loadData(shapes_file)
    data = shapes.getData()
    data.sort_values(by=['SHAPE_Leng'], inplace=True)
    # Pick up the 10 shortest shapes, and the one longest
    d = pd.concat([data[:10], data[-1:]])
    data_regression.check(d.to_json())


def test_load_file_and_save(cached_data_regression, shapes_toolkit, shapes_file):
    shapes_toolkit.loadData(shapes_file, regionName="Narnia", saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB)
    # Now try to get the region from the toolkit
    shape = shapes_toolkit.getShape("Narnia")
    cached_data_regression.check(shape)
