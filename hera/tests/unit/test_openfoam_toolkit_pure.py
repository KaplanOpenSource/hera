"""OFToolkit: the filesystem/parsing helpers that don't need a real solver run."""
import numpy
import pytest

from hera.simulations.openFoam.toolkit import OFToolkit


def _write_points_file(path, points):
    lines = ["FoamFile", "{", "}", str(len(points)), "("]
    lines += [f"({x} {y} {z})" for x, y, z in points]
    lines.append(")")
    path.write_text("\n".join(lines))


@pytest.mark.unit
class TestProcessorList:
    def test_no_processor_directories_gives_an_empty_list(self, tmp_path):
        assert OFToolkit.processorList(None, str(tmp_path)) == []

    def test_processor_directories_are_listed_by_basename(self, tmp_path):
        (tmp_path / "processor0").mkdir()
        (tmp_path / "processor1").mkdir()
        (tmp_path / "notAProcessor").mkdir()
        result = OFToolkit.processorList(None, str(tmp_path))
        assert sorted(result) == ["processor0", "processor1"]


@pytest.mark.unit
class TestReadPointsFile:
    def test_it_parses_coordinates_after_the_count_line(self, tmp_path):
        path = tmp_path / "points"
        _write_points_file(path, [(0, 0, 0), (1, 0, 0), (1, 1, 0)])
        result = OFToolkit.read_points_file(None, str(path))
        assert numpy.array_equal(result, numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))

    def test_it_stops_at_the_closing_paren(self, tmp_path):
        path = tmp_path / "points"
        _write_points_file(path, [(0, 0, 0)])
        result = OFToolkit.read_points_file(None, str(path))
        assert len(result) == 1


@pytest.mark.unit
class TestGetMeshExtent:
    def test_it_raises_when_the_points_file_is_missing(self, tmp_path):
        toolkit = OFToolkit.__new__(OFToolkit)
        with pytest.raises(FileNotFoundError):
            toolkit.getMeshExtent(str(tmp_path))

    def test_it_returns_the_bounding_box_of_the_points(self, tmp_path):
        case = tmp_path
        mesh_dir = case / "constant" / "polyMesh"
        mesh_dir.mkdir(parents=True)
        _write_points_file(mesh_dir / "points", [(0, 0, 0), (2, 3, 1), (1, 1, 1)])
        toolkit = OFToolkit.__new__(OFToolkit)
        bounds = toolkit.getMeshExtent(str(case))
        assert bounds == {"x": (0.0, 2.0), "y": (0.0, 3.0), "z": (0.0, 1.0)}


@pytest.mark.unit
class TestGetTimeListSingleProcessorIsBroken:
    @pytest.fixture()
    def case(self, tmp_path):
        (tmp_path / "0").mkdir()
        (tmp_path / "10").mkdir()
        (tmp_path / "constant").mkdir()
        (tmp_path / "system").mkdir()
        return tmp_path

    @staticmethod
    def _toolkit_pointing_at_directory(directory):
        toolkit = OFToolkit.__new__(OFToolkit)
        toolkit.getWorkflowDocumentFromDB = lambda **kw: []
        return toolkit

    def test_it_lists_the_time_step_directories(self, case):
        toolkit = self._toolkit_pointing_at_directory(case)
        result = toolkit.getTimeList(str(case), singleProcessor=True)
        assert result == [0.0, 10.0]

    def test_a_missing_case_raises(self, tmp_path):
        toolkit = self._toolkit_pointing_at_directory(tmp_path)
        with pytest.raises(ValueError, match="not found"):
            toolkit.getTimeList(str(tmp_path / "does_not_exist"))
