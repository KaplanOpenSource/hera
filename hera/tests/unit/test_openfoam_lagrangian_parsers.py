"""openFoam/lagrangian/abstractLagrangianSolver.py: the module-level
pure-ish file parsers. readLagrangianRecord needs several real OpenFOAM
case files at once (globalPositions, origId, procId, ...) and is left
uncovered."""
import os

import pandas
import pytest

from hera.simulations.openFoam.lagrangian.abstractLagrangianSolver import (
    readEulerianConcentration,
    robustOpenFOAMFileValuesParser,
)


@pytest.mark.unit
class TestRobustOpenFOAMFileValuesParser:
    def test_a_missing_file_returns_an_empty_float_dataframe(self, tmp_path):
        result = robustOpenFOAMFileValuesParser(str(tmp_path / "nofile"), ["x", "y", "z"])
        assert len(result) == 0
        assert list(result.columns) == ["x", "y", "z"]
        assert all(dt == float for dt in result.dtypes)

    def test_the_uniform_shorthand_is_tiled_to_n_rows(self, tmp_path):
        path = tmp_path / "uniform.txt"
        path.write_text("header\n" * 17 + "3{1 2 3}\n")
        result = robustOpenFOAMFileValuesParser(str(path), ["x", "y", "z"])
        assert len(result) == 3
        assert list(result.iloc[0]) == [1.0, 2.0, 3.0]
        assert list(result.iloc[2]) == [1.0, 2.0, 3.0]

    def test_a_normal_openfoam_vector_list_is_parsed_per_line(self, tmp_path):
        path = tmp_path / "normal.txt"
        lines = ["h"] * 18 + ["(1 2 3)", "(4 5 6)", "(7 8 9)"]
        path.write_text("\n".join(lines) + "\n")
        result = robustOpenFOAMFileValuesParser(str(path), ["x", "y", "z"])
        assert list(result.iloc[1]) == [4.0, 5.0, 6.0]


@pytest.mark.unit
class TestReadEulerianConcentration:
    def test_a_missing_time_directory_returns_an_empty_frame(self, tmp_path):
        result = readEulerianConcentration("0.5", str(tmp_path))
        assert len(result) == 0
        assert list(result.columns) == ["x", "y", "z", "C", "datetime"]

    def test_a_present_file_is_read_and_stamped_with_the_time(self, tmp_path):
        time_dir = tmp_path / "0.5"
        time_dir.mkdir()
        (time_dir / "kinematicCloudEulerConcentrations").write_text("x,y,z,C\n1,2,3,0.5\n")
        result = readEulerianConcentration("0.5", str(tmp_path))
        assert result["datetime"].iloc[0] == pytest.approx(0.5)
        assert result["C"].iloc[0] == pytest.approx(0.5)

    def test_an_empty_file_returns_an_empty_frame_not_an_error(self, tmp_path):
        time_dir = tmp_path / "0.5"
        time_dir.mkdir()
        (time_dir / "kinematicCloudEulerConcentrations").write_text("")
        result = readEulerianConcentration("0.5", str(tmp_path))
        assert len(result) == 0

    def test_a_custom_cloud_name_is_used_in_the_filename(self, tmp_path):
        time_dir = tmp_path / "0.5"
        time_dir.mkdir()
        (time_dir / "myCloudEulerConcentrations").write_text("x,y,z,C\n1,2,3,0.5\n")
        result = readEulerianConcentration("0.5", str(tmp_path), cloudName="myCloud")
        assert len(result) == 1
