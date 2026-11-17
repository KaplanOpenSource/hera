"""OpenFOAM field I/O after the move from PyFoam to foamlib (#1169).

Each scenario in ``_openfoam_golden.py`` runs on the current code and is
compared with what the PyFoam implementation produced for it
(``hera/tests/expected/openfoam_fields/pyfoam``):

* written files: both are parsed with foamlib and compared value by value,
  and the ``uniform`` / ``nonuniform List<T>`` keywords are counted in the
  raw text, because a 3-long scalar list and a uniform vector parse to the
  same array;
* returned DataFrames: compared with ``compare_outputs``;
* scenarios that raised with PyFoam (``<name>.error``) are bugs fixed on
  purpose: each is listed in ``FIXED`` and has its own test of the values.
"""
import collections
import os
import re

import numpy
import pytest
from foamlib import FoamFieldFile

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome
from hera.tests.conftest import compare_outputs
from hera.tests.unit._openfoam_golden import (
    FRAME_SCENARIOS,
    PYFOAM_DIR,
    WRITE_SCENARIOS,
    frameToCSV,
    readCSV,
)

pytestmark = pytest.mark.unit

# Where the PyFoam output was wrong, so the foamlib output differs on purpose.
# Each one is checked by its own test below.
BAD_INTERNAL_FIELD = {
    "documented_api": "PyFoam wrote 'internalField (5 0 0);' without the 'uniform' keyword",
    "boundaries_from_case_single_internal_value": "PyFoam wrote 'internalField uniform 3(5 0 0);'",
}
DROPPED_CATCH_ALL = {
    "scalar_from_dataframe": "PyFoam's DictProxy leaves regex keys out of items(), so merging the "
                             "existing boundaryField dropped the '\".*\"' zeroGradient entry",
}

# Scenarios that raised with PyFoam and are fixed now: name -> the bug.
_SCALAR_ROWS = "a scalar field could not become a DataFrame: atleast_2d turned a list into one row"
_SCALAR_LEN = "a scalar field could not become a DataFrame: a uniform patch value has no len()"
FIXED = {
    "extract_p": _SCALAR_ROWS,
    "getDataFrame_scalar_from_dataframe": _SCALAR_ROWS,
    "readFieldAsDataFrame_parallel_p": _SCALAR_ROWS,
    "extract_T": _SCALAR_LEN,
    "getDataFrame_single_T": _SCALAR_LEN,
    "getDataFrame_parallel_U": "B290: readFromCase read logger.DEBUG and raised on every decomposed case",
    "roundtrip_parallel_U": "B290",
    "roundtrip_parallel_p": "B290",
}

_KEYWORD = re.compile(r"nonuniform\s+List<\w+>|(?<!non)uniform")


def _golden(kind, name):
    return os.path.join(PYFOAM_DIR, kind, name)


def _goldenError(kind, name):
    path = _golden(kind, name) + ".error"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return f.read().strip()


def _filesUnder(root):
    return sorted(os.path.relpath(os.path.join(d, f), root) for d, _, files in os.walk(root) for f in files)


def _read(path):
    with open(path) as f:
        return f.read()


def _assertSameValue(new, old, where):
    if isinstance(old, dict):
        assert isinstance(new, dict), where
        assert list(new) == list(old), f"{where}: keys {list(new)} != {list(old)}"
        for key in old:
            _assertSameValue(new[key], old[key], f"{where}/{key}")
    elif isinstance(old, str):
        assert new == old, f"{where}: {new!r} != {old!r}"
    else:
        newArray, oldArray = numpy.asarray(new, dtype=float), numpy.asarray(old, dtype=float)
        assert newArray.shape == oldArray.shape, f"{where}: shape {newArray.shape} != {oldArray.shape}"
        assert numpy.allclose(newArray, oldArray), f"{where}: {new!r} != {old!r}"


def _parse(path):
    content = FoamFieldFile(path).as_dict(include_header=True)
    header = content.pop("FoamFile")
    return header, content


def _cases(scenarios, kind):
    return [name for name in sorted(scenarios) if _goldenError(kind, name) is None]


@pytest.fixture()
def home():
    return OFObjectHome()


# --------------------------------------------------------------------------- #
#  Written files                                                               #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", _cases(WRITE_SCENARIOS, "write"))
def test_written_files_hold_the_values_pyfoam_wrote(name, home, tmp_path):
    WRITE_SCENARIOS[name](home, str(tmp_path))
    golden = _golden("write", name)
    assert _filesUnder(str(tmp_path)) == _filesUnder(golden)

    for rel in _filesUnder(golden):
        newHeader, newContent = _parse(os.path.join(tmp_path, rel))
        oldHeader, oldContent = _parse(os.path.join(golden, rel))
        assert (newHeader["class"], newHeader["object"]) == (oldHeader["class"], oldHeader["object"])
        if name in BAD_INTERNAL_FIELD:
            oldContent.pop("internalField")
            newContent.pop("internalField")
        if name in DROPPED_CATCH_ALL:
            newContent["boundaryField"].pop('".*"')
        _assertSameValue(newContent, oldContent, f"{name}/{rel}")

        if name not in BAD_INTERNAL_FIELD:
            newKeywords = collections.Counter(_KEYWORD.findall(_read(os.path.join(tmp_path, rel))))
            oldKeywords = collections.Counter(_KEYWORD.findall(_read(os.path.join(golden, rel))))
            assert newKeywords == oldKeywords, f"{name}/{rel}"


def test_the_documented_uniform_internal_value_is_written_with_the_uniform_keyword(home, tmp_path):
    WRITE_SCENARIOS["documented_api"](home, str(tmp_path))
    header, content = _parse(os.path.join(tmp_path, "0", "U"))
    assert re.search(r"internalField\s+uniform\s+\(5(\.0)? 0(\.0)? 0(\.0)?\)", _read(os.path.join(tmp_path, "0", "U")))
    assert numpy.allclose(content["internalField"], [5, 0, 0])


def test_an_internal_value_read_from_a_case_is_written_as_a_uniform_vector(home, tmp_path):
    WRITE_SCENARIOS["boundaries_from_case_single_internal_value"](home, str(tmp_path))
    _, content = _parse(os.path.join(tmp_path, "0", "U"))
    assert isinstance(content["internalField"], numpy.ndarray)
    assert numpy.allclose(content["internalField"], [5, 0, 0])


def test_setting_values_from_a_dataframe_keeps_the_catch_all_boundary(home, tmp_path):
    WRITE_SCENARIOS["scalar_from_dataframe"](home, str(tmp_path))
    _, content = _parse(os.path.join(tmp_path, "0", "p"))
    assert content["boundaryField"]['".*"'] == {"type": "zeroGradient"}


def test_regex_patch_names_are_quoted(home, tmp_path):
    """OpenFOAM matches a boundaryField key as a regex only when it is quoted."""
    WRITE_SCENARIOS["empty_vector_parallel"](home, str(tmp_path))
    text = _read(os.path.join(tmp_path, "processor0", "0", "U"))
    assert '"proc.*"' in text and '".*"' in text
    assert '""' not in text


def test_a_three_long_scalar_list_is_not_written_as_a_vector(home, tmp_path):
    WRITE_SCENARIOS["scalar_from_dataframe"](home, str(tmp_path))
    text = _read(os.path.join(tmp_path, "0", "p"))
    assert re.search(r"nonuniform\s+List<scalar>\s+3\s*\(\s*7", text)


# --------------------------------------------------------------------------- #
#  Returned DataFrames                                                         #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", _cases(FRAME_SCENARIOS, "frames"))
def test_dataframes_match_pyfoam(name, home):
    result = readCSV(frameToCSV(FRAME_SCENARIOS[name](home)))
    expected = readCSV(_read(_golden("frames", name) + ".csv"))
    assert compare_outputs(result, expected, "dataframe"), f"{name}:\n{result}\n!=\n{expected}"


def test_every_scenario_that_raised_with_pyfoam_is_fixed_and_listed():
    """A scenario with a <name>.error golden has no values to compare with, so it
    must be in FIXED, which points at the test that checks its values."""
    raised = ({name for name in WRITE_SCENARIOS if _goldenError("write", name)}
              | {name for name in FRAME_SCENARIOS if _goldenError("frames", name)})
    assert raised == set(FIXED)


# --------------------------------------------------------------------------- #
#  Scalar fields as DataFrames (raised with PyFoam)                            #
# --------------------------------------------------------------------------- #
def _values(df, region, column, boundary=None):
    rows = df[df["region"] == region]
    if boundary is not None:
        rows = rows[rows["boundary"] == boundary]
    return rows.sort_values("processorIndex")[column].tolist()


def test_a_nonuniform_scalar_field_gives_one_row_per_cell_and_face(home):
    df = FRAME_SCENARIOS["extract_p"](home)
    assert _values(df, "internalField", "p") == [0.5, 0.25, -0.125, 1e-05]
    assert _values(df, "boundaryField", "p", "outlet") == [0.0, 0.5, 1.0]
    for patch in ("inlet", "walls", "sym", "frontAndBack"):
        assert numpy.isnan(_values(df, "boundaryField", "p", patch)).all()


def test_a_uniform_scalar_field_gives_one_row_per_region(home):
    df = FRAME_SCENARIOS["extract_T"](home)
    assert _values(df, "internalField", "T") == [300.0]
    assert _values(df, "boundaryField", "T", "inlet") == [310.0]
    assert _values(df, "boundaryField", "T", "outlet") == [300.0]


def test_a_scalar_field_read_from_a_case_gives_the_same_frame_as_extractFieldFile(home):
    viaField = readCSV(frameToCSV(FRAME_SCENARIOS["getDataFrame_single_T"](home)))
    viaFile = readCSV(frameToCSV(FRAME_SCENARIOS["extract_T"](home)))
    assert compare_outputs(viaField, viaFile, "dataframe")


def test_a_scalar_field_set_from_a_dataframe_gives_the_values_back(home):
    df = FRAME_SCENARIOS["getDataFrame_scalar_from_dataframe"](home)
    assert _values(df, "internalField", "p") == [1.0, 2.0, 3.0, 4.0]
    assert _values(df, "boundaryField", "p", "inlet") == [5.0, 6.0]
    assert _values(df, "boundaryField", "p", "outlet") == [7.0, 8.0, 9.0]


def test_a_decomposed_scalar_field_is_read_per_processor(home):
    df = FRAME_SCENARIOS["readFieldAsDataFrame_parallel_p"](home)
    assert _values(df[df["processor"] == 0], "internalField", "p") == [0.5, 0.25]
    assert _values(df[df["processor"] == 1], "internalField", "p") == [-0.125, 1e-05]
    assert _values(df[df["processor"] == 1], "boundaryField", "p", "outlet") == [0.0, 0.5, 1.0]
    assert _values(df[df["processor"] == 0], "boundaryField", "p", "procBoundary0to1") == [-0.125, 1e-05]


# --------------------------------------------------------------------------- #
#  Decomposed cases (B290, raised with PyFoam)                                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name, fileName", [("roundtrip_parallel_U", "U"), ("roundtrip_parallel_p", "p")])
def test_a_decomposed_field_read_and_written_back_keeps_its_values(name, fileName, home, tmp_path):
    WRITE_SCENARIOS[name](home, str(tmp_path))
    from hera.tests.unit._openfoam_golden import INPUT_CASE
    for processor in ("processor0", "processor1"):
        _, written = _parse(os.path.join(tmp_path, processor, "0", fileName))
        _, original = _parse(os.path.join(INPUT_CASE, processor, "0", fileName))
        _assertSameValue(written, original, f"{name}/{processor}")


def test_a_decomposed_field_gives_the_same_frame_as_reading_each_processor(home):
    viaField = FRAME_SCENARIOS["getDataFrame_parallel_U"](home)
    viaFiles = FRAME_SCENARIOS["readFieldAsDataFrame_parallel_U"](home)
    viaFiles = viaFiles[~viaFiles["boundary"].fillna("").str.startswith("proc")]
    assert compare_outputs(readCSV(frameToCSV(viaField)), readCSV(frameToCSV(viaFiles)).drop(columns="time"),
                           "dataframe")
