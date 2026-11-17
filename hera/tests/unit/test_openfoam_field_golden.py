"""OpenFOAM field I/O after the move from PyFoam to foamlib (#1169).

Each scenario in ``_openfoam_golden.py`` runs on the current code and is
compared with what the PyFoam implementation produced for it
(``hera/tests/expected/openfoam_fields/pyfoam``):

* written files: both are parsed with foamlib and compared value by value,
  and the ``uniform`` / ``nonuniform List<T>`` keywords are counted in the
  raw text, because a 3-long scalar list and a uniform vector parse to the
  same array;
* returned DataFrames: compared with ``compare_outputs``;
* scenarios that raised with PyFoam (``<name>.error``) must still raise the
  same exception, unless the bug was fixed on purpose (``FIXED`` below,
  each with its own test of the corrected values).
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
FIXED = {}

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


def _failedCases(scenarios, kind):
    return [name for name in sorted(scenarios) if _goldenError(kind, name) is not None and name not in FIXED]


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


@pytest.mark.parametrize("name", _failedCases(WRITE_SCENARIOS, "write"))
def test_write_scenarios_that_failed_with_pyfoam_still_fail(name, home, tmp_path):
    with pytest.raises(Exception) as caught:
        WRITE_SCENARIOS[name](home, str(tmp_path))
    assert type(caught.value).__name__ == _goldenError("write", name)


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


@pytest.mark.parametrize("name", _failedCases(FRAME_SCENARIOS, "frames"))
def test_frame_scenarios_that_failed_with_pyfoam_still_fail(name, home):
    with pytest.raises(Exception) as caught:
        FRAME_SCENARIOS[name](home)
    assert type(caught.value).__name__ == _goldenError("frames", name)
