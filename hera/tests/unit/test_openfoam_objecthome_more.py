"""OFObjectHome: the three case-reading methods batches 11-31 didn't reach.

``test_openfoam_objects.py`` covers the dimension vector, the field-definition
registry, ``getEmptyField`` and the ``pandasToFoamFormat`` defect.  What was
left uncovered is everything that takes a *case directory*:

* ``getEmptyFieldFromCase`` -- builds an uninitialized ``OFField`` and lets it
  read the patch list out of the case;
* ``readFieldFromCase`` -- the same, plus the field values;
* ``readFieldAsDataFrame`` -- the deprecated reader that walks
  ``processor*``/``<time>`` and concatenates one ``extractFieldFile`` per pair.

PyFoam is the conftest's bare ``MagicMock`` module (see ``_stubs.py``), so
``ParsedParameterFile``/``BoundaryDict``/``WriteParameterFile`` carry no real
storage and invent no values.  Nothing here asserts a value a stub produced.
The first two methods are therefore asserted on what *they* decide -- the
field's name, its on-disk file name, its type, the dimension dictionary
selected for the flow type, whether the field came out parallel, and which
processor keys it holds -- all of which OFObjectHome computes itself from its
catalog and from the directories on disk.  ``readFieldAsDataFrame`` is
asserted through a monkeypatched ``extractFieldFile`` seam: which paths it
asks for, in which order, with which column names and which per-call keyword
arguments.  Its return value is the concatenation of whatever that seam
returns, so only the shape of the concatenation is asserted.

Bugs pinned here (each with a strict xfail for the intended behaviour and a
passing characterisation of what happens today):

* B290 ``readFieldFromCase`` cannot read a decomposed case at all:
  ``OFField.readFromCase``'s parallel branch tests
  ``logger.getEffectiveLevel() == logger.DEBUG``, and ``logger`` is a
  ``logging.Logger`` instance, which has no ``DEBUG`` attribute (the constant
  lives on the ``logging`` module).  Every parallel read dies with
  AttributeError before a single processor file is opened.
* B291 ``readFieldAsDataFrame(readParallel=False, times=None)`` finds no
  time steps unless the process happens to be sitting inside the case: the
  serial branch calls ``os.path.isdir(x)`` on the bare directory *name*
  instead of on ``os.path.join(finalCasePath, x)`` -- the parallel branch two
  blocks above joins it correctly.  With nothing discovered, the trailing
  ``pandas.concat([])`` raises "No objects to concatenate".
"""
import os
import sys

import pandas
import pytest

from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome

_SCALAR_TIME_DIRS = ("0", "100", "constant")


@pytest.fixture()
def home():
    return OFObjectHome()


@pytest.fixture()
def home_module():
    """The OFObjectHome *module*.

    ``preprocessOFObjects/__init__.py`` does ``from .OFObjectHome import
    OFObjectHome``, which shadows the module attribute on the package with the
    class of the same name, so attribute access cannot reach the module.
    """
    import hera.simulations.openFoam.preprocessOFObjects.OFObjectHome  # noqa: F401

    return sys.modules["hera.simulations.openFoam.preprocessOFObjects.OFObjectHome"]


@pytest.fixture()
def serialCase(tmp_path):
    """A case with time directories but no processor* directories."""
    case = tmp_path / "serialCase"
    for name in _SCALAR_TIME_DIRS:
        (case / name).mkdir(parents=True)
    return case


@pytest.fixture()
def parallelCase(tmp_path):
    """A two-processor decomposed case with two time steps each."""
    case = tmp_path / "parallelCase"
    for processor in ("processor0", "processor1"):
        for name in _SCALAR_TIME_DIRS:
            (case / processor / name).mkdir(parents=True)
    return case


@pytest.fixture()
def extractRecorder(home_module, monkeypatch):
    """Replace extractFieldFile with a recorder; return the call list."""
    calls = []

    def recorder(path, columnNames=None, **kwargs):
        calls.append(dict(path=path, columnNames=columnNames, **kwargs))
        return pandas.DataFrame({(columnNames or ["value"])[0]: [1.0]})

    monkeypatch.setattr(home_module, "extractFieldFile", recorder)
    return calls


# ---------------------------------------------------------------------------
# getEmptyFieldFromCase
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetEmptyFieldFromCase:
    def test_the_field_carries_the_name_it_was_asked_for(self, home, serialCase):
        field = home.getEmptyFieldFromCase(fieldName="p", flowType="incompressible",
                                           caseDirectory=str(serialCase))
        assert field.name == "p"

    def test_a_field_whose_file_is_named_differently_keeps_both_names(self, home, serialCase):
        """"cellCenters" is written to a file called "C"."""
        field = home.getEmptyFieldFromCase(fieldName="cellCenters", flowType="incompressible",
                                           caseDirectory=str(serialCase))
        assert (field.name, field.fileName) == ("cellCenters", "C")
        assert field.componentNames == ["Cx", "Cy", "Cz"]

    def test_the_field_type_comes_from_the_catalog(self, home, serialCase):
        vector = home.getEmptyFieldFromCase(fieldName="U", flowType="incompressible",
                                            caseDirectory=str(serialCase))
        scalar = home.getEmptyFieldFromCase(fieldName="T", flowType="incompressible",
                                            caseDirectory=str(serialCase))
        assert (vector.fieldType, scalar.fieldType) == ("vector", "scalar")

    @pytest.mark.parametrize("flowType, expected", [
        ("incompressible", {"m": 2, "s": -2}),
        ("compressible", {"kg": 1, "m": -1, "s": -2}),
    ])
    def test_the_pressure_dimensions_follow_the_flow_type(self, home, serialCase,
                                                          flowType, expected):
        field = home.getEmptyFieldFromCase(fieldName="p", flowType=flowType,
                                           caseDirectory=str(serialCase))
        assert field.dimensions == expected

    def test_a_field_with_one_set_of_units_ignores_the_flow_type(self, home, serialCase):
        incompressible = home.getEmptyFieldFromCase(fieldName="U", flowType="incompressible",
                                                    caseDirectory=str(serialCase))
        compressible = home.getEmptyFieldFromCase(fieldName="U", flowType="compressible",
                                                  caseDirectory=str(serialCase))
        assert incompressible.dimensions == compressible.dimensions == {"m": 1, "s": -1}

    def test_a_case_without_processor_directories_comes_back_single_processor(
            self, home, serialCase):
        field = home.getEmptyFieldFromCase(fieldName="p", flowType="incompressible",
                                           caseDirectory=str(serialCase))
        assert field.parallel is False
        assert list(field.data) == ["singleProcessor"]

    def test_a_decomposed_case_gets_one_entry_per_processor(self, home, parallelCase):
        field = home.getEmptyFieldFromCase(fieldName="U", flowType="incompressible",
                                           caseDirectory=str(parallelCase))
        assert field.parallel is True
        assert sorted(field.data) == ["processor0", "processor1"]

    def test_read_parallel_false_reads_a_decomposed_case_as_a_single_one(
            self, home, parallelCase):
        field = home.getEmptyFieldFromCase(fieldName="U", flowType="incompressible",
                                           caseDirectory=str(parallelCase),
                                           readParallel=False)
        assert field.parallel is False
        assert list(field.data) == ["singleProcessor"]

    def test_an_unknown_field_lists_the_catalog(self, home, serialCase):
        with pytest.raises(ValueError, match="Field noSuchField not found"):
            home.getEmptyFieldFromCase(fieldName="noSuchField", flowType="incompressible",
                                       caseDirectory=str(serialCase))

    def test_the_catalog_is_named_in_the_error_so_the_caller_can_pick_one(self, home, serialCase):
        with pytest.raises(ValueError, match="U,p,p_rgh"):
            home.getEmptyFieldFromCase(fieldName="noSuchField", flowType="incompressible",
                                       caseDirectory=str(serialCase))

    def test_the_field_computation_is_taken_from_the_catalog(self, home, serialCase):
        field = home.getEmptyFieldFromCase(fieldName="p", flowType="incompressible",
                                           caseDirectory=str(serialCase))
        assert field.fieldComputation == "eulerian"


# ---------------------------------------------------------------------------
# readFieldFromCase
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestReadFieldFromCase:
    def test_a_serial_case_is_read_into_a_single_processor_entry(self, home, serialCase):
        field = home.readFieldFromCase(fieldName="p", flowType="incompressible",
                                       caseDirectory=str(serialCase))
        assert field.parallel is False
        assert list(field.data) == ["singleProcessor"]

    def test_the_metadata_is_resolved_exactly_as_for_an_empty_field(self, home, serialCase):
        field = home.readFieldFromCase(fieldName="cellCenters", flowType="incompressible",
                                       caseDirectory=str(serialCase))
        assert (field.name, field.fileName, field.fieldType) == ("cellCenters", "C", "vector")
        assert field.dimensions == {"m": 1}

    def test_an_unknown_field_is_refused_before_the_case_is_touched(self, home, tmp_path):
        with pytest.raises(ValueError, match="Field noSuchField not found"):
            home.readFieldFromCase(fieldName="noSuchField", flowType="incompressible",
                                   caseDirectory=str(tmp_path / "doesNotExist"))

    def test_read_parallel_false_reads_a_decomposed_case_as_a_single_one(
            self, home, parallelCase):
        field = home.readFieldFromCase(fieldName="p", flowType="incompressible",
                                       caseDirectory=str(parallelCase), readParallel=False)
        assert field.parallel is False
        assert list(field.data) == ["singleProcessor"]

    @pytest.mark.xfail(
        strict=True,
        reason="B290: OFField.readFromCase's parallel branch reads "
               "`logger.DEBUG`, but `logger` is a logging.Logger instance and "
               "the DEBUG constant lives on the logging module, so every "
               "attempt to read a decomposed case raises AttributeError "
               "before any processor file is opened.  (The `itrObj` the line "
               "computes is then never used either -- the loop below it "
               "iterates procPaths.)  See the consolidated findings issue.",
    )
    def test_a_decomposed_case_should_be_readable(self, home, parallelCase):
        field = home.readFieldFromCase(fieldName="p", flowType="incompressible",
                                       caseDirectory=str(parallelCase))
        assert sorted(field.data) == ["processor0", "processor1"]

    def test_a_decomposed_case_currently_raises_attributeerror(self, home, parallelCase):
        """Characterisation of B290."""
        with pytest.raises(AttributeError, match="'Logger' object has no attribute 'DEBUG'"):
            home.readFieldFromCase(fieldName="p", flowType="incompressible",
                                   caseDirectory=str(parallelCase))

    def test_the_same_case_read_as_a_single_processor_succeeds(self, home, parallelCase):
        """Characterisation of B290: only the parallel branch is broken."""
        assert home.readFieldFromCase(fieldName="p", flowType="incompressible",
                                      caseDirectory=str(parallelCase),
                                      readParallel=False).parallel is False


# ---------------------------------------------------------------------------
# readFieldAsDataFrame
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestReadFieldAsDataFrameParallel:
    def test_one_file_is_read_per_processor_and_time_pair(self, home, parallelCase,
                                                          extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=[0, 100])
        assert sorted(call["path"] for call in extractRecorder) == sorted(
            os.path.join(str(parallelCase), processor, time, "p")
            for processor in ("processor0", "processor1") for time in ("0", "100"))

    def test_the_processor_column_carries_the_index_as_an_integer(self, home, parallelCase,
                                                                  extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=0)
        assert sorted(call["processor"] for call in extractRecorder) == [0, 1]
        assert all(isinstance(call["processor"], int) for call in extractRecorder)

    def test_the_time_is_forwarded_as_it_was_given(self, home, parallelCase, extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=100)
        assert {call["time"] for call in extractRecorder} == {100}

    def test_the_column_names_are_the_field_components(self, home, parallelCase,
                                                       extractRecorder):
        home.readFieldAsDataFrame(fieldName="cellCenters", caseDirectory=str(parallelCase),
                                  times=0)
        assert {tuple(call["columnNames"]) for call in extractRecorder} == {("Cx", "Cy", "Cz")}

    def test_the_file_name_of_the_field_is_used_not_its_logical_name(self, home, parallelCase,
                                                                     extractRecorder):
        home.readFieldAsDataFrame(fieldName="cellCenters", caseDirectory=str(parallelCase),
                                  times=0)
        assert all(os.path.basename(call["path"]) == "C" for call in extractRecorder)

    def test_times_none_discovers_the_numeric_time_directories_only(self, home, parallelCase,
                                                                    extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=None)
        discovered = {os.path.basename(os.path.dirname(call["path"])) for call in extractRecorder}
        assert discovered == {"0", "100"}
        assert "constant" not in discovered

    def test_the_result_is_the_concatenation_of_every_pair(self, home, parallelCase,
                                                           extractRecorder):
        frame = home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase),
                                          times=[0, 100])
        # 2 processors x 2 times x 1 row from the seam
        assert len(frame) == 4

    def test_filter_internal_patches_is_forwarded(self, home, parallelCase, extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=0,
                                  filterInternalPatches=True)
        assert all(call["filterInternalPatches"] is True for call in extractRecorder)

    def test_the_default_keeps_the_internal_processor_patches(self, home, parallelCase,
                                                             extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(parallelCase), times=0)
        assert all(call["filterInternalPatches"] is False for call in extractRecorder)

    def test_a_case_with_no_processor_directories_says_so(self, home, serialCase,
                                                          extractRecorder):
        with pytest.raises(ValueError, match="no processor\\* directories"):
            home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase))
        assert extractRecorder == []

    def test_a_relative_case_directory_is_made_absolute(self, home, parallelCase,
                                                        extractRecorder, monkeypatch):
        monkeypatch.chdir(parallelCase.parent)
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=parallelCase.name, times=0)
        assert all(os.path.isabs(call["path"]) for call in extractRecorder)


@pytest.mark.unit
class TestReadFieldAsDataFrameSerial:
    def test_an_explicit_time_is_read_from_the_case_root(self, home, serialCase,
                                                         extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase), times=0,
                                  readParallel=False)
        assert [call["path"] for call in extractRecorder] == \
            [os.path.join(str(serialCase), "0", "p")]

    def test_no_processor_column_is_added(self, home, serialCase, extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase), times=0,
                                  readParallel=False)
        assert "processor" not in extractRecorder[0]

    def test_several_times_become_several_reads(self, home, serialCase, extractRecorder):
        frame = home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase),
                                          times=[0, 100], readParallel=False)
        assert len(extractRecorder) == 2
        assert len(frame) == 2

    @pytest.mark.xfail(
        strict=True,
        reason="B291: with times=None the serial branch tests "
               "`os.path.isdir(x)` on the bare entry name rather than on "
               "os.path.join(finalCasePath, x) -- the parallel branch six "
               "lines above joins it correctly -- so no time directory is "
               "ever recognised unless the process cwd happens to be the case "
               "itself, and the trailing pandas.concat([]) raises 'No objects "
               "to concatenate'.  See the consolidated findings issue.",
    )
    def test_times_none_should_discover_the_time_directories(self, home, serialCase,
                                                             extractRecorder):
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase), times=None,
                                  readParallel=False)
        assert {os.path.basename(os.path.dirname(call["path"])) for call in extractRecorder} == \
            {"0", "100"}

    def test_times_none_currently_finds_nothing_to_concatenate(self, home, serialCase,
                                                               extractRecorder):
        """Characterisation of B291."""
        with pytest.raises(ValueError, match="No objects to concatenate"):
            home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase),
                                      times=None, readParallel=False)
        assert extractRecorder == []

    def test_times_none_only_works_from_inside_the_case(self, home, serialCase,
                                                        extractRecorder, monkeypatch):
        """Characterisation of B291: the unqualified isdir() resolves
        against the cwd, so the very same call succeeds there."""
        monkeypatch.chdir(serialCase)
        home.readFieldAsDataFrame(fieldName="p", caseDirectory=str(serialCase), times=None,
                                  readParallel=False)
        assert {os.path.basename(os.path.dirname(call["path"])) for call in extractRecorder} == \
            {"0", "100"}
