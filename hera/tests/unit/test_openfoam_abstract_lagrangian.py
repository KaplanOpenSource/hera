"""openFoam/lagrangian/abstractLagrangianSolver.py:
``absractStochasticLagrangianSolver_toolkitExtension`` and its ``analysis``
companion -- the dispersion-case builder, the case-result reader/cacher and
the full-mesh concentration binner.

Already covered elsewhere and deliberately not repeated here:
``robustOpenFOAMFileValuesParser`` / ``readEulerianConcentration``
(``test_openfoam_lagrangian_parsers.py``), the ``makeSource_*`` generators,
``sourcesTypeList`` and ``writeParticlePositionFile``
(``test_openfoam_lagrangian_sources.py``),
``analysis.calcConcentrationPointWise``
(``test_openfoam_lagrangian_analysis.py``) and the concrete solver's
``createDispersionFlowField`` override
(``test_openfoam_stochasticlagrangiansolver.py``).

Construction
------------
The extension is built with ``Cls.__new__(Cls)`` and only the attributes the
methods under test read are set by hand.  This is the established technique
in this suite for this hierarchy: ``__init__`` opens a
``dask.distributed.Client()``, and although the unit conftest stubs
``dask.distributed`` with a MagicMock (so no scheduler actually starts), a
MagicMock client cannot be asserted against without fabricating results.
``ext.toolkit`` is therefore a *real* ``OFToolkit`` from
``unit_toolkit_factory``, so every DB assertion below is an assertion about
what really landed in mongomock, and ``ext.daskClient`` is a hand-written
recorder whose ``map`` returns real ``dask.delayed`` objects -- only the
scheduler is replaced, so ``from_delayed`` builds a genuine dask DataFrame
and the recorded argument list is the code's own.  ``__init__`` itself is
covered separately, with ``Client`` replaced on the module by a hand-written
stand-in class, so that no scheduler can start even where a real
``dask.distributed`` is importable.

Covered deeply (path assembly, branch selection, what was written, what
reached mongomock, what was raised):
    ``__init__``, ``_toTimeFormat``, ``_locateOriginalFlowCase``,
    ``_detectParallelAndTimesteps``, ``_buildTimeMapping``,
    ``_checkAndRegisterDispersionInDB``, ``_prepareDispersionDirectory``,
    ``_copyOrLinkTimestep``, ``_finalizeDispersionInDB``,
    ``createDispersionFlowField`` (the ``useDBSupport=False`` orchestration,
    end to end on a real case tree), ``_resolveDirectoryConflict``,
    ``_syncDispersionToDB`` (its add branch; its update branch runs into
    B330), ``createAndLinkDispersionCaseDirectory``,
    ``_resolveCaseDescriptorName``, ``_checkCache``,
    ``_locateCaseDirectory``, ``_discoverTimesteps``,
    ``_loadCaseDataViaDask``, ``_saveToCacheParquet``,
    ``_saveToCacheNetCDF``, ``getCaseResults``,
    ``getOriginalFlowFieldExtent``, ``analysis.__init__``,
    ``analysis.calcConcentrationTimeStepFullMesh``,
    ``analysis._removeConcentrationCache``,
    ``analysis._createConcentrationCacheDoc``,
    ``analysis._processConcentrationPartition``,
    ``analysis._padRemainingTimesteps``,
    ``analysis._writeConcentrationPartition``.

Covered shallowly, with the reason:
    ``getDispersionDocument`` and ``getOriginalFlowFieldMesh`` -- each body
    is a single forwarding call, so what is asserted is the arguments that
    crossed the seam and the object that came back, nothing more.
    ``_resolveDispersionFlowField``, ``_checkDBConsistency`` and
    ``createDispersionCaseDirectory`` -- their first DB statement calls a
    method that does not exist anywhere in hera (B330), so only the
    argument validation that precedes it and the resulting
    ``AttributeError`` can be observed.
    ``getDispersionFlowDocument`` / ``getOriginalFlowDocument`` -- the two
    reachable branches are asserted on the forwarded query; the third is
    B327.
    ``getMassFromLog``, ``getCaseConcentrationsEulerian``,
    ``getOriginalFlowFieldExtentAsDict``,
    ``analysis.calcConcentrationFieldFullMesh``,
    ``analysis.calcDocumentConcentrationPointWise`` and
    ``analysis.getConcentrationField`` -- each dies on a defect pinned
    below, so only the failure and the pre-failure behaviour are asserted.

Not covered:
    ``analysis.calcConcentrationPointWise`` -- see above.
    ``readLagrangianRecord`` -- it needs five real OpenFOAM per-timestep
    field files at once and its own parser is already covered in
    ``test_openfoam_lagrangian_parsers.py``; the empty-input path it would
    exercise here proves nothing about the field assembly.

Bugs pinned below (each with a strict xfail for the intended behaviour and a
passing characterisation of what happens today):

* B324 ``getCaseConcentrationsEulerian`` passes five positional arguments
  to the four-parameter ``_loadCaseDataViaDask``.
* B325 ``getOriginalFlowFieldExtentAsDict`` reads ``zmax`` off the *y*
  column.
* B326 ``analysis.calcConcentrationFieldFullMesh``'s first statement calls
  a method that lives on the other class.
* B327 ``getDispersionFlowDocument`` / ``getOriginalFlowDocument`` build a
  type-error message and never raise it.
* B328 ``getMassFromLog`` reads ``self._datalayer``, which this class never
  defines.
* B329 ``analysis.calcDocumentConcentrationPointWise`` /
  ``analysis.getConcentrationField`` call the cache API on the wrong object.
* B330 four call sites use ``toolkit.getCaseListDocumentFromDB``, a method
  that does not exist anywhere in hera.
* B331 ``_checkDBConsistency`` calls ``toolkit.compareWorkflowsObj``; the
  real name has no ``s``.
* B332 the module's star-import shadows the builtin ``min`` with the unum
  minute unit, killing ``_buildTimeMapping``'s nearest-timestep selection.
"""
import glob
import os

import dask.dataframe
import numpy
import pandas
import pytest
import xarray
from dask.delayed import delayed

from hera import toolkitHome
from hera.datalayer import datatypes
from hera.simulations.openFoam import FLOWTYPE_INCOMPRESSIBLE
from hera.simulations.openFoam.OFWorkflow import workflow_StochasticLagrangianSolver
from hera.simulations.openFoam.preprocessOFObjects.OFObjectHome import OFObjectHome
from hera.simulations.openFoam.lagrangian import abstractLagrangianSolver as mod
from hera.simulations.openFoam.lagrangian.abstractLagrangianSolver import (
    absractStochasticLagrangianSolver_toolkitExtension as Extension,
    analysis,
)


# ---------------------------------------------------------------------------
# Hand-written stand-ins.  Real classes rather than MagicMocks, so that the
# production code's isinstance() dispatch is genuinely exercised.
# ---------------------------------------------------------------------------

class _RecordingDaskClient:
    """Stand-in for ``dask.distributed.Client`` that records ``map`` calls.

    ``map`` returns real ``dask.delayed`` objects, so the caller's
    ``dask.dataframe.from_delayed`` builds a genuine dask DataFrame; only the
    distributed scheduler is replaced.
    """

    def __init__(self):
        self.mapped = []

    def map(self, fn, items):
        items = list(items)
        self.mapped.append(items)
        return [delayed(fn)(item) for item in items]

    @property
    def lastMapped(self):
        return self.mapped[-1]


class _StandInWorkflow(workflow_StochasticLagrangianSolver):
    """A real subclass of the workflow class the module dispatches on.

    Every field the module under test reads is a *property* on the hermes
    base class, backed by a workflow JSON this stand-in has none of, so each
    one is overridden here rather than assigned in ``__init__``.
    """

    def __init__(self, name=None, dispersionFlowFieldName=None,
                 originalFlowFieldName=None, dispersionDuration=None):
        self._name = name
        self._dispersionFlowFieldName = dispersionFlowFieldName
        self._originalFlowFieldName = originalFlowFieldName
        self._dispersionDuration = dispersionDuration

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def dispersionFlowFieldName(self):
        return self._dispersionFlowFieldName

    @dispersionFlowFieldName.setter
    def dispersionFlowFieldName(self, value):
        self._dispersionFlowFieldName = value

    @property
    def originalFlowFieldName(self):
        return self._originalFlowFieldName

    @property
    def dispersionDuration(self):
        return self._dispersionDuration

    @property
    def workflowType(self):
        return "standInType"

    @property
    def json(self):
        return {"a": 1}

    @property
    def parametersJSON(self):
        return {"b": 2}


class _StandInMesh:
    """Stand-in for the object ``toolkit.getMesh`` returns."""

    def __init__(self, frame):
        self._frame = frame

    def getDataFrame(self):
        return self._frame


class _StandInField:
    """Stand-in for an OFField; records where it was asked to write."""

    def __init__(self):
        self.writes = []

    def writeToCase(self, caseDirectory, timeOrLocation):
        self.writes.append((caseDirectory, timeOrLocation))
        os.makedirs(os.path.join(caseDirectory, str(timeOrLocation)), exist_ok=True)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def toolkit(unit_toolkit_factory):
    """A real OFToolkit backed by mongomock and a per-test files directory."""
    return unit_toolkit_factory(toolkitHome.SIMULATIONS_OPENFOAM)


@pytest.fixture()
def ext(toolkit):
    """The extension without its Client-opening ``__init__`` -- see the
    module docstring."""
    obj = Extension.__new__(Extension)
    obj.toolkit = toolkit
    obj.daskClient = _RecordingDaskClient()
    return obj


@pytest.fixture()
def ana(ext):
    """The real ``analysis`` constructor, on the real toolkit."""
    return analysis(ext)


@pytest.fixture()
def bareAnalysis(ext):
    """``analysis`` without its constructor, for the pure-computation members."""
    obj = analysis.__new__(analysis)
    obj.datalayer = ext
    return obj


def _makeCase(root, timeDirs=("0", "10"), processors=None, extra=("constant", "system")):
    """Write a minimal OpenFOAM case tree and return its path."""
    root = str(root)
    os.makedirs(root, exist_ok=True)
    for sub in extra:
        os.makedirs(os.path.join(root, sub), exist_ok=True)
        with open(os.path.join(root, sub, f"{sub}File"), "w") as handle:
            handle.write(sub)
    os.makedirs(os.path.join(root, "constant", "polyMesh"), exist_ok=True)
    with open(os.path.join(root, "constant", "polyMesh", "points"), "w") as handle:
        handle.write("points")
    if processors is None:
        for timeDir in timeDirs:
            os.makedirs(os.path.join(root, timeDir), exist_ok=True)
            with open(os.path.join(root, timeDir, "U"), "w") as handle:
                handle.write(timeDir)
    else:
        # The dispersion builder copies a root-level "0" even from a
        # decomposed case, so a realistic parallel tree has one.
        os.makedirs(os.path.join(root, "0"), exist_ok=True)
        with open(os.path.join(root, "0", "U"), "w") as handle:
            handle.write("0")
        for proc in processors:
            os.makedirs(os.path.join(root, proc, "constant", "polyMesh"), exist_ok=True)
            with open(os.path.join(root, proc, "constant", "polyMesh", "points"), "w") as handle:
                handle.write("points")
            for timeDir in timeDirs:
                os.makedirs(os.path.join(root, proc, timeDir), exist_ok=True)
                with open(os.path.join(root, proc, timeDir, "U"), "w") as handle:
                    handle.write(timeDir)
    return root


def _particleFrame(times=(0.0, 1.0)):
    """A tiny Lagrangian-shaped frame: one particle per timestep."""
    return pandas.DataFrame({
        "datetime": [float(t) for t in times],
        "x": [0.25] * len(times),
        "y": [0.25] * len(times),
        "z": [0.25] * len(times),
        "mass": [2.0] * len(times),
    })


# ===========================================================================
# __init__
# ===========================================================================

@pytest.mark.unit
class TestConstructor:
    """``Client`` is replaced on the MODULE with a hand-written stand-in, so
    no scheduler can start even if a real dask.distributed is importable."""

    @pytest.fixture()
    def clientSeam(self, monkeypatch):
        built = []

        class _Client:
            def __init__(self):
                built.append(self)

        monkeypatch.setattr(mod, "Client", _Client)
        return built

    def test_the_toolkit_is_kept_as_a_back_reference(self, toolkit, clientSeam):
        assert Extension(toolkit).toolkit is toolkit

    def test_exactly_one_dask_client_is_opened(self, toolkit, clientSeam):
        instance = Extension(toolkit)
        assert len(clientSeam) == 1
        assert instance.daskClient is clientSeam[0]

    def test_the_analysis_layer_is_built_and_points_back_at_the_extension(
            self, toolkit, clientSeam):
        instance = Extension(toolkit)
        assert isinstance(instance.analysis, analysis)
        assert instance.analysis.datalayer is instance

    def test_the_presentation_layer_is_left_unwired(self, toolkit, clientSeam):
        """Characterisation: the class attribute stays None."""
        assert Extension(toolkit).presentation is None


# ===========================================================================
# _toTimeFormat
# ===========================================================================

@pytest.mark.unit
class TestToTimeFormat:
    def test_a_whole_number_string_becomes_an_int(self):
        result = Extension._toTimeFormat("10")
        assert result == 10 and isinstance(result, int)

    def test_a_trailing_zero_decimal_also_becomes_an_int(self):
        result = Extension._toTimeFormat("10.0")
        assert result == 10 and isinstance(result, int)

    def test_a_fractional_value_stays_a_float(self):
        result = Extension._toTimeFormat("0.5")
        assert result == pytest.approx(0.5) and isinstance(result, float)

    def test_a_negative_fractional_value_stays_a_float(self):
        result = Extension._toTimeFormat(-1.5)
        assert result == pytest.approx(-1.5) and isinstance(result, float)

    def test_a_non_numeric_name_raises(self):
        with pytest.raises(ValueError):
            Extension._toTimeFormat("constant")


# ===========================================================================
# _locateOriginalFlowCase
# ===========================================================================

@pytest.mark.unit
class TestLocateOriginalFlowCase:
    def test_a_source_that_is_neither_in_the_db_nor_a_directory_raises(self, ext, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist in the DB"):
            ext._locateOriginalFlowCase({"source": str(tmp_path / "nothing")})

    def test_a_directory_without_system_and_constant_is_rejected(self, ext, tmp_path):
        plain = tmp_path / "plain"
        plain.mkdir()
        with pytest.raises(ValueError, match="not a case directory"):
            ext._locateOriginalFlowCase({"source": str(plain)})

    def test_a_case_directory_is_returned_with_its_basename_as_the_group(self, ext, tmp_path):
        case = _makeCase(tmp_path / "myCase")
        assert ext._locateOriginalFlowCase({"source": case}) == (case, "myCase")

    def test_a_hermes_workflow_document_is_rebased_onto_the_workflow_name(self, ext, tmp_path):
        resource = str(tmp_path / "runs" / "someDirectory")
        ext.toolkit.addSimulationsDocument(
            resource=resource, dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="FLOW_1", groupName="FLOW"))
        casePath, group = ext._locateOriginalFlowCase({"source": "FLOW_1"})
        assert casePath == str(tmp_path / "runs" / "FLOW_1")
        assert group == "FLOW_1"

    def test_a_document_of_another_type_is_not_seen_by_this_lookup(self, ext, tmp_path):
        """Characterisation: the DB search filters on the hermesWorkflow type,
        so the method's ``type != "hermesWorkflow"`` branch is unreachable
        from here and a foreign document falls through to the disk check."""
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "runs" / "someDirectory"),
            dataFormat=datatypes.STRING, type="someOtherType",
            desc=dict(workflowName="FLOW_2", groupName="FLOW"))
        with pytest.raises(FileNotFoundError):
            ext._locateOriginalFlowCase({"source": "FLOW_2"})

    def test_two_matching_documents_raise(self, ext, tmp_path):
        for index in (1, 2):
            ext.toolkit.addSimulationsDocument(
                resource=str(tmp_path / f"r{index}"), dataFormat=datatypes.STRING,
                type=ext.toolkit.DOCTYPE_WORKFLOW,
                desc=dict(workflowName="FLOW_DUP", groupName="FLOW"))
        with pytest.raises(ValueError, match="more than one simulation"):
            ext._locateOriginalFlowCase({"source": "FLOW_DUP"})


# ===========================================================================
# _detectParallelAndTimesteps
# ===========================================================================

@pytest.mark.unit
class TestDetectParallelAndTimesteps:
    def test_a_serial_case_reports_its_own_time_directories(self, ext, tmp_path):
        case = _makeCase(tmp_path / "serial", timeDirs=("0", "5", "10"))
        isParallel, timesteps = ext._detectParallelAndTimesteps(case)
        assert isParallel is False
        assert timesteps == [0.0, 5.0, 10.0]

    def test_a_parallel_case_reads_the_times_from_processor0_only(self, ext, tmp_path):
        case = _makeCase(tmp_path / "par", timeDirs=("0", "20"),
                         processors=("processor0", "processor1"))
        isParallel, timesteps = ext._detectParallelAndTimesteps(case)
        assert isParallel is True
        assert timesteps == [0.0, 20.0]

    def test_fractional_time_directories_are_accepted(self, ext, tmp_path):
        case = _makeCase(tmp_path / "frac", timeDirs=("0", "0.5", "1.5"))
        _, timesteps = ext._detectParallelAndTimesteps(case)
        assert timesteps == [0.0, 0.5, 1.5]

    def test_non_numeric_directories_are_filtered_out(self, ext, tmp_path):
        case = _makeCase(tmp_path / "mixed", timeDirs=("0", "3"))
        _, timesteps = ext._detectParallelAndTimesteps(case)
        assert timesteps == [0.0, 3.0]

    def test_the_result_is_sorted_numerically_not_lexically(self, ext, tmp_path):
        case = _makeCase(tmp_path / "sorted", timeDirs=("100", "2", "30"))
        _, timesteps = ext._detectParallelAndTimesteps(case)
        assert timesteps == [2.0, 30.0, 100.0]


# ===========================================================================
# _buildTimeMapping
# ===========================================================================

@pytest.mark.unit
class TestBuildTimeMapping:
    def test_steady_state_without_a_timestep_freezes_the_last_one(self, ext):
        timeList, chosen = ext._buildTimeMapping(
            ext.toolkit.TIME_STEADYSTATE, [0.0, 100.0, 400.0], None, 60)
        assert chosen == 400.0
        assert timeList == [("400", "0"), ("400", "60")]


    def test_dynamic_without_a_timestep_starts_at_the_first_one(self, ext):
        timeList, chosen = ext._buildTimeMapping(
            ext.toolkit.TIME_DYNAMIC, [10.0, 20.0, 30.0], None, 20)
        assert chosen == 10.0
        assert timeList == [("10", "0"), ("20", "10"), ("30", "20")]


    def test_dynamic_extends_the_last_timestep_to_the_dispersion_duration(self, ext):
        timeList, _ = ext._buildTimeMapping(
            ext.toolkit.TIME_DYNAMIC, [0.0, 5.0], None, 100)
        assert timeList[-1] == ("5", "100")
        assert len(timeList) == 3

    def test_dynamic_does_not_extend_when_the_flow_outlasts_the_dispersion(self, ext):
        timeList, _ = ext._buildTimeMapping(
            ext.toolkit.TIME_DYNAMIC, [0.0, 5.0, 500.0], None, 100)
        assert timeList == [("0", "0"), ("5", "5"), ("500", "500")]

    def test_a_requested_timestep_is_currently_unusable(self, ext):
        """Characterisation of B332."""
        with pytest.raises(TypeError, match="not callable"):
            ext._buildTimeMapping(
                ext.toolkit.TIME_STEADYSTATE, [0.0, 100.0, 400.0], 120, 60)

    def test_an_unknown_temporal_type_raises_and_names_both_valid_ones(self, ext):
        with pytest.raises(ValueError) as caught:
            ext._buildTimeMapping("someNonsense", [0.0], None, 10)
        message = str(caught.value)
        assert ext.toolkit.TIME_STEADYSTATE in message
        assert ext.toolkit.TIME_DYNAMIC in message


# ===========================================================================
# B332: the star-import shadows the builtin min
# ===========================================================================

@pytest.mark.unit
class TestBuiltinMinIsShadowed:
    def test_the_module_global_named_min_is_a_unit_not_the_builtin(self):
        """Characterisation of B332, at its source."""
        import builtins

        shadowed = vars(mod)["min"]
        assert shadowed is not builtins.min
        assert not callable(shadowed)

    @pytest.mark.xfail(
        strict=True,
        reason="B332: the module does `from hera.utils.unitHandler import *`, "
               "and that module re-exports the unum unit table, whose entry "
               "for the minute is literally named `min`. It therefore binds a "
               "unum.Unum over the builtin `min` in this module's globals. "
               "_buildTimeMapping's only two uses of it -- "
               "TS[min(range(len(TS)), key=lambda i: abs(TS[i] - timeStep))], "
               "once in the steadyState branch and once in the dynamic one -- "
               "are exactly the code that honours an explicit `timeStep` in "
               "flowData['originalFlow'], so passing one raises "
               "TypeError: 'Unum' object is not callable and only the "
               "first/last-timestep defaults can ever be used. See the "
               "consolidated findings issue.")
    def test_steady_state_should_snap_to_the_nearest_available_timestep(self, ext):
        timeList, chosen = ext._buildTimeMapping(
            ext.toolkit.TIME_STEADYSTATE, [0.0, 100.0, 400.0], 120, 60)
        assert chosen == 100.0
        assert timeList == [("100", "0"), ("100", "60")]

    @pytest.mark.xfail(
        strict=True,
        reason="B332, second call site: the dynamic branch's nearest-timestep "
               "selection uses the same shadowed `min`. See the consolidated "
               "findings issue.")
    def test_dynamic_should_drop_timesteps_before_the_requested_start(self, ext):
        timeList, _ = ext._buildTimeMapping(
            ext.toolkit.TIME_DYNAMIC, [0.0, 10.0, 20.0], 10, 10)
        assert timeList == [("10", "0"), ("20", "10")]

    def test_the_dynamic_branch_raises_the_same_way(self, ext):
        """Characterisation of B332, second call site."""
        with pytest.raises(TypeError, match="not callable"):
            ext._buildTimeMapping(
                ext.toolkit.TIME_DYNAMIC, [0.0, 10.0, 20.0], 10, 10)

    def test_a_flow_with_an_explicit_timestep_cannot_be_prepared_at_all(
            self, ext, tmp_path):
        """Characterisation of B332, reached through the public entry point."""
        origin = _makeCase(tmp_path / "myFlow", timeDirs=("0", "400"))
        flowData = {
            "originalFlow": {
                "time": {"temporalType": ext.toolkit.TIME_STEADYSTATE},
                "timeStep": 400,
            },
            "dispersionFields": {},
        }
        with pytest.raises(TypeError, match="not callable"):
            ext.createDispersionFlowField(
                flowName="a", flowData=flowData, OriginalFlowField=origin,
                dispersionDuration=60, useDBSupport=False)


# ===========================================================================
# _checkAndRegisterDispersionInDB
# ===========================================================================

_FLOW_DATA = {"dispersionFields": {"ustar": 0.25}}
_ORIGINAL_FLOW = {"source": "FLOW_1", "time": {"temporalType": "steadyState"}}


def _dispersionDesc(name, group, flowData=None, duration=60, originalFlow=None):
    return dict(
        groupName=group,
        workflowName=name,
        flowParameters=dict(
            flowFields=(flowData or _FLOW_DATA)["dispersionFields"],
            dispersionDuration=duration,
            originalFlow=originalFlow or _ORIGINAL_FLOW,
        ),
    )


@pytest.mark.unit
class TestCheckAndRegisterDispersionInDB:
    def test_without_db_support_it_short_circuits_to_none(self, ext):
        assert ext._checkAndRegisterDispersionInDB(
            False, "NAME", "GROUP", _ORIGINAL_FLOW, _FLOW_DATA, 60, False) is None

    def test_an_empty_database_yields_none(self, ext):
        assert ext._checkAndRegisterDispersionInDB(
            True, "GROUP_DFF_a", "GROUP_DFF", _ORIGINAL_FLOW, _FLOW_DATA, 60, False) is None

    def test_an_identical_record_is_returned_for_the_caller_to_reject(self, ext, tmp_path):
        desc = _dispersionDesc("GROUP_DFF_a", "GROUP_DFF")
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "existing"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION, desc=desc)
        found = ext._checkAndRegisterDispersionInDB(
            True, "GROUP_DFF_a", "GROUP_DFF", _ORIGINAL_FLOW, _FLOW_DATA, 60, False)
        assert found is not None
        assert found.resource == str(tmp_path / "existing")

    def test_the_same_name_with_different_parameters_is_also_returned(self, ext, tmp_path):
        desc = _dispersionDesc("GROUP_DFF_a", "GROUP_DFF", duration=999)
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "other"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION, desc=desc)
        found = ext._checkAndRegisterDispersionInDB(
            True, "GROUP_DFF_a", "GROUP_DFF", _ORIGINAL_FLOW, _FLOW_DATA, 60, False)
        assert found is not None
        assert found.desc["flowParameters"]["dispersionDuration"] == 999

    def test_overwrite_deletes_the_matching_document_and_returns_none(self, ext, tmp_path):
        desc = _dispersionDesc("GROUP_DFF_a", "GROUP_DFF")
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "existing"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION, desc=desc)
        assert ext._checkAndRegisterDispersionInDB(
            True, "GROUP_DFF_a", "GROUP_DFF", _ORIGINAL_FLOW, _FLOW_DATA, 60, True) is None
        assert len(ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)) == 0

    def test_duplicate_records_are_refused_rather_than_guessed_at(self, ext, tmp_path):
        for index in (1, 2):
            ext.toolkit.addSimulationsDocument(
                resource=str(tmp_path / f"d{index}"), dataFormat=datatypes.STRING,
                type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION,
                desc=_dispersionDesc("GROUP_DFF_a", "GROUP_DFF"))
        with pytest.raises(ValueError, match="Fix manually"):
            ext._checkAndRegisterDispersionInDB(
                True, "GROUP_DFF_a", "GROUP_DFF", _ORIGINAL_FLOW, _FLOW_DATA, 60, False)


# ===========================================================================
# _copyOrLinkTimestep
# ===========================================================================

@pytest.mark.unit
class TestCopyOrLinkTimestep:
    def test_linked_data_becomes_one_symlink_per_file_in_the_timestep(self, ext, tmp_path):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",))
        dest = str(tmp_path / "dest")
        ext._copyOrLinkTimestep(origin, dest, "10", "0", False, True, True)
        link = os.path.join(dest, "0", "U")
        assert os.path.islink(link)
        assert os.path.realpath(link) == os.path.join(origin, "10", "U")

    def test_unlinked_data_is_copied_as_real_files(self, ext, tmp_path):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",))
        dest = str(tmp_path / "dest")
        ext._copyOrLinkTimestep(origin, dest, "10", "0", False, False, True)
        assert os.path.isfile(os.path.join(dest, "0", "U"))
        assert not os.path.islink(os.path.join(dest, "0", "U"))

    def test_a_parallel_timestep_lands_under_the_processor_subdirectory(self, ext, tmp_path):
        case = _makeCase(tmp_path / "par", timeDirs=("10",), processors=("processor3",))
        dest = str(tmp_path / "dest")
        ext._copyOrLinkTimestep(os.path.join(case, "processor3"), dest,
                                "10", "0", True, True, True)
        assert os.path.isfile(os.path.join(dest, "processor3", "0", "U"))

    def test_a_linked_mesh_becomes_a_symlink_to_the_original_polymesh(self, ext, tmp_path):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",))
        dest = str(tmp_path / "dest")
        ext._copyOrLinkTimestep(origin, dest, "10", "0", False, True, True)
        meshLink = os.path.join(dest, "constant", "polyMesh")
        assert os.path.islink(meshLink)
        assert os.path.realpath(meshLink) == os.path.join(origin, "constant", "polyMesh")

    def test_an_unlinked_mesh_copies_the_whole_constant_directory(self, ext, tmp_path):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",))
        dest = str(tmp_path / "dest")
        ext._copyOrLinkTimestep(origin, dest, "10", "0", False, True, False)
        target = os.path.join(dest, "constant", "polyMesh", "points")
        assert os.path.isfile(target) and not os.path.islink(target)

    def test_an_existing_destination_timestep_is_replaced_not_merged(self, ext, tmp_path):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",))
        dest = str(tmp_path / "dest")
        os.makedirs(os.path.join(dest, "0"))
        with open(os.path.join(dest, "0", "stale"), "w") as handle:
            handle.write("stale")
        ext._copyOrLinkTimestep(origin, dest, "10", "0", False, True, True)
        assert not os.path.exists(os.path.join(dest, "0", "stale"))
        assert os.path.exists(os.path.join(dest, "0", "U"))


# ===========================================================================
# _prepareDispersionDirectory
# ===========================================================================

@pytest.fixture()
def fieldSeam(monkeypatch):
    """Record every ``getEmptyFieldFromCase`` call, on the CLASS."""
    calls = []
    field = _StandInField()

    def recorder(self, **kwargs):
        calls.append(kwargs)
        return field

    monkeypatch.setattr(OFObjectHome, "getEmptyFieldFromCase", recorder)
    return calls, field


@pytest.mark.unit
class TestPrepareDispersionDirectory:
    def test_constant_system_and_zero_are_copied_from_the_original_case(
            self, ext, tmp_path, fieldSeam):
        origin = _makeCase(tmp_path / "orig", timeDirs=("0", "10"))
        dest = str(tmp_path / "dff")
        ext._prepareDispersionDirectory(dest, origin, False, [("10", "0")],
                                       True, True, {}, FLOWTYPE_INCOMPRESSIBLE)
        assert os.path.isfile(os.path.join(dest, "constant", "constantFile"))
        assert os.path.isfile(os.path.join(dest, "system", "systemFile"))
        assert os.path.isdir(os.path.join(dest, "0"))

    def test_every_destination_time_gets_a_root_level_directory(
            self, ext, tmp_path, fieldSeam):
        origin = _makeCase(tmp_path / "orig", timeDirs=("0", "10"))
        dest = str(tmp_path / "dff")
        ext._prepareDispersionDirectory(dest, origin, False,
                                       [("10", "0"), ("10", "60")],
                                       True, True, {}, FLOWTYPE_INCOMPRESSIBLE)
        assert os.path.isdir(os.path.join(dest, "0"))
        assert os.path.isdir(os.path.join(dest, "60"))

    def test_a_parallel_original_is_walked_processor_by_processor(
            self, ext, tmp_path, fieldSeam):
        origin = _makeCase(tmp_path / "orig", timeDirs=("10",),
                           processors=("processor0", "processor1"))
        dest = str(tmp_path / "dff")
        ext._prepareDispersionDirectory(dest, origin, True, [("10", "0")],
                                       True, True, {}, FLOWTYPE_INCOMPRESSIBLE)
        assert os.path.isfile(os.path.join(dest, "processor0", "0", "U"))
        assert os.path.isfile(os.path.join(dest, "processor1", "0", "U"))

    def test_each_dispersion_field_is_requested_once_per_destination_time(
            self, ext, tmp_path, fieldSeam):
        calls, field = fieldSeam
        origin = _makeCase(tmp_path / "orig", timeDirs=("0", "10"))
        dest = str(tmp_path / "dff")
        ext._prepareDispersionDirectory(
            dest, origin, False, [("10", "0"), ("10", "60")], True, True,
            {"dispersionFields": {"ustar": 0.25, "Hmix": 1000}},
            FLOWTYPE_INCOMPRESSIBLE)
        assert len(calls) == 4
        assert [c["fieldName"] for c in calls] == ["ustar", "Hmix", "ustar", "Hmix"]
        assert [c["internalValue"] for c in calls] == [0.25, 1000, 0.25, 1000]
        assert {c["caseDirectory"] for c in calls} == {dest}
        assert {c["flowType"] for c in calls} == {FLOWTYPE_INCOMPRESSIBLE}
        assert field.writes == [(dest, "0"), (dest, "0"), (dest, "60"), (dest, "60")]

    def test_no_dispersion_fields_means_no_field_requests(self, ext, tmp_path, fieldSeam):
        calls, _ = fieldSeam
        origin = _makeCase(tmp_path / "orig", timeDirs=("0", "10"))
        ext._prepareDispersionDirectory(str(tmp_path / "dff"), origin, False,
                                       [("10", "0")], True, True, {},
                                       FLOWTYPE_INCOMPRESSIBLE)
        assert calls == []


# ===========================================================================
# _finalizeDispersionInDB
# ===========================================================================

@pytest.mark.unit
class TestFinalizeDispersionInDB:
    def test_without_db_support_nothing_is_written_and_none_comes_back(self, ext, tmp_path):
        assert ext._finalizeDispersionInDB(
            False, None, "NAME", "GROUP", str(tmp_path), _ORIGINAL_FLOW,
            _FLOW_DATA, 60) is None
        assert len(ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)) == 0

    def test_a_new_flow_is_registered_and_its_directory_returned(self, ext, tmp_path):
        directory = str(tmp_path / "dff")
        result = ext._finalizeDispersionInDB(
            True, None, "GROUP_DFF_a", "GROUP_DFF", directory,
            _ORIGINAL_FLOW, _FLOW_DATA, 60)
        assert result == directory
        docs = ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)
        assert len(docs) == 1
        assert docs[0].resource == directory
        assert docs[0].dataFormat == datatypes.STRING
        assert docs[0].desc["workflowName"] == "GROUP_DFF_a"
        assert docs[0].desc["groupName"] == "GROUP_DFF"
        assert docs[0].desc["flowParameters"]["dispersionDuration"] == 60
        assert docs[0].desc["flowParameters"]["flowFields"] == {"ustar": 0.25}
        assert docs[0].desc["flowParameters"]["originalFlow"] == _ORIGINAL_FLOW

    def test_an_existing_document_is_updated_in_place_and_its_resource_returned(
            self, ext, tmp_path):
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "old"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION,
            desc=_dispersionDesc("GROUP_DFF_a", "GROUP_DFF", duration=1))
        doc = ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)[0]
        result = ext._finalizeDispersionInDB(
            True, doc, "GROUP_DFF_a", "GROUP_DFF", str(tmp_path / "ignored"),
            _ORIGINAL_FLOW, _FLOW_DATA, 60)
        assert result == str(tmp_path / "old")
        reread = ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)
        assert len(reread) == 1
        assert reread[0].desc["flowParameters"]["dispersionDuration"] == 60


# ===========================================================================
# createDispersionFlowField (orchestration)
# ===========================================================================

@pytest.mark.unit
class TestCreateDispersionFlowField:
    def test_the_whole_steady_state_tree_is_built_without_db_support(
            self, ext, tmp_path, fieldSeam, monkeypatch):
        calls, _ = fieldSeam
        origin = _makeCase(tmp_path / "myFlow", timeDirs=("0", "400"))
        monkeypatch.setattr(type(ext.toolkit), "FilesDirectory",
                            property(lambda self: str(tmp_path / "files")))
        flowData = {
            "originalFlow": {"time": {"temporalType": ext.toolkit.TIME_STEADYSTATE}},
            "dispersionFields": {"ustar": 0.25},
        }
        result = ext.createDispersionFlowField(
            flowName="a", flowData=flowData, OriginalFlowField=origin,
            dispersionDuration=60, useDBSupport=False)
        assert result is None
        built = str(tmp_path / "files" / "myFlow_DFF_a")
        assert os.path.isdir(os.path.join(built, "0"))
        assert os.path.isdir(os.path.join(built, "60"))
        assert os.path.isfile(os.path.join(built, "system", "systemFile"))
        assert [c["fieldName"] for c in calls] == ["ustar", "ustar"]

    def test_the_group_name_gets_the_dff_suffix_and_the_flow_name_appended(
            self, ext, tmp_path, fieldSeam, monkeypatch):
        origin = _makeCase(tmp_path / "baseFlow", timeDirs=("0",))
        monkeypatch.setattr(type(ext.toolkit), "FilesDirectory",
                            property(lambda self: str(tmp_path / "files")))
        flowData = {
            "originalFlow": {"time": {"temporalType": ext.toolkit.TIME_STEADYSTATE}},
            "dispersionFields": {},
        }
        ext.createDispersionFlowField(
            flowName="run7", flowData=flowData, OriginalFlowField=origin,
            dispersionDuration=10, useDBSupport=False)
        assert os.path.isdir(str(tmp_path / "files" / "baseFlow_DFF_run7"))

    def test_an_existing_registration_without_overwrite_is_refused(
            self, ext, tmp_path, fieldSeam):
        origin = _makeCase(tmp_path / "myFlow", timeDirs=("0", "400"))
        flowData = {
            "originalFlow": {"time": {"temporalType": ext.toolkit.TIME_STEADYSTATE}},
            "dispersionFields": {"ustar": 0.25},
        }
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "already"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION,
            desc=dict(workflowName="myFlow_DFF_a", groupName="myFlow_DFF"))
        with pytest.raises(FileExistsError, match="already exists"):
            ext.createDispersionFlowField(
                flowName="a", flowData=flowData, OriginalFlowField=origin,
                dispersionDuration=60, useDBSupport=True)

    def test_the_registered_document_carries_the_resolved_timestep_and_base_directory(
            self, ext, tmp_path, fieldSeam, monkeypatch):
        origin = _makeCase(tmp_path / "myFlow", timeDirs=("0", "400"))
        monkeypatch.setattr(type(ext.toolkit), "FilesDirectory",
                            property(lambda self: str(tmp_path / "files")))
        flowData = {
            "originalFlow": {"time": {"temporalType": ext.toolkit.TIME_STEADYSTATE}},
            "dispersionFields": {},
        }
        ext.createDispersionFlowField(
            flowName="a", flowData=flowData, OriginalFlowField=origin,
            dispersionDuration=60, useDBSupport=True)
        docs = ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION)
        assert len(docs) == 1
        original = docs[0].desc["flowParameters"]["originalFlow"]
        assert original["baseFlowDirectory"] == origin
        assert original["timeStep"] == 400.0
        assert original["source"] == origin


# ===========================================================================
# _resolveDirectoryConflict
# ===========================================================================

@pytest.mark.unit
class TestResolveDirectoryConflict:
    def test_a_missing_directory_is_not_a_conflict(self, ext, tmp_path):
        assert ext._resolveDirectoryConflict(str(tmp_path / "absent"), "anything", False) is None

    def test_a_consistent_existing_directory_is_left_alone(self, ext, tmp_path):
        directory = tmp_path / "case"
        directory.mkdir()
        ext._resolveDirectoryConflict(str(directory), None, False)
        assert directory.exists()

    def test_an_inconsistent_directory_without_rewrite_raises(self, ext, tmp_path):
        directory = tmp_path / "case"
        directory.mkdir()
        with pytest.raises(ValueError, match="rewrite=True"):
            ext._resolveDirectoryConflict(str(directory), "someDifference", False)
        assert directory.exists()

    def test_an_inconsistent_directory_with_rewrite_is_removed(self, ext, tmp_path):
        directory = tmp_path / "case"
        directory.mkdir()
        (directory / "f").write_text("x")
        ext._resolveDirectoryConflict(str(directory), "someDifference", True)
        assert not directory.exists()


# ===========================================================================
# _syncDispersionToDB
# ===========================================================================

@pytest.mark.unit
class TestSyncDispersionToDB:
    def test_nothing_happens_when_neither_flag_is_set(self, ext, tmp_path):
        workflow = _StandInWorkflow(name="GROUP_1")
        ext._syncDispersionToDB(workflow, str(tmp_path), "FF", False, False)
        assert len(ext.toolkit.getSimulationsDocuments(
            type=ext.toolkit.DOCTYPE_WORKFLOW)) == 0

    def test_a_new_workflow_is_added_with_group_and_id_split_from_its_name(
            self, ext, tmp_path):
        workflow = _StandInWorkflow(name="GROUP_17")
        ext._syncDispersionToDB(workflow, str(tmp_path / "case"), "THE_FLOW_FIELD",
                                True, False)
        docs = ext.toolkit.getSimulationsDocuments(type=ext.toolkit.DOCTYPE_WORKFLOW)
        assert len(docs) == 1
        assert docs[0].resource == str(tmp_path / "case")
        assert docs[0].desc["groupName"] == "GROUP"
        assert docs[0].desc["groupID"] == "17"
        assert docs[0].desc["workflowType"] == "standInType"
        assert docs[0].desc["workflow"] == {"a": 1}
        assert docs[0].desc["parameters"] == {"b": 2}

    def test_the_added_document_records_the_flow_field_name_not_the_workflow_name(
            self, ext, tmp_path):
        """Characterisation: the desc's ``workflowName`` is the *flow field*."""
        workflow = _StandInWorkflow(name="GROUP_17")
        ext._syncDispersionToDB(workflow, str(tmp_path / "case"), "THE_FLOW_FIELD",
                                True, False)
        docs = ext.toolkit.getSimulationsDocuments(type=ext.toolkit.DOCTYPE_WORKFLOW)
        assert docs[0].desc["workflowName"] == "THE_FLOW_FIELD"

    def test_a_name_without_an_underscore_cannot_be_split_into_group_and_id(
            self, ext, tmp_path):
        """Characterisation: the naming convention is mandatory."""
        workflow = _StandInWorkflow(name="NoUnderscore")
        with pytest.raises(IndexError):
            ext._syncDispersionToDB(workflow, str(tmp_path), "FF", True, False)


# ===========================================================================
# createAndLinkDispersionCaseDirectory
# ===========================================================================

@pytest.mark.unit
class TestCreateAndLinkDispersionCaseDirectory:
    def test_a_flow_directory_that_does_not_exist_is_refused(self, ext, tmp_path):
        with pytest.raises(ValueError, match="not a directory"):
            ext.createAndLinkDispersionCaseDirectory(
                str(tmp_path / "disp"), str(tmp_path / "noSuchFlow"))

    def test_constant_and_system_are_copied_from_the_flow_field(self, ext, tmp_path):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",))
        disp = str(tmp_path / "disp")
        ext.createAndLinkDispersionCaseDirectory(disp, flow)
        assert os.path.isfile(os.path.join(disp, "constant", "constantFile"))
        assert os.path.isfile(os.path.join(disp, "system", "systemFile"))

    def test_the_root_gets_a_rootcase_symlink_and_a_zero_time_directory(
            self, ext, tmp_path):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",))
        disp = str(tmp_path / "disp")
        ext.createAndLinkDispersionCaseDirectory(disp, flow)
        assert os.path.islink(os.path.join(disp, "rootCase"))
        assert os.path.realpath(os.path.join(disp, "rootCase")) == flow
        assert os.path.isdir(os.path.join(disp, "0"))

    def test_each_processor_gets_a_polymesh_link_a_rootcase_link_and_a_zero(
            self, ext, tmp_path):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",),
                         processors=("processor0", "processor1"))
        disp = str(tmp_path / "disp")
        ext.createAndLinkDispersionCaseDirectory(disp, flow)
        for proc in ("processor0", "processor1"):
            mesh = os.path.join(disp, proc, "constant", "polyMesh")
            assert os.path.islink(mesh)
            assert os.path.realpath(mesh) == os.path.join(flow, proc, "constant", "polyMesh")
            assert os.path.islink(os.path.join(disp, proc, "rootCase"))
            assert os.path.isdir(os.path.join(disp, proc, "0"))

    def test_stale_constant_and_system_directories_are_replaced(self, ext, tmp_path):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",))
        disp = tmp_path / "disp"
        (disp / "constant").mkdir(parents=True)
        (disp / "constant" / "stale").write_text("x")
        (disp / "system").mkdir()
        (disp / "system" / "stale").write_text("x")
        ext.createAndLinkDispersionCaseDirectory(str(disp), flow)
        assert not (disp / "constant" / "stale").exists()
        assert (disp / "constant" / "constantFile").exists()

    def test_it_is_idempotent_because_every_link_is_guarded(self, ext, tmp_path):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",), processors=("processor0",))
        disp = str(tmp_path / "disp")
        ext.createAndLinkDispersionCaseDirectory(disp, flow)
        ext.createAndLinkDispersionCaseDirectory(disp, flow)
        assert os.path.islink(os.path.join(disp, "processor0", "rootCase"))

    def test_relative_paths_are_made_absolute_before_anything_is_touched(
            self, ext, tmp_path, monkeypatch):
        flow = _makeCase(tmp_path / "flow", timeDirs=("0",))
        monkeypatch.chdir(tmp_path)
        ext.createAndLinkDispersionCaseDirectory("disp", "flow")
        assert os.path.realpath(os.path.join(str(tmp_path), "disp", "rootCase")) == flow


# ===========================================================================
# B330: getCaseListDocumentFromDB does not exist
# ===========================================================================

@pytest.mark.unit
class TestGetCaseListDocumentFromDBDoesNotExist:
    def test_the_method_the_module_calls_is_absent_from_the_toolkit(self, ext):
        """Characterisation of B330, at its source."""
        assert not hasattr(ext.toolkit, "getCaseListDocumentFromDB")
        assert hasattr(ext.toolkit, "getWorkflowListDocumentFromDB")

    @pytest.mark.xfail(
        strict=True,
        reason="B330: _resolveDispersionFlowField calls "
               "self.toolkit.getCaseListDocumentFromDB, a name defined nowhere "
               "in hera -- the only four occurrences in the repository are the "
               "call sites in this module. The real name is "
               "getWorkflowListDocumentFromDB. See the consolidated findings "
               "issue.")
    def test_resolving_a_registered_dispersion_flow_field_should_work(self, ext, tmp_path):
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "dff"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_OF_FLOWDISPERSION,
            desc=dict(workflowName="GROUP_DFF_a", groupName="GROUP_DFF"))
        name, directory = ext._resolveDispersionFlowField(
            _StandInWorkflow(dispersionFlowFieldName="GROUP_DFF_a"))
        assert name == "GROUP_DFF_a"

    def test_resolving_a_dispersion_flow_field_currently_raises(self, ext):
        """Characterisation of B330."""
        with pytest.raises(AttributeError, match="getCaseListDocumentFromDB"):
            ext._resolveDispersionFlowField(
                _StandInWorkflow(dispersionFlowFieldName="anything"))

    def test_checking_db_consistency_currently_raises(self, ext):
        """Characterisation of B330, at the second call site."""
        with pytest.raises(AttributeError, match="getCaseListDocumentFromDB"):
            ext._checkDBConsistency(_StandInWorkflow(name="GROUP_1"),
                                    True, False, False)


def _supplyMissingCaseListLookup(ext, monkeypatch):
    """Simulate the fix for B330 so that B331 becomes reachable.

    Both are patched on the toolkit CLASS, never on the instance.
    """
    monkeypatch.setattr(
        type(ext.toolkit), "getCaseListDocumentFromDB",
        lambda self, name, **query: self.getWorkflowListDocumentFromDB(name, **query),
        raising=False)
    monkeypatch.setattr(
        type(ext.toolkit), "getHemresWorkflowFromDocument",
        lambda self, documentList, **kwargs: _StandInWorkflow(name="GROUP_1"))


@pytest.mark.unit
class TestCompareWorkflowsObjIsMisspelled:
    def test_only_the_singular_name_exists_on_the_toolkit(self, ext):
        """Characterisation of B331, at its source."""
        assert not hasattr(ext.toolkit, "compareWorkflowsObj")
        assert hasattr(ext.toolkit, "compareWorkflowObj")

    @pytest.mark.xfail(
        strict=True,
        reason="B331: _checkDBConsistency calls "
               "self.toolkit.compareWorkflowsObj([...]); the method on "
               "hermesWorkflowToolkit is compareWorkflowObj, without the 's'. "
               "Even with B330's getCaseListDocumentFromDB supplied, the "
               "same-name branch cannot run. See the consolidated findings "
               "issue.")
    def test_the_same_name_branch_should_compare_the_two_workflows(
            self, ext, tmp_path, monkeypatch):
        _supplyMissingCaseListLookup(ext, monkeypatch)
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "c"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="GROUP_1", groupName="GROUP"))
        workflow, updateWorkflow, needDBAdd, inconsistency = ext._checkDBConsistency(
            _StandInWorkflow(name="GROUP_1"), True, False, False)
        assert workflow.name == "GROUP_1"

    def test_the_same_name_branch_currently_raises_on_the_misspelling(
            self, ext, tmp_path, monkeypatch):
        """Characterisation of B331."""
        _supplyMissingCaseListLookup(ext, monkeypatch)
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "c"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="GROUP_1", groupName="GROUP"))
        with pytest.raises(AttributeError, match="compareWorkflowsObj"):
            ext._checkDBConsistency(_StandInWorkflow(name="GROUP_1"),
                                    True, False, False)


@pytest.mark.unit
class TestCreateDispersionCaseDirectoryValidation:
    """Only the argument validation ahead of B330 is observable."""

    def test_updatedb_and_exportfromdb_are_mutually_exclusive(self, ext):
        with pytest.raises(ValueError, match="both updateDB and exportFromDB"):
            ext.createDispersionCaseDirectory(_StandInWorkflow(name="GROUP_1"),
                                              updateDB=True, exportFromDB=True)

    def test_a_workflow_without_a_name_is_refused(self, ext):
        with pytest.raises(ValueError, match="name property"):
            ext.createDispersionCaseDirectory(_StandInWorkflow(name=None),
                                              updateDB=False)

    def test_a_valid_call_still_dies_inside_bzzz_g(self, ext, unit_files_directory):
        """Characterisation of B330, reached through the public entry point."""
        with pytest.raises(AttributeError, match="getCaseListDocumentFromDB"):
            ext.createDispersionCaseDirectory(_StandInWorkflow(name="GROUP_1"),
                                              updateDB=False)


# ===========================================================================
# getMassFromLog -- B328
# ===========================================================================

_SOLVER_LOG = "\n".join([
    "Exec   : StochasticLagrangianSolver",
    "Date   : today",
    "Case   : /somewhere",
    "Time = 1",
    "Solving 3-D cloud kinematicCloud",
    "kinematicCloud",
    "sourceA:",
    "    parcels added = 10",
    "    mass added = 0.5",
    "Parcel fate: sourceA",
    "    escape = 1, 0.1",
    "    stick  = 2, 0.2",
    "ExecutionTime = 1 s",
    "End",
    "",
])


@pytest.mark.unit
class TestGetMassFromLog:
    def test_the_attribute_the_parser_reads_is_defined_nowhere_in_the_class(self):
        """Characterisation of B328, at its source."""
        assert not hasattr(Extension, "_datalayer")
        assert hasattr(Extension, "toolkit")

    @pytest.mark.xfail(
        strict=True,
        reason="B328: getMassFromLog reads self._datalayer.cloudName, but it "
               "is defined on absractStochasticLagrangianSolver_toolkitExtension, "
               "whose attributes are toolkit/analysis/presentation/daskClient -- "
               "no class in the hierarchy defines _datalayer. (The near-identical "
               "copy on OFLSMToolkit.Analysis does have one; this one is a copy "
               "that lost its owner.) Every call raises AttributeError before a "
               "single line of the log is parsed. See the consolidated findings "
               "issue.")
    def test_it_should_parse_release_and_parcel_fate_records(self, ext, tmp_path):
        logFile = tmp_path / "solver.log"
        logFile.write_text(_SOLVER_LOG)
        frame = ext.getMassFromLog(str(logFile))
        assert sorted(frame["action"].unique()) == ["escape", "release", "stick"]

    def test_it_currently_raises_before_reading_any_record(self, ext, tmp_path):
        """Characterisation of B328."""
        logFile = tmp_path / "solver.log"
        logFile.write_text(_SOLVER_LOG)
        with pytest.raises(AttributeError, match="_datalayer"):
            ext.getMassFromLog(str(logFile))

    def test_even_a_log_with_no_timestep_block_cannot_be_parsed(self, ext, tmp_path):
        """Characterisation of B328: the frame it builds at the end reads
        the same missing attribute, so no input reaches a return."""
        logFile = tmp_path / "empty.log"
        logFile.write_text("Exec   : StochasticLagrangianSolver\na\nb\nEnd\n")
        with pytest.raises(AttributeError):
            ext.getMassFromLog(str(logFile))


# ===========================================================================
# _resolveCaseDescriptorName
# ===========================================================================

@pytest.mark.unit
class TestResolveCaseDescriptorName:
    def test_a_string_is_returned_unchanged(self, ext):
        assert ext._resolveCaseDescriptorName("SOME_CASE") == "SOME_CASE"

    def test_a_document_yields_its_workflow_name(self, ext, tmp_path):
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "c"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="FROM_DESC", groupName="G"))
        doc = ext.toolkit.getSimulationsDocuments(type=ext.toolkit.DOCTYPE_WORKFLOW)[0]
        assert ext._resolveCaseDescriptorName(doc) == "FROM_DESC"

    def test_anything_else_is_refused(self, ext):
        with pytest.raises(ValueError, match="string or MetadataFrame"):
            ext._resolveCaseDescriptorName(17)


# ===========================================================================
# _checkCache
# ===========================================================================

def _writeParquet(path, frame=None):
    frame = _particleFrame() if frame is None else frame
    dask.dataframe.from_pandas(frame, npartitions=1).to_parquet(str(path))
    return str(path)


@pytest.mark.unit
class TestCheckCache:
    def test_an_empty_cache_reports_a_miss(self, ext):
        assert ext._checkCache("NO_SUCH", ext.DOCTYPE_LAGRANGIAN_CACHE, False) == (None, None)

    def test_duplicate_cache_entries_are_refused(self, ext, tmp_path):
        for index in (1, 2):
            ext.toolkit.addCacheDocument(
                resource=str(tmp_path / f"c{index}"), dataFormat=datatypes.PARQUET,
                type=ext.DOCTYPE_LAGRANGIAN_CACHE,
                desc=dict(workflowName="DUP", cloudName="kinematicCloud"))
        with pytest.raises(ValueError, match="Remove duplicates"):
            ext._checkCache("DUP", ext.DOCTYPE_LAGRANGIAN_CACHE, False)

    def test_a_warm_cache_returns_the_stored_data(self, ext, tmp_path):
        resource = _writeParquet(tmp_path / "cached.parquet")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="WARM", cloudName="kinematicCloud"))
        doc, data = ext._checkCache("WARM", ext.DOCTYPE_LAGRANGIAN_CACHE, False)
        assert doc is not None
        assert data is not None
        assert len(data.compute()) == 2

    def test_overwrite_removes_a_cached_directory_and_reports_a_miss(self, ext, tmp_path):
        resource = _writeParquet(tmp_path / "cached.parquet")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="HOT", cloudName="kinematicCloud"))
        doc, data = ext._checkCache("HOT", ext.DOCTYPE_LAGRANGIAN_CACHE, True)
        assert data is None
        assert doc is not None
        assert not os.path.exists(resource)

    def test_overwrite_removes_a_cached_single_file_too(self, ext, tmp_path):
        resource = str(tmp_path / "cached.nc")
        with open(resource, "w") as handle:
            handle.write("x")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.STRING,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="HOTFILE", cloudName="kinematicCloud"))
        doc, data = ext._checkCache("HOTFILE", ext.DOCTYPE_LAGRANGIAN_CACHE, True)
        assert data is None and doc is not None
        assert not os.path.exists(resource)

    def test_the_document_is_kept_by_overwrite_so_the_path_can_be_reused(self, ext, tmp_path):
        resource = _writeParquet(tmp_path / "cached.parquet")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="KEEP", cloudName="kinematicCloud"))
        ext._checkCache("KEEP", ext.DOCTYPE_LAGRANGIAN_CACHE, True)
        assert len(ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)) == 1

    def test_a_stale_document_whose_file_vanished_is_deleted(self, ext, tmp_path):
        ext.toolkit.addCacheDocument(
            resource=str(tmp_path / "gone.parquet"), dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="STALE", cloudName="kinematicCloud"))
        assert ext._checkCache("STALE", ext.DOCTYPE_LAGRANGIAN_CACHE, False) == (None, None)
        assert len(ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)) == 0


# ===========================================================================
# _locateCaseDirectory
# ===========================================================================

@pytest.mark.unit
class TestLocateCaseDirectory:
    def test_a_case_not_in_the_db_falls_back_to_the_filesystem(self, ext, tmp_path):
        case = _makeCase(tmp_path / "onDisk")
        assert ext._locateCaseDirectory(case, case) == (case, case)

    def test_a_relative_path_is_resolved_against_the_working_directory(
            self, ext, tmp_path, monkeypatch):
        _makeCase(tmp_path / "onDisk")
        monkeypatch.chdir(tmp_path)
        path, name = ext._locateCaseDirectory("onDisk", "onDisk")
        assert path == str(tmp_path / "onDisk")
        assert name == "onDisk"

    def test_a_name_that_is_neither_in_the_db_nor_on_disk_raises(self, ext, tmp_path):
        with pytest.raises(FileNotFoundError, match="not in the DB"):
            ext._locateCaseDirectory("NOPE", "NOPE")

    def test_a_hermes_workflow_document_is_rebased_onto_the_workflow_name(
            self, ext, tmp_path):
        ext.toolkit.addSimulationsDocument(
            resource=str(tmp_path / "runs" / "stale"), dataFormat=datatypes.STRING,
            type=ext.toolkit.DOCTYPE_WORKFLOW,
            desc=dict(workflowName="DISP_1", groupName="DISP"))
        path, name = ext._locateCaseDirectory("DISP_1", "DISP_1")
        assert path == str(tmp_path / "runs" / "DISP_1")
        assert name == "DISP_1"


# ===========================================================================
# _discoverTimesteps
# ===========================================================================

@pytest.mark.unit
class TestDiscoverTimesteps:
    def test_time_directories_come_back_sorted_numerically(self, ext, tmp_path):
        for name in ("0", "2", "10", "100"):
            (tmp_path / name).mkdir()
        assert ext._discoverTimesteps(str(tmp_path)) == ["0", "2", "10", "100"]

    def test_the_openfoam_bookkeeping_directories_are_skipped(self, ext, tmp_path):
        for name in ("0", "constant", "system", "rootCase", "VTK", "processor0"):
            (tmp_path / name).mkdir()
        assert ext._discoverTimesteps(str(tmp_path)) == ["0"]

    def test_plain_files_are_not_mistaken_for_time_directories(self, ext, tmp_path):
        (tmp_path / "0").mkdir()
        (tmp_path / "5").write_text("not a directory")
        assert ext._discoverTimesteps(str(tmp_path)) == ["0"]

    def test_a_processor_subdirectory_is_scanned_when_named(self, ext, tmp_path):
        (tmp_path / "processor2").mkdir()
        for name in ("0", "7"):
            (tmp_path / "processor2" / name).mkdir()
        assert ext._discoverTimesteps(str(tmp_path), "processor2") == ["0", "7"]

    def test_fractional_time_directories_are_dropped_here(self, ext, tmp_path):
        """Characterisation: ``isdigit`` rejects ``"0.5"``, unlike the sibling
        ``_detectParallelAndTimesteps`` in the same class, which accepts it."""
        for name in ("0", "0.5", "1"):
            (tmp_path / name).mkdir()
        assert ext._discoverTimesteps(str(tmp_path)) == ["0", "1"]


# ===========================================================================
# _loadCaseDataViaDask
# ===========================================================================

def _stubLoader(timeName):
    """A loader that stamps the name it was handed into the frame."""
    return pandas.DataFrame({"datetime": [float(os.path.basename(timeName))],
                             "x": [0.0], "name": [timeName]})


@pytest.mark.unit
class TestLoadCaseDataViaDask:
    def test_a_serial_case_maps_the_loader_over_its_own_time_directories(
            self, ext, tmp_path):
        for name in ("0", "1", "2"):
            (tmp_path / name).mkdir()
        result = ext._loadCaseDataViaDask(str(tmp_path), _stubLoader, None, False)
        assert ext.daskClient.lastMapped == ["0", "1", "2"]
        assert isinstance(result, dask.dataframe.DataFrame)

    def test_a_parallel_case_maps_over_every_processor_time_pair(self, ext, tmp_path):
        case = _makeCase(tmp_path / "par", timeDirs=("0", "1"),
                         processors=("processor0", "processor1"))
        ext._loadCaseDataViaDask(case, _stubLoader, None, False)
        assert sorted(ext.daskClient.lastMapped) == [
            os.path.join("processor0", "0"), os.path.join("processor0", "1"),
            os.path.join("processor1", "0"), os.path.join("processor1", "1"),
        ]

    def test_forcing_a_single_processor_reads_the_root_time_directories(
            self, ext, tmp_path):
        case = _makeCase(tmp_path / "par", timeDirs=("0", "1"),
                         processors=("processor0",))
        os.makedirs(os.path.join(case, "5"))
        ext._loadCaseDataViaDask(case, _stubLoader, None, True)
        assert ext.daskClient.lastMapped == ["0", "5"]

    def test_an_explicit_time_list_is_used_verbatim_in_the_serial_case(self, ext, tmp_path):
        ext._loadCaseDataViaDask(str(tmp_path), _stubLoader, ["3", "9"], False)
        assert ext.daskClient.lastMapped == ["3", "9"]

    def test_an_explicit_time_list_is_crossed_with_the_processors(self, ext, tmp_path):
        case = _makeCase(tmp_path / "par", timeDirs=("0",), processors=("processor0",))
        ext._loadCaseDataViaDask(case, _stubLoader, ["4"], False)
        assert ext.daskClient.lastMapped == [os.path.join("processor0", "4")]

    def test_the_frame_that_comes_back_contains_what_the_loader_produced(
            self, ext, tmp_path):
        for name in ("0", "1"):
            (tmp_path / name).mkdir()
        result = ext._loadCaseDataViaDask(str(tmp_path), _stubLoader, None, False)
        assert sorted(result.compute()["name"]) == ["0", "1"]


# ===========================================================================
# _saveToCacheParquet
# ===========================================================================

@pytest.mark.unit
class TestSaveToCacheParquet:
    def test_a_cold_cache_builds_the_path_from_the_files_directory(self, ext):
        data = dask.dataframe.from_pandas(_particleFrame(), npartitions=1)
        result = ext._saveToCacheParquet(
            data, None, "MY_CASE", "MY_CASE", "kinematicCloud",
            ext.DOCTYPE_LAGRANGIAN_CACHE, datatypes.PARQUET)
        expected = os.path.join(ext.toolkit.filesDirectory, "cachedLagrangianData",
                                "MY_CASE", "kinematicCloud.parquet")
        assert os.path.exists(expected)
        assert isinstance(result, dask.dataframe.DataFrame)

    def test_a_cold_cache_registers_exactly_one_document(self, ext):
        data = dask.dataframe.from_pandas(_particleFrame(), npartitions=1)
        ext._saveToCacheParquet(data, None, "MY_CASE", "MY_CASE", "myCloud",
                                ext.DOCTYPE_LAGRANGIAN_CACHE, datatypes.PARQUET)
        docs = ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)
        assert len(docs) == 1
        assert docs[0].desc["workflowName"] == "MY_CASE"
        assert docs[0].desc["cloudName"] == "myCloud"
        assert docs[0].dataFormat == datatypes.PARQUET
        assert docs[0].resource.endswith(os.path.join("MY_CASE", "myCloud.parquet"))

    def test_the_written_data_is_indexed_by_datetime(self, ext):
        data = dask.dataframe.from_pandas(_particleFrame(times=(0.0, 1.0, 2.0)),
                                          npartitions=1)
        result = ext._saveToCacheParquet(
            data, None, "MY_CASE", "MY_CASE", "kinematicCloud",
            ext.DOCTYPE_LAGRANGIAN_CACHE, datatypes.PARQUET)
        computed = result.compute()
        assert computed.index.name == "datetime"
        assert sorted(computed.index) == [0.0, 1.0, 2.0]

    def test_a_warm_document_is_reused_without_registering_another(self, ext, tmp_path):
        ext.toolkit.addCacheDocument(
            resource=str(tmp_path / "existing.parquet"), dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="MY_CASE", cloudName="kinematicCloud"))
        doc = ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)[0]
        data = dask.dataframe.from_pandas(_particleFrame(), npartitions=1)
        ext._saveToCacheParquet(data, doc, "MY_CASE", "MY_CASE", "kinematicCloud",
                                ext.DOCTYPE_LAGRANGIAN_CACHE, datatypes.PARQUET)
        assert os.path.exists(str(tmp_path / "existing.parquet"))
        assert len(ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)) == 1


# ===========================================================================
# _saveToCacheNetCDF
# ===========================================================================

def _eulerianFrame():
    return pandas.DataFrame({
        "datetime": [0.0, 0.0, 1.0],
        "x": [0.0, 0.0, 0.0],
        "y": [0.0, 0.0, 0.0],
        "z": [0.0, 0.0, 0.0],
        "C": [1.0, 2.0, 4.0],
    })


@pytest.mark.unit
class TestSaveToCacheNetCDF:
    def test_a_cold_cache_writes_a_netcdf_named_after_the_cloud(self, ext):
        data = dask.dataframe.from_pandas(_eulerianFrame(), npartitions=1)
        ext._saveToCacheNetCDF(data, None, "MY_CASE", "MY_CASE", "kinematicCloud",
                               ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        expected = os.path.join(ext.toolkit.filesDirectory, "cachedLagrangianData",
                                "MY_CASE", "kinematicCloudConcentrationEulerian.nc")
        assert os.path.exists(expected)

    def test_the_document_records_the_netcdf_format(self, ext):
        data = dask.dataframe.from_pandas(_eulerianFrame(), npartitions=1)
        ext._saveToCacheNetCDF(data, None, "MY_CASE", "MY_CASE", "myCloud",
                               ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        docs = ext.toolkit.getCacheDocuments(
            type=ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        assert len(docs) == 1
        assert docs[0].dataFormat == datatypes.NETCDF_XARRAY
        assert docs[0].desc["workflowName"] == "MY_CASE"
        assert docs[0].desc["cloudName"] == "myCloud"

    def test_rows_sharing_a_cell_and_time_are_summed(self, ext):
        data = dask.dataframe.from_pandas(_eulerianFrame(), npartitions=1)
        result = ext._saveToCacheNetCDF(
            data, None, "MY_CASE", "MY_CASE", "kinematicCloud",
            ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        assert result["C"].sel(datetime=0.0).item() == pytest.approx(3.0)
        assert result["C"].sel(datetime=1.0).item() == pytest.approx(4.0)

    def test_the_result_is_indexed_by_datetime_x_y_and_z(self, ext):
        data = dask.dataframe.from_pandas(_eulerianFrame(), npartitions=1)
        result = ext._saveToCacheNetCDF(
            data, None, "MY_CASE", "MY_CASE", "kinematicCloud",
            ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        assert set(result.dims) == {"datetime", "x", "y", "z"}

    def test_a_warm_document_is_written_to_its_own_resource(self, ext, tmp_path):
        resource = str(tmp_path / "warm.nc")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.NETCDF_XARRAY,
            type=ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE,
            desc=dict(workflowName="MY_CASE", cloudName="kinematicCloud"))
        doc = ext.toolkit.getCacheDocuments(
            type=ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)[0]
        data = dask.dataframe.from_pandas(_eulerianFrame(), npartitions=1)
        ext._saveToCacheNetCDF(data, doc, "MY_CASE", "MY_CASE", "kinematicCloud",
                               ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)
        assert os.path.exists(resource)
        assert len(ext.toolkit.getCacheDocuments(
            type=ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE)) == 1


# ===========================================================================
# getCaseResults
# ===========================================================================

@pytest.fixture()
def lagrangianLoaderSeam(monkeypatch):
    """Replace the module-level record reader; record the kwargs it was given."""
    calls = []

    def recorder(timeName, **kwargs):
        calls.append(dict(timeName=timeName, **kwargs))
        return pandas.DataFrame({
            "datetime": [float(os.path.basename(timeName))],
            "x": [0.25], "y": [0.25], "z": [0.25], "mass": [2.0],
        })

    monkeypatch.setattr(mod, "readLagrangianRecord", recorder)
    return calls


@pytest.mark.unit
class TestGetCaseResults:
    def test_the_loader_receives_the_resolved_case_path_and_every_flag(
            self, ext, tmp_path, lagrangianLoaderSeam):
        case = _makeCase(tmp_path / "case", timeDirs=("0", "1"))
        result = ext.getCaseResults(case, withVelocity=False, withReleaseTimes=True,
                                    withMass=False, cloudName="myCloud", cache=False)
        result.compute()
        assert len(lagrangianLoaderSeam) > 0
        first = lagrangianLoaderSeam[0]
        assert first["casePath"] == case
        assert first["withVelocity"] is False
        assert first["withReleaseTimes"] is True
        assert first["withMass"] is False
        assert first["cloudName"] == "myCloud"

    def test_the_timesteps_discovered_on_disk_are_the_ones_loaded(
            self, ext, tmp_path, lagrangianLoaderSeam):
        case = _makeCase(tmp_path / "case", timeDirs=("0", "1", "2"))
        ext.getCaseResults(case, cache=False)
        assert ext.daskClient.lastMapped == ["0", "1", "2"]

    def test_an_explicit_time_list_overrides_discovery(
            self, ext, tmp_path, lagrangianLoaderSeam):
        case = _makeCase(tmp_path / "case", timeDirs=("0", "1", "2"))
        ext.getCaseResults(case, timeList=["1"], cache=False)
        assert ext.daskClient.lastMapped == ["1"]

    def test_caching_writes_a_parquet_and_registers_a_cache_document(
            self, ext, tmp_path, lagrangianLoaderSeam):
        case = _makeCase(tmp_path / "case", timeDirs=("0", "1"))
        ext.getCaseResults(case, cache=True, cloudName="myCloud")
        docs = ext.toolkit.getCacheDocuments(type=ext.DOCTYPE_LAGRANGIAN_CACHE)
        assert len(docs) == 1
        assert docs[0].desc["cloudName"] == "myCloud"
        assert os.path.exists(docs[0].resource)

    def test_a_warm_cache_is_returned_without_touching_the_loader(
            self, ext, tmp_path, lagrangianLoaderSeam):
        resource = _writeParquet(tmp_path / "warm.parquet")
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.PARQUET,
            type=ext.DOCTYPE_LAGRANGIAN_CACHE,
            desc=dict(workflowName="WARM_CASE", cloudName="kinematicCloud"))
        result = ext.getCaseResults("WARM_CASE")
        assert len(result.compute()) == 2
        assert lagrangianLoaderSeam == []

    def test_a_case_that_is_neither_cached_nor_on_disk_raises(
            self, ext, lagrangianLoaderSeam):
        with pytest.raises(FileNotFoundError):
            ext.getCaseResults("NO_SUCH_CASE", cache=False)


# ===========================================================================
# getDispersionDocument / getDispersionFlowDocument / getOriginalFlowDocument
# ===========================================================================

@pytest.fixture()
def workflowDocumentSeam(monkeypatch):
    """Record every ``getWorkflowDocumentFromDB`` call, on the toolkit CLASS."""
    calls = []
    documents = []

    def recorder(self, nameOrWorkflowFileOrJSONOrResource, **kwargs):
        calls.append(dict(name=nameOrWorkflowFileOrJSONOrResource, **kwargs))
        return list(documents)

    def install(toolkit, docs):
        documents[:] = docs
        monkeypatch.setattr(type(toolkit), "getWorkflowDocumentFromDB", recorder)
        return calls

    return install


@pytest.mark.unit
class TestGetDispersionDocument:
    """Delegation only: the body is a single forwarding call."""

    def test_the_name_is_forwarded_verbatim_and_the_result_returned(
            self, ext, workflowDocumentSeam):
        sentinel = object()
        calls = workflowDocumentSeam(ext.toolkit, [sentinel])
        assert ext.getDispersionDocument("SOME_NAME") == [sentinel]
        assert calls == [{"name": "SOME_NAME"}]


@pytest.mark.unit
class TestGetDispersionFlowDocument:
    def test_a_workflow_object_supplies_its_own_flow_field_name(
            self, ext, workflowDocumentSeam):
        sentinel = object()
        calls = workflowDocumentSeam(ext.toolkit, [sentinel])
        result = ext.getDispersionFlowDocument(
            _StandInWorkflow(dispersionFlowFieldName="THE_DFF"))
        assert result is sentinel
        assert calls[0]["name"] == "THE_DFF"
        assert calls[0]["doctype"] == ext.toolkit.DOCTYPE_OF_FLOWDISPERSION

    def test_nothing_found_yields_none(self, ext, workflowDocumentSeam):
        workflowDocumentSeam(ext.toolkit, [])
        assert ext.getDispersionFlowDocument(
            _StandInWorkflow(dispersionFlowFieldName="THE_DFF")) is None

    def test_a_string_is_resolved_through_the_workflow_in_the_db(
            self, ext, workflowDocumentSeam, monkeypatch):
        sentinel = object()
        calls = workflowDocumentSeam(ext.toolkit, [sentinel])
        monkeypatch.setattr(
            type(ext.toolkit), "getHermesWorkflowFromDB",
            lambda self, name, **kw: _StandInWorkflow(dispersionFlowFieldName="RESOLVED"))
        assert ext.getDispersionFlowDocument("SOME_NAME") is sentinel
        assert calls[0]["name"] == "RESOLVED"


@pytest.mark.unit
class TestGetOriginalFlowDocument:
    def test_a_workflow_object_supplies_its_original_flow_field_name(
            self, ext, workflowDocumentSeam):
        sentinel = object()
        calls = workflowDocumentSeam(ext.toolkit, [sentinel])
        result = ext.getOriginalFlowDocument(
            _StandInWorkflow(originalFlowFieldName="THE_OFF"))
        assert result is sentinel
        assert calls[0]["name"] == "THE_OFF"
        assert calls[0]["doctype"] == ext.toolkit.DOCTYPE_WORKFLOW

    def test_nothing_found_yields_none(self, ext, workflowDocumentSeam):
        workflowDocumentSeam(ext.toolkit, [])
        assert ext.getOriginalFlowDocument(
            _StandInWorkflow(originalFlowFieldName="THE_OFF")) is None


@pytest.mark.unit
class TestTypeErrorIsBuiltButNeverRaised:
    @pytest.mark.xfail(
        strict=True,
        reason="B327: getDispersionFlowDocument's else branch assigns the "
               "message to a local named err and never raises it, so execution "
               "falls through to logger.info(f'... {dffname}') with dffname "
               "bound in neither of the two taken branches. The caller gets "
               "UnboundLocalError instead of the TypeError the message "
               "describes. getOriginalFlowDocument has the identical defect. "
               "See the consolidated findings issue.")
    def test_a_bad_type_should_raise_the_error_the_message_describes(self, ext):
        with pytest.raises(TypeError):
            ext.getDispersionFlowDocument(17)

    def test_a_bad_type_currently_raises_unboundlocalerror(self, ext):
        """Characterisation of B327."""
        with pytest.raises(UnboundLocalError, match="dffname"):
            ext.getDispersionFlowDocument(17)

    def test_the_original_flow_variant_has_the_same_defect(self, ext):
        """Characterisation of B327, at the second site."""
        with pytest.raises(UnboundLocalError, match="dffname"):
            ext.getOriginalFlowDocument(17)


# ===========================================================================
# getOriginalFlowFieldMesh / getOriginalFlowFieldExtent[AsDict]
# ===========================================================================

_MESH_FRAME = pandas.DataFrame({
    "Cx": [0.0, 10.0, 5.0],
    "Cy": [-3.0, 3.0, 0.0],
    "Cz": [1.0, 7.0, 4.0],
})


@pytest.fixture()
def meshSeam(monkeypatch):
    """Seam over ``getOriginalFlowDocument`` and ``toolkit.getMesh``, both on
    their CLASS. Records what crossed."""
    calls = {}

    class _Document:
        def getData(self):
            return "theResourceData"

    def install(ext, frame=_MESH_FRAME):
        def getOriginalFlowDocument(self, name):
            calls["document"] = name
            return _Document()

        monkeypatch.setattr(Extension, "getOriginalFlowDocument",
                            getOriginalFlowDocument)

        def getMesh(self, resource, **kwargs):
            calls["getMesh"] = resource
            return _StandInMesh(frame)

        monkeypatch.setattr(type(ext.toolkit), "getMesh", getMesh)
        return calls

    return install


@pytest.mark.unit
class TestGetOriginalFlowFieldMesh:
    """Delegation only: document -> getData -> toolkit.getMesh."""

    def test_the_document_data_is_what_reaches_get_mesh(self, ext, meshSeam):
        calls = meshSeam(ext)
        mesh = ext.getOriginalFlowFieldMesh("SOME_NAME")
        assert calls["document"] == "SOME_NAME"
        assert calls["getMesh"] == "theResourceData"
        assert isinstance(mesh, _StandInMesh)


@pytest.mark.unit
class TestGetOriginalFlowFieldExtent:
    def test_it_aggregates_the_cell_centres_to_a_min_max_frame(self, ext, meshSeam):
        meshSeam(ext)
        extent = ext.getOriginalFlowFieldExtent("SOME_NAME")
        assert list(extent.index) == ["min", "max"]
        assert list(extent.columns) == ["Cx", "Cy", "Cz"]
        assert extent.Cx.loc["min"] == 0.0 and extent.Cx.loc["max"] == 10.0
        assert extent.Cz.loc["min"] == 1.0 and extent.Cz.loc["max"] == 7.0


@pytest.mark.unit
class TestGetOriginalFlowFieldExtentAsDict:
    def test_the_horizontal_bounds_and_the_floor_are_metre_quantities(self, ext, meshSeam):
        meshSeam(ext)
        limits = ext.getOriginalFlowFieldExtentAsDict("SOME_NAME")
        assert set(limits) == {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}
        assert limits["xmin"].magnitude == 0.0
        assert limits["xmax"].magnitude == 10.0
        assert limits["ymin"].magnitude == -3.0
        assert limits["ymax"].magnitude == 3.0
        assert limits["zmin"].magnitude == 1.0
        assert str(limits["xmin"].units) == "meter"

    @pytest.mark.xfail(
        strict=True,
        reason="B325: getOriginalFlowFieldExtentAsDict builds zmax from "
               "lims.Cy.loc['max'] instead of lims.Cz.loc['max']. The other "
               "five keys all read their own axis (zmin reads Cz), so this is "
               "a copy-paste slip, and it silently truncates the vertical "
               "extent of every concentration mesh built from it. See the "
               "consolidated findings issue.")
    def test_zmax_should_come_from_the_z_column(self, ext, meshSeam):
        meshSeam(ext)
        limits = ext.getOriginalFlowFieldExtentAsDict("SOME_NAME")
        assert limits["zmax"].magnitude == 7.0

    def test_zmax_currently_repeats_ymax(self, ext, meshSeam):
        """Characterisation of B325."""
        meshSeam(ext)
        limits = ext.getOriginalFlowFieldExtentAsDict("SOME_NAME")
        assert limits["zmax"].magnitude == 3.0
        assert limits["zmax"] == limits["ymax"]


# ===========================================================================
# getCaseConcentrationsEulerian -- B324
# ===========================================================================

@pytest.mark.unit
class TestGetCaseConcentrationsEulerian:
    def test_a_warm_cache_is_returned_before_the_broken_call_is_reached(
            self, ext, tmp_path):
        frame = _eulerianFrame().groupby(
            ["datetime", "x", "y", "z"]).sum()
        resource = str(tmp_path / "warm.nc")
        frame.to_xarray().to_netcdf(resource)
        ext.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.NETCDF_XARRAY,
            type=ext.DOCTYPE_CONCENTRATIONEULERIAN_CACHE,
            desc=dict(workflowName="WARM_CASE", cloudName="kinematicCloud"))
        result = ext.getCaseConcentrationsEulerian("WARM_CASE")
        assert isinstance(result, xarray.Dataset)
        assert "C" in result

    @pytest.mark.xfail(
        strict=True,
        reason="B324: getCaseConcentrationsEulerian calls "
               "self._loadCaseDataViaDask(finalCasePath, loader, timeList, "
               "forceSingleProcessor, self.daskClient) -- five positional "
               "arguments to a method declared as "
               "_loadCaseDataViaDask(self, casePath, loader, timeList, "
               "forceSingleProcessor), which reads self.daskClient itself. The "
               "sibling getCaseResults calls the same helper with four. Every "
               "cold-cache call raises TypeError. See the consolidated findings "
               "issue.")
    def test_a_cold_cache_should_read_the_case_from_disk(self, ext, tmp_path):
        case = _makeCase(tmp_path / "case", timeDirs=("0",))
        ext.getCaseConcentrationsEulerian(case, cache=False)

    def test_a_cold_cache_currently_raises_typeerror_on_the_extra_argument(
            self, ext, tmp_path):
        """Characterisation of B324."""
        case = _makeCase(tmp_path / "case", timeDirs=("0",))
        with pytest.raises(TypeError, match="_loadCaseDataViaDask"):
            ext.getCaseConcentrationsEulerian(case, cache=False)


# ===========================================================================
# analysis.__init__
# ===========================================================================

@pytest.mark.unit
class TestAnalysisConstructor:
    def test_the_datalayer_reference_is_kept(self, ext):
        instance = analysis(ext)
        assert instance.datalayer is ext

    def test_the_two_counters_are_seeded_in_the_project_configuration(self, ext):
        analysis(ext)
        config = ext.toolkit.getConfig()
        assert config["analysisFullMeshCounter"] == 0
        assert config["analysisPointWiseCounter"] == 0

    def test_the_cartesian_mesh_counter_becomes_usable(self, ext):
        instance = analysis(ext)
        assert ext.toolkit.addCounter("cartesianMeshCounter") == 1
        assert ext.toolkit.addCounter("cartesianMeshCounter") == 2


# ===========================================================================
# analysis.calcConcentrationTimeStepFullMesh
# ===========================================================================

_EXTENTS = dict(xmin=0, xmax=2, ymin=0, ymax=2, zmin=0, zmax=2)


@pytest.mark.unit
class TestCalcConcentrationTimeStepFullMesh:
    def test_an_integer_time_yields_an_all_zero_grid_stamped_with_that_time(
            self, bareAnalysis):
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            7, extents=_EXTENTS, dxdydz=1)
        assert list(result.datetime.values) == [7]
        assert float(result.sum()) == 0.0

    def test_the_grid_spans_the_extents_at_the_requested_resolution(self, bareAnalysis):
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            0, extents=_EXTENTS, dxdydz=1)
        assert list(result.xI.values) == [0, 1, 2]
        assert list(result.yI.values) == [0, 1, 2]
        assert list(result.zI.values) == [0, 1, 2]

    def test_a_coarser_cell_size_yields_fewer_grid_points(self, bareAnalysis):
        fine = bareAnalysis.calcConcentrationTimeStepFullMesh(
            0, extents=dict(xmin=0, xmax=10, ymin=0, ymax=10, zmin=0, zmax=10), dxdydz=1)
        coarse = bareAnalysis.calcConcentrationTimeStepFullMesh(
            0, extents=dict(xmin=0, xmax=10, ymin=0, ymax=10, zmin=0, zmax=10), dxdydz=5)
        assert len(fine.xI) > len(coarse.xI)

    def test_particle_mass_is_binned_into_its_cell_and_divided_by_the_volume(
            self, bareAnalysis):
        frame = pandas.DataFrame({
            "datetime": [0.0, 0.0],
            "x": [0.25, 0.75], "y": [0.25, 0.25], "z": [0.25, 0.25],
            "mass": [1.0, 3.0],
        })
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            frame, extents=_EXTENTS, dxdydz=1)
        assert result.sel(xI=0, yI=0, zI=0).item() == pytest.approx(4.0)

    def test_the_field_units_are_recorded_in_the_attributes(self, bareAnalysis):
        frame = pandas.DataFrame({
            "datetime": [0.0], "x": [0.25], "y": [0.25], "z": [0.25], "mass": [1.0]})
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            frame, extents=_EXTENTS, dxdydz=1)
        assert result.attrs["field"] == "1*kg/m**3"

    def test_the_timestamp_comes_from_the_data_not_the_caller(self, bareAnalysis):
        frame = pandas.DataFrame({
            "datetime": [42.0], "x": [0.25], "y": [0.25], "z": [0.25], "mass": [1.0]})
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            frame, extents=_EXTENTS, dxdydz=1)
        assert list(result.datetime.values) == [42.0]

    def test_the_coordinate_column_names_are_configurable(self, bareAnalysis):
        frame = pandas.DataFrame({
            "datetime": [0.0], "east": [0.25], "north": [0.25], "up": [0.25],
            "mass": [5.0]})
        result = bareAnalysis.calcConcentrationTimeStepFullMesh(
            frame, extents=_EXTENTS, dxdydz=1,
            xfield="east", yfield="north", zfield="up")
        assert result.sel(xI=0, yI=0, zI=0).item() == pytest.approx(5.0)


# ===========================================================================
# analysis._writeConcentrationPartition
# ===========================================================================

def _oneTimestepGrid(time=0.0, value=1.0):
    array = xarray.DataArray(
        numpy.full((2, 2, 2), value),
        coords=dict(xI=[0, 1], yI=[0, 1], zI=[0, 1]),
        dims=["xI", "yI", "zI"])
    return array.expand_dims(dict(datetime=[time]), axis=-1)


@pytest.mark.unit
class TestWriteConcentrationPartition:
    def test_the_file_name_is_the_partition_id_zero_padded_to_four(self, tmp_path):
        analysis._writeConcentrationPartition([_oneTimestepGrid()], 7, str(tmp_path))
        assert os.path.exists(str(tmp_path / "Concentrations0007.nc"))

    def test_the_grid_indices_are_renamed_to_x_y_and_z(self, tmp_path):
        analysis._writeConcentrationPartition([_oneTimestepGrid()], 0, str(tmp_path))
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0000.nc"))
        try:
            assert set(stored.dims) == {"x", "y", "z", "datetime"}
        finally:
            stored.close()

    def test_the_variable_is_named_c(self, tmp_path):
        analysis._writeConcentrationPartition([_oneTimestepGrid()], 0, str(tmp_path))
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0000.nc"))
        try:
            assert list(stored.data_vars) == ["C"]
        finally:
            stored.close()

    def test_the_dimension_order_is_y_x_z_datetime(self, tmp_path):
        analysis._writeConcentrationPartition([_oneTimestepGrid()], 0, str(tmp_path))
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0000.nc"))
        try:
            assert stored["C"].dims == ("y", "x", "z", "datetime")
        finally:
            stored.close()

    def test_every_timestep_handed_in_is_concatenated_along_datetime(self, tmp_path):
        analysis._writeConcentrationPartition(
            [_oneTimestepGrid(0.0), _oneTimestepGrid(1.0), _oneTimestepGrid(2.0)],
            3, str(tmp_path))
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0003.nc"))
        try:
            assert list(stored.datetime.values) == [0.0, 1.0, 2.0]
        finally:
            stored.close()


# ===========================================================================
# analysis._processConcentrationPartition / _padRemainingTimesteps
# ===========================================================================

@pytest.mark.unit
class TestProcessConcentrationPartition:
    def test_each_datetime_group_becomes_one_slice_of_the_written_file(
            self, bareAnalysis, tmp_path):
        frame = pandas.DataFrame({
            "datetime": [0.0, 1.0, 1.0],
            "x": [0.25, 0.25, 0.75], "y": [0.25] * 3, "z": [0.25] * 3,
            "mass": [1.0, 2.0, 3.0],
        })
        partition = dask.dataframe.from_pandas(frame, npartitions=1).partitions[0]
        lastTime = bareAnalysis._processConcentrationPartition(
            partition, 0, str(tmp_path), _EXTENTS, 1, "x", "y", "z")
        assert lastTime == 1.0
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0000.nc"))
        try:
            assert list(stored.datetime.values) == [0.0, 1.0]
            assert stored["C"].sel(x=0, y=0, z=0, datetime=1.0).item() == pytest.approx(5.0)
        finally:
            stored.close()

    def test_an_empty_partition_writes_nothing_and_reports_time_zero(
            self, bareAnalysis, tmp_path):
        frame = pandas.DataFrame({
            "datetime": [], "x": [], "y": [], "z": [], "mass": []}).astype(float)
        partition = dask.dataframe.from_pandas(frame, npartitions=1).partitions[0]
        lastTime = bareAnalysis._processConcentrationPartition(
            partition, 0, str(tmp_path), _EXTENTS, 1, "x", "y", "z")
        assert lastTime == 0
        assert glob.glob(os.path.join(str(tmp_path), "*.nc")) == []


@pytest.mark.unit
class TestPadRemainingTimesteps:
    def test_the_gap_up_to_the_duration_is_filled_with_zero_fields(
            self, bareAnalysis, tmp_path):
        bareAnalysis._padRemainingTimesteps(
            2, 6, 4, str(tmp_path), _EXTENTS, 1, "x", "y", "z")
        stored = xarray.open_dataset(str(tmp_path / "Concentrations0004.nc"))
        try:
            assert list(stored.datetime.values) == [3, 4, 5]
            assert float(stored["C"].sum()) == 0.0
        finally:
            stored.close()

    def test_nothing_is_written_when_the_data_already_reaches_the_duration(
            self, bareAnalysis, tmp_path):
        bareAnalysis._padRemainingTimesteps(
            10, 6, 1, str(tmp_path), _EXTENTS, 1, "x", "y", "z")
        assert glob.glob(os.path.join(str(tmp_path), "*.nc")) == []

    def test_the_padding_lands_in_its_own_partition_file(self, bareAnalysis, tmp_path):
        bareAnalysis._padRemainingTimesteps(
            0, 2, 12, str(tmp_path), _EXTENTS, 1, "x", "y", "z")
        assert os.path.exists(str(tmp_path / "Concentrations0012.nc"))


# ===========================================================================
# analysis._removeConcentrationCache / _createConcentrationCacheDoc
# ===========================================================================

@pytest.mark.unit
class TestRemoveConcentrationCache:
    def test_every_file_matching_the_resource_glob_is_removed_with_its_directory(
            self, ana, tmp_path):
        cacheDir = tmp_path / "cachedLagrangianData" / "CASE_fullMeshCache_1"
        cacheDir.mkdir(parents=True)
        for index in (0, 1):
            (cacheDir / f"Concentrations000{index}.nc").write_text("x")
        resource = str(cacheDir / "Concentrations*.nc")
        ana.datalayer.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.NETCDF_XARRAY,
            type=ana.DOCTYPE_CONCENTRATION, desc=dict(caseDescriptorName="CASE"))
        doc = ana.datalayer.toolkit.getCacheDocuments(type=ana.DOCTYPE_CONCENTRATION)[0]
        ana._removeConcentrationCache(doc)
        assert not cacheDir.exists()
        assert len(ana.datalayer.toolkit.getCacheDocuments(
            type=ana.DOCTYPE_CONCENTRATION)) == 0

    def test_a_cache_whose_files_are_already_gone_is_still_unregistered(
            self, ana, tmp_path):
        resource = str(tmp_path / "gone" / "Concentrations*.nc")
        ana.datalayer.toolkit.addCacheDocument(
            resource=resource, dataFormat=datatypes.NETCDF_XARRAY,
            type=ana.DOCTYPE_CONCENTRATION, desc=dict(caseDescriptorName="CASE"))
        doc = ana.datalayer.toolkit.getCacheDocuments(type=ana.DOCTYPE_CONCENTRATION)[0]
        ana._removeConcentrationCache(doc)
        assert len(ana.datalayer.toolkit.getCacheDocuments(
            type=ana.DOCTYPE_CONCENTRATION)) == 0


@pytest.mark.unit
class TestCreateConcentrationCacheDoc:
    def test_the_resource_path_is_built_from_the_case_name_and_the_counter(
            self, ana):
        doc = ana._createConcentrationCacheDoc("MY_CASE", dict(dxdydz=1))
        expected = os.path.join(
            ana.datalayer.toolkit.filesDirectory, "cachedLagrangianData",
            "MY_CASE_fullMeshCache_1", "Concentrations*.nc")
        assert doc.resource == expected

    def test_the_counter_advances_so_two_caches_never_collide(self, ana):
        first = ana._createConcentrationCacheDoc("MY_CASE", dict(dxdydz=1))
        second = ana._createConcentrationCacheDoc("MY_CASE", dict(dxdydz=2))
        assert first.resource != second.resource
        assert second.resource.endswith(
            os.path.join("MY_CASE_fullMeshCache_2", "Concentrations*.nc"))

    def test_the_document_is_registered_with_the_metadata_as_its_description(self, ana):
        ana._createConcentrationCacheDoc("MY_CASE", dict(dxdydz=1, extents=None))
        docs = ana.datalayer.toolkit.getCacheDocuments(type=ana.DOCTYPE_CONCENTRATION)
        assert len(docs) == 1
        assert docs[0].desc["dxdydz"] == 1
        assert docs[0].dataFormat == datatypes.NETCDF_XARRAY


# ===========================================================================
# B326: calcConcentrationFieldFullMesh calls the other class's method
# ===========================================================================

@pytest.mark.unit
class TestCalcConcentrationFieldFullMesh:
    def test_the_helper_it_calls_lives_on_the_other_class(self):
        """Characterisation of B326, at its source."""
        assert not hasattr(analysis, "_resolveCaseDescriptorName")
        assert hasattr(Extension, "_resolveCaseDescriptorName")

    @pytest.mark.xfail(
        strict=True,
        reason="B326: calcConcentrationFieldFullMesh's first statement is "
               "self._resolveCaseDescriptorName(caseDescriptor), but that "
               "method is defined on "
               "absractStochasticLagrangianSolver_toolkitExtension, not on "
               "analysis; every other datalayer call in the same method goes "
               "through self.datalayer(.toolkit). The method raises "
               "AttributeError on its first line, so no concentration field "
               "can ever be computed. See the consolidated findings issue.")
    def test_it_should_at_least_get_past_resolving_the_case_name(self, ana, tmp_path):
        ana.calcConcentrationFieldFullMesh("SOME_CASE", dxdydz=1)

    def test_it_currently_raises_on_its_first_statement(self, ana):
        """Characterisation of B326."""
        with pytest.raises(AttributeError, match="_resolveCaseDescriptorName"):
            ana.calcConcentrationFieldFullMesh("SOME_CASE", dxdydz=1)


# ===========================================================================
# B329: the cache API is called on the extension instead of the toolkit
# ===========================================================================

@pytest.mark.unit
class TestCacheApiIsCalledOnTheWrongObject:
    def test_the_extension_exposes_no_cache_api(self, ext):
        """Characterisation of B329, at its source."""
        assert not hasattr(ext, "getCacheDocuments")
        assert not hasattr(ext, "addCacheDocument")
        assert hasattr(ext.toolkit, "getCacheDocuments")
        assert hasattr(ext.toolkit, "addCacheDocument")

    @pytest.mark.xfail(
        strict=True,
        reason="B329: analysis.calcDocumentConcentrationPointWise and "
               "analysis.getConcentrationField call "
               "self.datalayer.getCacheDocuments(...) / "
               "self.datalayer.addCacheDocument(...), but self.datalayer is an "
               "absractStochasticLagrangianSolver_toolkitExtension, which "
               "defines neither -- the cache API belongs to the Project the "
               "toolkit inherits from. The sibling "
               "calcConcentrationFieldFullMesh in the same class correctly "
               "writes self.datalayer.toolkit.getCacheDocuments. See the "
               "consolidated findings issue.")
    def test_looking_up_a_concentration_field_should_query_the_cache(self, ana):
        class _Doc:
            id = "abc123"

        assert ana.getConcentrationField(_Doc()) is None

    def test_looking_up_a_concentration_field_currently_raises(self, ana):
        """Characterisation of B329."""
        class _Doc:
            id = "abc123"

        with pytest.raises(AttributeError, match="getCacheDocuments"):
            ana.getConcentrationField(_Doc())

    def test_the_pointwise_document_calculator_has_the_same_defect(self, ana):
        """Characterisation of B329, at the second site."""
        class _Doc:
            id = "abc123"

        with pytest.raises(AttributeError, match="getCacheDocuments"):
            ana.calcDocumentConcentrationPointWise(_Doc(), dxdydz=1)
