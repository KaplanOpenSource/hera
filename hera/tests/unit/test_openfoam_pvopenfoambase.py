"""openFoam/postProcess/pvOpenFOAMBase.py: ``paraviewOpenFOAM``.

What the environment allows
--------------------------
Neither ``paraview`` nor ``vtkmodules`` can be installed here (paraview
ships as a binary, and ``vtkmodules`` is not in requirements.txt).  The
module's whole VTK import block is wrapped in ``try/except ImportError``,
so under the unit layer it fails on its **first** line and the module ends
up with *none* of ``dsa``, ``vtkMultiBlockDataSet``, ``pvsimple``,
``servermanager``, ``Proxy`` or ``ProxyProperty`` defined -- every method
that touches them raises ``NameError`` as imported.  (The ``paraview`` /
``paraview.simple`` entries in ``_stubs.py`` are bare namespace modules,
not MagicMocks, and they are never reached because
``import vtkmodules...`` on the line above them raises first.)

The tests therefore inject those six names as module globals with
``monkeypatch.setattr(module, name, value, raising=False)``.  No production
code is modified; the injection is per-test and restored.

Covered *deeply* -- real numpy / pandas / xarray / dask, no mock in the
assertion path, because these methods are pure Python over data structures
the test supplies:

  ``__init__`` (component-name map), ``_parsePointSet``, ``_parseTable``,
  ``_getBlockName``, ``_assignBlockName``, ``_parseMultiBlockDataSet``,
  ``_parseVTKData`` (branch dispatch), ``_resolveTimeList``,
  ``_removeOldOutputs``, ``_ensureOutputDirs``, ``_collectTmpFiles``,
  ``_atomicReplace``, ``_cleanupTmpFiles``, ``_writeTimeStepBlocks``
  (block arithmetic), ``writeList`` (parquet branch, written to and read
  back from ``tmp_path``), ``_mergeParquet`` (real parquet round-trip).

Covered *shallowly* -- these methods exist only to hand arguments to
paraview / VTK / zarr, so the honest test is an argument-forwarding and
call-order test against the stand-ins.  Nothing numerical is asserted,
because a stand-in would happily fabricate any number asked of it:

  ``initializeReader``   -> ``pvsimple.OpenFOAMReader`` args + reader wiring
  ``readTimeSteps``      -> ``pvsimple.FindSource`` lookup, the ``Proxy``
                            assertion, ``None``-skipping, and the nested
                            ``debug_proxy_data`` walk
  ``_readTimeStep``      -> ``UpdatePipeline(t)`` then ``servermanager.Fetch``
  ``writeCase``          -> which helper runs, in which order
  ``_mergeToFinalOutput``-> zarr/parquet dispatch
  ``_mergeZarr`` and ``writeList``'s regularMesh branch -> ``zarr`` is not
                            installed either, so ``xarray`` is replaced by a
                            recording stand-in and only the arguments the
                            method builds itself (paths, ``dim="time"``,
                            ``mode='w'``, the ``NotImplementedError``
                            fallback) are asserted.
  the six ``@deprecated`` aliases -> forwarded positional arguments.

Not covered: nothing in the class is skipped.

Four defects are pinned, each provable by reading the source rather than by
watching a stand-in: B310 (``__init__``'s ``name`` argument is discarded),
B311 (``readTimeSteps``'s ``timelist=None`` default is iterated unguarded),
B312 (``_ensureOutputDirs`` calls ``os.makedirs('')``) and B313
(``_parsePointSet``'s ``squeeze()`` cannot tell one point from no points).
Each has a strict xfail for the intended behaviour plus passing
characterisation tests of what actually happens today.
"""
import importlib
import os
import types
from unittest.mock import MagicMock

import numpy
import pandas
import pytest
import xarray

PV_MODULE = "hera.simulations.openFoam.postProcess.pvOpenFOAMBase"


# ---------------------------------------------------------------------------
# Stand-ins for the VTK / paraview objects the module expects
# ---------------------------------------------------------------------------

class FakeVTKNoneArray:
    """Stands in for ``dsa.VTKNoneArray`` -- the "field is absent" marker."""


class FakeVTKArray(numpy.ndarray):
    """Stands in for ``dsa.VTKArray``, which really is an ndarray subclass."""


def vtkarray(values):
    return numpy.asarray(values).view(FakeVTKArray)


class FakeCompositeArray:
    """Neither None nor VTKArray, so the code takes its ``GetArrays`` path."""

    def __init__(self, arrays):
        self._arrays = list(arrays)

    def GetArrays(self):
        return self._arrays


class FakePointSet:
    """Stands in for ``dsa.PointSet``."""

    def __init__(self, points, pointData=None):
        self.Points = points
        self.PointData = dict(pointData or {})


class FakePolyData(FakePointSet):
    pass


class FakeUnstructuredGrid(FakePointSet):
    pass


class FakeTable:
    def __init__(self, rowData):
        self.RowData = dict(rowData)


class FakeCompositeDataSet:
    """Stands in for ``dsa.CompositeDataSet``; unwrapped via ``.VTKObject``."""

    def __init__(self, vtkObject):
        self.VTKObject = vtkObject


class FakeMetaData:
    def __init__(self, values):
        self._values = dict(values)

    def Has(self, key):
        return key in self._values

    def Get(self, key):
        return self._values[key]


class FakeMultiBlockDataSet:
    """Stands in for ``vtkMultiBlockDataSet``."""

    def __init__(self, blocks=(), metadata=None):
        self._blocks = list(blocks)
        self._metadata = list(metadata) if metadata is not None else [None] * len(self._blocks)

    def GetNumberOfBlocks(self):
        return len(self._blocks)

    def GetBlock(self, i):
        return self._blocks[i]

    def HasMetaData(self, i):
        return self._metadata[i] is not None

    def GetMetaData(self, i):
        return self._metadata[i]

    def NAME(self):
        return "NAME"


def named_block(name):
    return FakeMetaData({"NAME": name})


class FakeProxyProperty:
    def __init__(self, data):
        self._data = data

    def GetData(self):
        return self._data


class FakeProxy:
    """Stands in for ``paraview.servermanager.Proxy``."""

    def __init__(self, properties=None):
        self._properties = dict(properties or {})
        self.listPropertiesCalls = 0
        self.getPropertyCalls = []

    def ListProperties(self):
        self.listPropertiesCalls += 1
        return list(self._properties)

    def GetProperty(self, name):
        self.getPropertyCalls.append(name)
        return self._properties[name]


class FakeMeshRegions(list):
    """``reader.MeshRegions`` is both iterable and has ``SelectAll()``."""

    def __init__(self, items):
        super().__init__(items)
        self.selectAllCalls = 0

    def SelectAll(self):
        self.selectAllCalls += 1


class FakeReader:
    def __init__(self, regions=("internalMesh",), timesteps=()):
        self.MeshRegions = FakeMeshRegions(regions)
        self.TimestepValues = list(timesteps)
        self.updatePipelineCalls = []

    def UpdatePipeline(self, time=None):
        self.updatePipelineCalls.append(time)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def pvmod():
    return importlib.import_module(PV_MODULE)


@pytest.fixture()
def pv(pvmod, monkeypatch):
    """Inject the six missing VTK/paraview module globals.

    Returns a namespace exposing the module plus every stand-in, so a test
    can inspect what the production code did with them.
    """
    dsa = types.SimpleNamespace(
        VTKNoneArray=FakeVTKNoneArray,
        VTKArray=FakeVTKArray,
        PointSet=FakePointSet,
        PolyData=FakePolyData,
        UnstructuredGrid=FakeUnstructuredGrid,
        Table=FakeTable,
        CompositeDataSet=FakeCompositeDataSet,
        WrapDataObject=MagicMock(side_effect=lambda data: data),
    )
    reader = FakeReader()
    pvsimple = MagicMock()
    pvsimple.OpenFOAMReader.return_value = reader
    servermanager = MagicMock()

    injected = dict(
        dsa=dsa,
        vtkMultiBlockDataSet=FakeMultiBlockDataSet,
        pvsimple=pvsimple,
        servermanager=servermanager,
        Proxy=FakeProxy,
        ProxyProperty=FakeProxyProperty,
    )
    for name, value in injected.items():
        monkeypatch.setattr(pvmod, name, value, raising=False)

    return types.SimpleNamespace(
        module=pvmod,
        cls=pvmod.paraviewOpenFOAM,
        dsa=dsa,
        pvsimple=pvsimple,
        servermanager=servermanager,
        reader=reader,
    )


@pytest.fixture()
def base(pvmod):
    """A ``paraviewOpenFOAM`` needing no paraview at all (servername=None)."""
    return pvmod.paraviewOpenFOAM(casePath="/cases/demo")


def fake_xarray(**overrides):
    """A recording stand-in for the module's ``xarray`` global.

    ``zarr`` is not installed, so the zarr write paths can only be checked
    at the argument-forwarding level.
    """
    namespace = types.SimpleNamespace(
        concat=MagicMock(),
        open_mfdataset=MagicMock(),
        Dataset=xarray.Dataset,
        DataArray=xarray.DataArray,
    )
    for key, value in overrides.items():
        setattr(namespace, key, value)
    return namespace


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConstruction:
    def test_it_stores_the_case_path_and_case_type(self, pvmod):
        obj = pvmod.paraviewOpenFOAM(casePath="/cases/demo", caseType="Reconstructed Case")
        assert obj.casePath == "/cases/demo"
        assert obj.caseType == "Reconstructed Case"

    def test_the_case_type_defaults_to_the_decomposed_constant(self, pvmod):
        from hera.simulations.openFoam import CASETYPE_DECOMPOSED

        assert pvmod.paraviewOpenFOAM(casePath="/c").caseType == CASETYPE_DECOMPOSED

    def test_it_builds_a_component_suffix_for_every_scalar_vector_and_tensor_index(self, base):
        names = base._componentsNames
        assert len(names) == 13
        assert names[()] == ""
        assert [names[(i,)] for i in range(3)] == ["_x", "_y", "_z"]
        assert names[(0, 0)] == "_xx"
        assert names[(1, 2)] == "_yz"
        assert names[(2, 0)] == "_zx"

    def test_the_component_map_is_per_instance_not_shared_on_the_class(self, pvmod):
        first = pvmod.paraviewOpenFOAM(casePath="/a")
        second = pvmod.paraviewOpenFOAM(casePath="/b")
        assert first._componentsNames is not second._componentsNames
        assert pvmod.paraviewOpenFOAM._componentsNames is None

    def test_it_does_not_connect_to_a_server_when_no_servername_is_given(self, pv):
        pv.cls(casePath="/c", servername=None)
        assert pv.pvsimple.Connect.call_count == 0

    def test_it_connects_to_the_named_paraview_server(self, pv):
        pv.cls(casePath="/c", servername="cs://host:11111")
        pv.pvsimple.Connect.assert_called_once_with("cs://host:11111")


# ---------------------------------------------------------------------------
# B310: __init__'s `name` argument is accepted and thrown away
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestConstructorNameArgument:
    @pytest.mark.xfail(
        strict=True,
        reason="B310: paraviewOpenFOAM.__init__ declares name='mainreader' "
               "but its body never reads it -- the parameter is bound and "
               "immediately discarded (source: lines 39-91 assign only "
               "_componentsNames, casePath and caseType). Callers who name "
               "the pipeline at construction time are silently ignored, and "
               "initializeReader falls back to its own unrelated default "
               "readerName='reader'. See the consolidated findings issue.",
    )
    def test_the_constructor_name_becomes_the_readers_default_pipeline_name(self, pv):
        obj = pv.cls(casePath="/cases/demo", name="myPipeline")
        obj.initializeReader()
        assert obj.readerName == "myPipeline"

    def test_the_constructor_name_is_not_recorded_anywhere_on_the_instance(self, pvmod):
        """Characterisation of B310."""
        obj = pvmod.paraviewOpenFOAM(casePath="/cases/demo", name="myPipeline")
        assert "myPipeline" not in vars(obj).values()
        assert not hasattr(obj, "name")
        assert sorted(vars(obj)) == ["_componentsNames", "casePath", "caseType"]

    def test_the_reader_name_comes_only_from_initialize_reader(self, pv):
        """Characterisation of B310."""
        obj = pv.cls(casePath="/cases/demo", name="myPipeline")
        obj.initializeReader()
        assert obj.readerName == "reader"
        assert pv.pvsimple.OpenFOAMReader.call_args.kwargs["guiName"] == "reader"


# ---------------------------------------------------------------------------
# initializeReader -- shallow: argument forwarding to pvsimple
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestInitializeReader:
    def test_it_asks_paraview_for_a_reader_on_the_cases_tmp_foam_file(self, pv):
        obj = pv.cls(casePath="/cases/demo", caseType="Reconstructed Case")
        obj.initializeReader(readerName="myReader")
        pv.pvsimple.OpenFOAMReader.assert_called_once_with(
            FileName="/cases/demo/tmp.foam",
            CaseType="Reconstructed Case",
            guiName="myReader",
        )

    def test_it_selects_every_mesh_region_and_writes_the_list_back(self, pv):
        pv.reader.MeshRegions = FakeMeshRegions(["internalMesh", "patch/inlet"])
        obj = pv.cls(casePath="/c")
        reader = obj.initializeReader()
        assert reader.MeshRegions == ["internalMesh", "patch/inlet"]
        assert not isinstance(reader.MeshRegions, FakeMeshRegions)

    def test_it_updates_the_pipeline_once_with_no_time_argument(self, pv):
        pv.cls(casePath="/c").initializeReader()
        assert pv.reader.updatePipelineCalls == [None]

    def test_it_returns_the_reader_and_remembers_it_with_its_name(self, pv):
        obj = pv.cls(casePath="/c")
        reader = obj.initializeReader(readerName="pipe1")
        assert reader is pv.reader
        assert obj.reader is pv.reader
        assert obj.readerName == "pipe1"


# ---------------------------------------------------------------------------
# readTimeSteps -- shallow: source lookup, Proxy guard, None-skipping
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestReadTimeSteps:
    def test_it_yields_one_dictionary_per_timestep_keyed_by_filter_name(self, pv, monkeypatch):
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(
            pv.cls, "_readTimeStep",
            lambda self, datasource, timeslice, fieldnames, regularMesh: f"data@{timeslice}",
        )
        obj = pv.cls(casePath="/c")
        result = list(obj.readTimeSteps({"slice1": "/out/a.parquet"}, timelist=[0.0, 1.0]))
        assert result == [{"slice1": "data@0.0"}, {"slice1": "data@1.0"}]

    def test_it_looks_each_filter_up_by_name_in_the_paraview_pipeline(self, pv, monkeypatch):
        monkeypatch.setattr(pv.pvsimple, "FindSource",
                            MagicMock(side_effect=lambda name: FakeProxy()))
        monkeypatch.setattr(pv.cls, "_readTimeStep", lambda *a, **k: "x")
        obj = pv.cls(casePath="/c")
        list(obj.readTimeSteps({"a": "/out/a", "b": "/out/b"}, timelist=[1.0]))
        assert [c.args[0] for c in pv.pvsimple.FindSource.call_args_list] == ["a", "b"]

    def test_it_omits_filters_whose_timestep_produced_no_data(self, pv, monkeypatch):
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(
            pv.cls, "_readTimeStep",
            lambda self, datasource, timeslice, fieldnames, regularMesh: None if timeslice == 1.0 else 7,
        )
        obj = pv.cls(casePath="/c")
        assert list(obj.readTimeSteps({"a": "/o/a"}, timelist=[0.0, 1.0])) == [{"a": 7}, {}]

    def test_it_refuses_a_pipeline_object_that_is_not_a_paraview_proxy(self, pv, monkeypatch):
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: object())
        obj = pv.cls(casePath="/c")
        with pytest.raises(AssertionError):
            list(obj.readTimeSteps({"a": "/o/a"}, timelist=[0.0]))

    def test_it_walks_nested_proxy_properties_before_computing_the_filter(self, pv, monkeypatch):
        """debug_proxy_data must visit plain properties, proxy properties
        holding data, and proxy properties holding another proxy."""
        nested = FakeProxy({"Radius": 3.0})
        source = FakeProxy({
            "PlainProperty": "just a value",
            "DataProperty": FakeProxyProperty("not a proxy"),
            "NestedProxy": FakeProxyProperty(nested),
        })
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: source)
        monkeypatch.setattr(pv.cls, "_readTimeStep", lambda *a, **k: "x")
        obj = pv.cls(casePath="/c")
        list(obj.readTimeSteps({"a": "/o/a"}, timelist=[0.0]))
        assert source.getPropertyCalls == ["PlainProperty", "DataProperty", "NestedProxy"]
        assert nested.getPropertyCalls == ["Radius"]

    def test_the_filters_output_paths_are_never_read_only_their_names(self, pv, monkeypatch):
        """readTimeSteps takes filtername -> path but only uses the keys."""
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(pv.cls, "_readTimeStep", lambda *a, **k: "x")
        obj = pv.cls(casePath="/c")
        assert list(obj.readTimeSteps({"a": None}, timelist=[0.0])) == [{"a": "x"}]


# ---------------------------------------------------------------------------
# B311: readTimeSteps declares timelist=None but iterates it unguarded
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestReadTimeStepsDefaultTimeList:
    @pytest.mark.xfail(
        strict=True,
        reason="B311: readTimeSteps(datasourcenamedict, timelist=None, ...) "
               "iterates `for timeslice in timelist:` with no None guard, so "
               "the declared default is unusable -- it raises TypeError: "
               "'NoneType' object is not iterable. The same class already "
               "shows the intended handling in _resolveTimeList "
               "('self.reader.TimestepValues if timeList is None else "
               "timeList') and writeCase documents 'None = all available "
               "from the reader'. Four deprecated aliases (to_pandas, "
               "to_xarray, to_dataFrame, to_dataArray) forward the same "
               "None default straight through. See the consolidated "
               "findings issue.",
    )
    def test_a_missing_time_list_means_every_timestep_the_reader_knows(self, pv, monkeypatch):
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(pv.cls, "_readTimeStep", lambda *a, **k: "x")
        obj = pv.cls(casePath="/c")
        obj.reader = FakeReader(timesteps=[0.0, 1.0, 2.0])
        assert list(obj.readTimeSteps({"a": "/o/a"}, timelist=None)) == [{"a": "x"}] * 3

    def test_a_missing_time_list_raises_a_type_error(self, pv):
        """Characterisation of B311."""
        obj = pv.cls(casePath="/c")
        obj.reader = FakeReader(timesteps=[0.0, 1.0])
        with pytest.raises(TypeError, match="not iterable"):
            list(obj.readTimeSteps({"a": "/o/a"}))

    @pytest.mark.parametrize("alias", ["to_pandas", "to_xarray", "to_dataFrame", "to_dataArray"])
    def test_the_deprecated_read_aliases_inherit_the_same_type_error(self, pv, alias):
        """Characterisation of B311."""
        obj = pv.cls(casePath="/c")
        with pytest.raises(TypeError, match="not iterable"):
            list(getattr(obj, alias)({"a": "/o/a"}))


# ---------------------------------------------------------------------------
# _parsePointSet -- deep: real numpy in, real pandas/xarray out
# ---------------------------------------------------------------------------

POINTS_3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]


@pytest.mark.unit
class TestParsePointSet:
    def test_a_filter_with_no_points_at_all_produces_nothing(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(points=FakeVTKNoneArray())
        assert obj._parsePointSet(pointSet, 1.0, None, False) is None

    def test_it_turns_the_point_coordinates_into_x_y_z_and_time_columns(self, pv):
        obj = pv.cls(casePath="/c")
        result = obj._parsePointSet(FakePointSet(vtkarray(POINTS_3)), 2.5, None, False)
        assert list(result.columns) == ["x", "y", "z", "time"]
        assert result.x.tolist() == [1.0, 4.0, 7.0]
        assert result.z.tolist() == [3.0, 6.0, 9.0]
        assert result.time.tolist() == [2.5, 2.5, 2.5]

    def test_it_rounds_the_coordinates_and_the_time_to_seven_decimals(self, pv):
        obj = pv.cls(casePath="/c")
        points = vtkarray([[1.123456789, 2.0, 3.0], [0.0, 0.0, 0.0]])
        result = obj._parsePointSet(FakePointSet(points), 0.123456789, None, False)
        assert result.x.tolist() == [1.1234568, 0.0]
        assert result.time.tolist() == [0.1234568, 0.1234568]

    def test_it_stitches_together_the_point_arrays_of_a_composite_data_set(self, pv):
        obj = pv.cls(casePath="/c")
        points = FakeCompositeArray([
            numpy.array([[1.0, 1.0, 1.0]]),
            numpy.array([[2.0, 2.0, 2.0]]),
        ])
        result = obj._parsePointSet(FakePointSet(points), 0.0, None, False)
        assert result.x.tolist() == [1.0, 2.0]

    def test_a_scalar_field_keeps_its_bare_name(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": vtkarray([300.0, 301.0, 302.0])})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert result["T"].tolist() == [300.0, 301.0, 302.0]

    def test_a_vector_field_is_split_into_underscore_x_y_z_columns(self, pv):
        obj = pv.cls(casePath="/c")
        velocity = vtkarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        pointSet = FakePointSet(vtkarray(POINTS_3), {"U": velocity})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert [c for c in result.columns if c.startswith("U")] == ["U_x", "U_y", "U_z"]
        assert result["U_y"].tolist() == [2.0, 5.0, 8.0]

    def test_a_tensor_field_is_split_into_all_nine_named_components(self, pv):
        obj = pv.cls(casePath="/c")
        tensor = vtkarray(numpy.arange(3 * 3 * 3, dtype=float).reshape(3, 3, 3))
        pointSet = FakePointSet(vtkarray(POINTS_3), {"R": tensor})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert [c for c in result.columns if c.startswith("R")] == [
            "R_xx", "R_xy", "R_xz", "R_yx", "R_yy", "R_yz", "R_zx", "R_zy", "R_zz",
        ]
        assert result["R_zy"].tolist() == tensor[:, 2, 1].tolist()

    def test_it_stitches_together_the_field_arrays_of_a_composite_data_set(self, pv):
        obj = pv.cls(casePath="/c")
        field = FakeCompositeArray([numpy.array([1.0, 2.0]), numpy.array([3.0])])
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": field})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert result["T"].tolist() == [1.0, 2.0, 3.0]

    def test_absent_sub_arrays_are_dropped_before_stitching_a_composite_field(self, pv):
        obj = pv.cls(casePath="/c")
        field = FakeCompositeArray([
            numpy.array([1.0, 2.0]), FakeVTKNoneArray(), numpy.array([3.0]),
        ])
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": field})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert result["T"].tolist() == [1.0, 2.0, 3.0]

    def test_a_field_that_is_absent_for_this_filter_is_skipped(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": FakeVTKNoneArray()})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert "T" not in result.columns

    def test_an_explicit_field_list_wins_over_everything_the_filter_carries(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray(POINTS_3), {
            "T": vtkarray([1.0, 2.0, 3.0]),
            "p": vtkarray([4.0, 5.0, 6.0]),
        })
        result = obj._parsePointSet(pointSet, 0.0, ["p"], False)
        assert "p" in result.columns
        assert "T" not in result.columns

    def test_a_field_whose_length_disagrees_with_the_point_count_is_omitted(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": vtkarray([1.0, 2.0])})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert "T" not in result.columns
        assert list(result.columns) == ["x", "y", "z", "time"]

    def test_an_integrated_filter_with_a_single_point_row_carries_only_time(self, pv):
        """Filters such as IntegrateVariables squeeze away the point axis."""
        obj = pv.cls(casePath="/c")
        result = obj._parsePointSet(FakePointSet(vtkarray([1.0, 2.0, 3.0])), 4.0, None, False)
        assert list(result.columns) == ["time"]
        assert result.time.tolist() == [4.0]

    def test_an_integrated_zero_dimensional_field_becomes_one_scalar_column(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([1.0, 2.0, 3.0]), {"Volume": vtkarray(12.5)})
        result = obj._parsePointSet(pointSet, 4.0, None, False)
        assert result["Volume"].tolist() == [12.5]

    def test_a_regular_mesh_is_indexed_by_time_and_the_three_coordinates(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray(POINTS_3), {"T": vtkarray([1.0, 2.0, 3.0])})
        result = obj._parsePointSet(pointSet, 0.0, None, True)
        assert isinstance(result, xarray.Dataset)
        assert list(result.dims) == ["time", "x", "y", "z"]
        assert "T" in result.data_vars

    def test_a_regular_mesh_without_coordinates_is_indexed_by_time_alone(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([1.0, 2.0, 3.0]), {"Volume": vtkarray(12.5)})
        result = obj._parsePointSet(pointSet, 4.0, None, True)
        assert isinstance(result, xarray.Dataset)
        assert list(result.dims) == ["time"]
        assert result["Volume"].values.tolist() == [12.5]


# ---------------------------------------------------------------------------
# B313: squeeze() cannot tell "no points" from "exactly one point"
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestParsePointSetSinglePoint:
    @pytest.mark.xfail(
        strict=True,
        reason="B313: _parsePointSet does "
               "points = numpy.array(pointSet.Points).squeeze() and then "
               "keeps the x/y/z columns only `if len(points.shape)==2`, with "
               "the comment 'it means that there are no points'. That "
               "inference is wrong for a filter with exactly ONE point (e.g. "
               "ProbeLocation / a one-point PointSource): numpy.squeeze "
               "drops every length-1 axis, so an (1, 3) point array becomes "
               "(3,) and the branch silently throws the probe's coordinates "
               "away, returning a row that carries only `time`. The point "
               "count is available un-squeezed and is what the guard should "
               "test. See the consolidated findings issue.",
    )
    def test_a_filter_with_exactly_one_point_keeps_that_points_coordinates(self, pv):
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([[1.0, 2.0, 3.0]]))
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert list(result.columns) == ["x", "y", "z", "time"]
        assert result.x.tolist() == [1.0]

    def test_a_single_point_filter_comes_back_with_only_a_time_column(self, pv):
        """Characterisation of B313."""
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([[1.0, 2.0, 3.0]]))
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert list(result.columns) == ["time"]

    def test_a_single_point_scalar_field_still_survives_as_a_scalar(self, pv):
        """Characterisation of B313: the field value is kept, its place is not."""
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([[1.0, 2.0, 3.0]]), {"T": vtkarray([300.0])})
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert result["T"].tolist() == [300.0]
        assert "x" not in result.columns

    def test_two_points_are_enough_for_the_coordinates_to_be_kept(self, pv):
        """Characterisation of B313: the boundary is exactly one point."""
        obj = pv.cls(casePath="/c")
        pointSet = FakePointSet(vtkarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        result = obj._parsePointSet(pointSet, 0.0, None, False)
        assert list(result.columns) == ["x", "y", "z", "time"]


# ---------------------------------------------------------------------------
# _parseTable -- deep
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestParseTable:
    def test_every_row_data_key_becomes_a_column_plus_a_time_column(self, pv):
        obj = pv.cls(casePath="/c")
        table = FakeTable({"a": [1, 2], "b": [3.0, 4.0]})
        result = obj._parseTable(table, 7.0, False)
        assert list(result.columns) == ["a", "b", "time"]
        assert result.a.tolist() == [1, 2]
        assert result.time.tolist() == [7.0, 7.0]

    def test_an_empty_table_still_gets_a_time_column(self, pv):
        obj = pv.cls(casePath="/c")
        result = obj._parseTable(FakeTable({}), 7.0, False)
        assert list(result.columns) == ["time"]
        assert len(result) == 0

    def test_a_regular_mesh_table_comes_back_as_an_xarray_dataset(self, pv):
        obj = pv.cls(casePath="/c")
        result = obj._parseTable(FakeTable({"a": [1, 2]}), 7.0, True)
        assert isinstance(result, xarray.Dataset)
        assert result["a"].values.tolist() == [1, 2]


# ---------------------------------------------------------------------------
# _getBlockName -- deep
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestGetBlockName:
    def test_a_block_with_no_metadata_has_no_name(self, pv):
        obj = pv.cls(casePath="/c")
        multiBlock = FakeMultiBlockDataSet(blocks=["b"], metadata=[None])
        assert numpy.isnan(obj._getBlockName(multiBlock, 0))

    def test_metadata_that_does_not_carry_the_name_key_yields_no_name(self, pv):
        obj = pv.cls(casePath="/c")
        multiBlock = FakeMultiBlockDataSet(blocks=["b"], metadata=[FakeMetaData({"Other": 1})])
        assert numpy.isnan(obj._getBlockName(multiBlock, 0))

    def test_it_reads_the_name_the_block_metadata_advertises(self, pv):
        obj = pv.cls(casePath="/c")
        multiBlock = FakeMultiBlockDataSet(
            blocks=["b0", "b1"], metadata=[named_block("inlet"), named_block("outlet")],
        )
        assert obj._getBlockName(multiBlock, 1) == "outlet"


# ---------------------------------------------------------------------------
# _assignBlockName -- deep
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAssignBlockName:
    def test_a_pandas_block_gets_the_name_as_a_broadcast_column(self, pv):
        obj = pv.cls(casePath="/c")
        frame = pandas.DataFrame({"x": [1, 2, 3]})
        result = obj._assignBlockName(False, "inlet", frame)
        assert result["blockName"].tolist() == ["inlet"] * 3

    def test_a_pandas_block_is_annotated_in_place_and_returned(self, pv):
        obj = pv.cls(casePath="/c")
        frame = pandas.DataFrame({"x": [1]})
        assert obj._assignBlockName(False, "inlet", frame) is frame

    def test_an_xarray_block_gets_the_name_as_a_variable_over_all_its_dims(self, pv):
        obj = pv.cls(casePath="/c")
        dataset = xarray.Dataset({"T": (("time", "x"), numpy.zeros((2, 3)))})
        result = obj._assignBlockName(True, "inlet", dataset)
        assert result["blockName"].dims == ("time", "x")
        assert result["blockName"].shape == (2, 3)
        assert set(numpy.unique(result["blockName"].values)) == {"inlet"}

    def test_an_xarray_data_array_block_is_promoted_to_a_dataset_first(self, pv):
        obj = pv.cls(casePath="/c")
        array = xarray.DataArray(numpy.zeros((2,)), dims=("time",), name="T")
        result = obj._assignBlockName(True, "inlet", array)
        assert isinstance(result, xarray.Dataset)
        assert set(result.data_vars) == {"T", "blockName"}

    def test_a_regular_mesh_block_that_is_not_xarray_is_rejected(self, pv):
        obj = pv.cls(casePath="/c")
        with pytest.raises(AssertionError, match="instead of xarray"):
            obj._assignBlockName(True, "inlet", pandas.DataFrame({"x": [1]}))

    def test_a_non_regular_mesh_block_that_is_not_a_data_frame_is_rejected(self, pv):
        obj = pv.cls(casePath="/c")
        with pytest.raises(AssertionError, match="instead of pandas.DataFrame"):
            obj._assignBlockName(False, "inlet", xarray.Dataset({"T": (("time",), [0.0])}))


# ---------------------------------------------------------------------------
# _parseMultiBlockDataSet -- deep
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestParseMultiBlockDataSet:
    def test_a_multi_block_with_no_blocks_produces_an_empty_list(self, pv):
        obj = pv.cls(casePath="/c")
        assert obj._parseMultiBlockDataSet(FakeMultiBlockDataSet(), 0.0, None, False) == []

    def test_a_single_block_is_returned_unwrapped_and_unnamed(self, pv):
        obj = pv.cls(casePath="/c")
        block = FakePointSet(vtkarray(POINTS_3))
        result = obj._parseMultiBlockDataSet(
            FakeMultiBlockDataSet([block], [named_block("only")]), 0.0, None, False,
        )
        assert isinstance(result, pandas.DataFrame)
        assert "blockName" not in result.columns

    def test_several_blocks_come_back_as_a_list_each_tagged_with_its_name(self, pv):
        obj = pv.cls(casePath="/c")
        blocks = [FakePointSet(vtkarray(POINTS_3)), FakePointSet(vtkarray(POINTS_3))]
        result = obj._parseMultiBlockDataSet(
            FakeMultiBlockDataSet(blocks, [named_block("inlet"), named_block("outlet")]),
            0.0, None, False,
        )
        assert len(result) == 2
        assert result[0]["blockName"].unique().tolist() == ["inlet"]
        assert result[1]["blockName"].unique().tolist() == ["outlet"]

    def test_a_nested_multi_block_is_flattened_rather_than_tagged_again(self, pv):
        obj = pv.cls(casePath="/c")
        inner = FakeMultiBlockDataSet(
            [FakePointSet(vtkarray(POINTS_3)), FakePointSet(vtkarray(POINTS_3))],
            [named_block("a"), named_block("b")],
        )
        outer = FakeMultiBlockDataSet(
            [inner, FakePointSet(vtkarray(POINTS_3))], [named_block("group"), named_block("c")],
        )
        result = obj._parseMultiBlockDataSet(outer, 0.0, None, False)
        assert len(result) == 3
        assert [frame["blockName"].unique()[0] for frame in result] == ["a", "b", "c"]

    def test_a_block_the_data_set_cannot_hand_over_is_a_runtime_error(self, pv):
        obj = pv.cls(casePath="/c")
        multiBlock = FakeMultiBlockDataSet(
            [FakePointSet(vtkarray(POINTS_3)), None], [named_block("a"), named_block("b")],
        )
        with pytest.raises(RuntimeError, match="Failed to get block number 1"):
            obj._parseMultiBlockDataSet(multiBlock, 0.0, None, False)


# ---------------------------------------------------------------------------
# _parseVTKData -- deep on the dispatch, shallow on the wrapping
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestParseVTKData:
    def test_raw_vtk_data_is_handed_to_the_numpy_wrapper_first(self, pv):
        obj = pv.cls(casePath="/c")
        raw = FakePointSet(vtkarray(POINTS_3))
        obj._parseVTKData(raw, 0.0, None, False)
        pv.dsa.WrapDataObject.assert_called_once_with(raw)

    def test_an_already_composite_data_set_is_not_wrapped_again(self, pv):
        obj = pv.cls(casePath="/c")
        obj._parseVTKData(FakeMultiBlockDataSet(), 0.0, None, False)
        assert pv.dsa.WrapDataObject.call_count == 0

    def test_a_wrapped_composite_data_set_is_unwrapped_back_to_its_vtk_object(self, pv):
        obj = pv.cls(casePath="/c")
        blocks = [FakePointSet(vtkarray(POINTS_3)), FakePointSet(vtkarray(POINTS_3))]
        multiBlock = FakeMultiBlockDataSet(blocks, [named_block("a"), named_block("b")])
        result = obj._parseVTKData(FakeCompositeDataSet(multiBlock), 0.0, None, False)
        assert [frame["blockName"].unique()[0] for frame in result] == ["a", "b"]

    @pytest.mark.parametrize("pointSetClass", [FakePointSet, FakePolyData, FakeUnstructuredGrid])
    def test_every_point_bearing_data_set_flavour_goes_to_the_point_set_parser(self, pv, pointSetClass):
        obj = pv.cls(casePath="/c")
        result = obj._parseVTKData(pointSetClass(vtkarray(POINTS_3)), 1.0, None, False)
        assert isinstance(result, pandas.DataFrame)
        assert result.x.tolist() == [1.0, 4.0, 7.0]

    def test_table_data_goes_to_the_table_parser(self, pv):
        obj = pv.cls(casePath="/c")
        result = obj._parseVTKData(FakeTable({"a": [1]}), 1.0, None, False)
        assert list(result.columns) == ["a", "time"]

    def test_an_unsupported_vtk_type_is_reported_as_not_implemented(self, pv):
        obj = pv.cls(casePath="/c")
        with pytest.raises(NotImplementedError, match="doesn't implement parsing of"):
            obj._parseVTKData(object(), 1.0, None, False)


# ---------------------------------------------------------------------------
# _readTimeStep -- shallow: pipeline update then Fetch
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestReadTimeStep:
    def test_it_advances_the_filter_to_the_requested_time_before_fetching(self, pv):
        obj = pv.cls(casePath="/c")
        datasource = FakeReader()
        pv.servermanager.Fetch.return_value = FakeTable({"a": [1]})
        obj._readTimeStep(datasource, 3.5)
        assert datasource.updatePipelineCalls == [3.5]
        pv.servermanager.Fetch.assert_called_once_with(datasource)

    def test_it_parses_whatever_the_server_hands_back(self, pv):
        obj = pv.cls(casePath="/c")
        pv.servermanager.Fetch.return_value = FakePointSet(vtkarray(POINTS_3))
        result = obj._readTimeStep(FakeReader(), 3.5)
        assert result.time.tolist() == [3.5, 3.5, 3.5]

    def test_the_field_list_and_mesh_flag_reach_the_parser(self, pv, monkeypatch):
        obj = pv.cls(casePath="/c")
        pv.servermanager.Fetch.return_value = "raw"
        seen = {}
        monkeypatch.setattr(
            pv.cls, "_parseVTKData",
            lambda self, data, timeslice, fieldnames, regularMesh: seen.update(
                data=data, timeslice=timeslice, fieldnames=fieldnames, regularMesh=regularMesh,
            ),
        )
        obj._readTimeStep(FakeReader(), 3.5, fieldnames=["T"], regularMesh=True)
        assert seen == dict(data="raw", timeslice=3.5, fieldnames=["T"], regularMesh=True)


# ---------------------------------------------------------------------------
# _resolveTimeList -- deep
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestResolveTimeList:
    def test_no_time_list_means_every_timestep_the_reader_advertises(self, base):
        base.reader = FakeReader(timesteps=[0.0, 1.0, 2.0])
        assert base._resolveTimeList(None, False) == [0.0, 1.0, 2.0]

    def test_a_caller_supplied_time_list_is_used_verbatim(self, base):
        base.reader = FakeReader(timesteps=[0.0, 1.0, 2.0])
        assert base._resolveTimeList([5.0, 6.0], False) == [5.0, 6.0]

    def test_asking_for_the_latest_timestamp_keeps_only_the_last_entry(self, base):
        base.reader = FakeReader(timesteps=[0.0, 1.0, 2.0])
        assert base._resolveTimeList(None, True) == [2.0]

    def test_asking_for_the_latest_timestamp_of_nothing_stays_empty(self, base):
        base.reader = FakeReader(timesteps=[])
        assert base._resolveTimeList(None, True) == []

    def test_the_reader_is_never_consulted_when_a_time_list_is_given(self, base):
        assert base._resolveTimeList([1.0], False) == [1.0]
        assert not hasattr(base, "reader")


# ---------------------------------------------------------------------------
# Filesystem helpers -- deep, against tmp_path
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestRemoveOldOutputs:
    def test_an_existing_output_file_is_deleted(self, base, tmp_path):
        target = tmp_path / "a.parquet"
        target.write_text("x")
        base._removeOldOutputs({"a": str(target)})
        assert not target.exists()

    def test_an_existing_output_directory_is_deleted_whole(self, base, tmp_path):
        target = tmp_path / "a.parquet"
        (target / "part").mkdir(parents=True)
        base._removeOldOutputs({"a": str(target)})
        assert not target.exists()

    def test_an_output_that_was_never_written_is_left_alone(self, base, tmp_path):
        base._removeOldOutputs({"a": str(tmp_path / "missing.parquet")})
        assert sorted(p.name for p in tmp_path.iterdir()) == []

    def test_every_filter_in_the_dictionary_is_cleaned(self, base, tmp_path):
        for name in ("a", "b"):
            (tmp_path / f"{name}.parquet").write_text("x")
        base._removeOldOutputs({n: str(tmp_path / f"{n}.parquet") for n in ("a", "b")})
        assert list(tmp_path.iterdir()) == []


@pytest.mark.unit
class TestEnsureOutputDirs:
    def test_the_directory_holding_the_output_file_is_created(self, base, tmp_path):
        base._ensureOutputDirs({"a": str(tmp_path / "nested" / "deep" / "a.parquet")})
        assert (tmp_path / "nested" / "deep").is_dir()

    def test_an_existing_directory_is_accepted_without_complaint(self, base, tmp_path):
        base._ensureOutputDirs({"a": str(tmp_path / "a.parquet")})
        base._ensureOutputDirs({"a": str(tmp_path / "a.parquet")})
        assert tmp_path.is_dir()

    def test_each_filter_may_live_in_its_own_directory(self, base, tmp_path):
        base._ensureOutputDirs({
            "a": str(tmp_path / "one" / "a.parquet"),
            "b": str(tmp_path / "two" / "b.parquet"),
        })
        assert (tmp_path / "one").is_dir() and (tmp_path / "two").is_dir()


# ---------------------------------------------------------------------------
# B312: _ensureOutputDirs cannot handle an output path with no directory
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestEnsureOutputDirsBareFilename:
    @pytest.mark.xfail(
        strict=True,
        reason="B312: _ensureOutputDirs does "
               "outputPath = os.path.dirname(outputFile) and then "
               "`if not os.path.isdir(outputPath): os.makedirs(outputPath)`. "
               "For a bare relative output name os.path.dirname returns '', "
               "os.path.isdir('') is False, and os.makedirs('') raises "
               "FileNotFoundError -- so the guard admits the one value "
               "makedirs cannot take. The directory in that case is the cwd, "
               "which already exists, so the correct behaviour is a no-op. "
               "See the consolidated findings issue.",
    )
    def test_an_output_name_relative_to_the_working_directory_needs_no_new_directory(
        self, base, tmp_path, monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        base._ensureOutputDirs({"a": "a.parquet"})

    def test_an_output_name_with_no_directory_part_raises_file_not_found(
        self, base, tmp_path, monkeypatch,
    ):
        """Characterisation of B312."""
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError):
            base._ensureOutputDirs({"a": "a.parquet"})

    def test_prefixing_the_same_name_with_dot_slash_makes_it_work(self, base, tmp_path, monkeypatch):
        """Characterisation of B312: the workaround."""
        monkeypatch.chdir(tmp_path)
        base._ensureOutputDirs({"a": "./a.parquet"})
        assert tmp_path.is_dir()


@pytest.mark.unit
class TestCollectTmpFiles:
    def _make(self, tmp_path, *names):
        for name in names:
            (tmp_path / name).write_text("x")

    def test_it_finds_every_numbered_temporary_block_of_one_filter(self, base, tmp_path):
        self._make(tmp_path, "tmp_slice1_000000.parquet", "tmp_slice1_000001.parquet")
        found = base._collectTmpFiles("slice1", str(tmp_path / "slice1.parquet"), "parquet")
        assert sorted(os.path.basename(f) for f in found) == [
            "tmp_slice1_000000.parquet", "tmp_slice1_000001.parquet",
        ]

    def test_temporary_blocks_of_other_filters_are_not_collected(self, base, tmp_path):
        self._make(tmp_path, "tmp_slice1_000000.parquet", "tmp_other_000000.parquet")
        found = base._collectTmpFiles("slice1", str(tmp_path / "slice1.parquet"), "parquet")
        assert [os.path.basename(f) for f in found] == ["tmp_slice1_000000.parquet"]

    def test_blocks_of_the_other_mesh_format_are_not_collected(self, base, tmp_path):
        self._make(tmp_path, "tmp_slice1_000000.parquet", "tmp_slice1_000001.zarr")
        found = base._collectTmpFiles("slice1", str(tmp_path / "slice1.zarr"), "zarr")
        assert [os.path.basename(f) for f in found] == ["tmp_slice1_000001.zarr"]

    def test_dots_in_the_filter_name_are_dashes_in_the_file_name(self, base, tmp_path):
        self._make(tmp_path, "tmp_a-b-c_000000.parquet")
        found = base._collectTmpFiles("a.b.c", str(tmp_path / "out.parquet"), "parquet")
        assert [os.path.basename(f) for f in found] == ["tmp_a-b-c_000000.parquet"]

    def test_it_searches_the_directory_of_the_final_output_file(self, base, tmp_path):
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        self._make(tmp_path, "tmp_slice1_000000.parquet")
        found = base._collectTmpFiles("slice1", str(elsewhere / "slice1.parquet"), "parquet")
        assert found == []

    def test_no_temporary_blocks_means_an_empty_list(self, base, tmp_path):
        assert base._collectTmpFiles("slice1", str(tmp_path / "s.parquet"), "parquet") == []


@pytest.mark.unit
class TestAtomicReplace:
    def test_the_staged_final_file_is_renamed_onto_the_target(self, base, tmp_path):
        target = tmp_path / "out.parquet"
        (tmp_path / "out.parquet.final").write_text("new")
        base._atomicReplace(str(target))
        assert target.read_text() == "new"
        assert not (tmp_path / "out.parquet.final").exists()

    def test_an_existing_target_file_is_removed_before_the_rename(self, base, tmp_path):
        target = tmp_path / "out.parquet"
        target.write_text("old")
        (tmp_path / "out.parquet.final").write_text("new")
        base._atomicReplace(str(target))
        assert target.read_text() == "new"

    def test_an_existing_target_directory_is_removed_before_the_rename(self, base, tmp_path):
        target = tmp_path / "out.parquet"
        (target / "part0").mkdir(parents=True)
        staged = tmp_path / "out.parquet.final"
        staged.mkdir()
        (staged / "part1").write_text("new")
        base._atomicReplace(str(target))
        assert sorted(p.name for p in target.iterdir()) == ["part1"]

    def test_a_missing_staged_file_is_an_error_not_a_silent_no_op(self, base, tmp_path):
        with pytest.raises(FileNotFoundError):
            base._atomicReplace(str(tmp_path / "out.parquet"))


@pytest.mark.unit
class TestCleanupTmpFiles:
    def test_temporary_files_are_removed(self, base, tmp_path):
        paths = []
        for name in ("t0", "t1"):
            path = tmp_path / name
            path.write_text("x")
            paths.append(str(path))
        base._cleanupTmpFiles(paths)
        assert list(tmp_path.iterdir()) == []

    def test_temporary_directories_are_removed_whole(self, base, tmp_path):
        block = tmp_path / "tmp_a_000000.parquet"
        (block / "part.0.parquet").mkdir(parents=True)
        base._cleanupTmpFiles([str(block)])
        assert not block.exists()

    def test_an_empty_list_of_temporaries_is_a_no_op(self, base, tmp_path):
        base._cleanupTmpFiles([])
        assert list(tmp_path.iterdir()) == []


# ---------------------------------------------------------------------------
# _writeTimeStepBlocks -- deep on the block arithmetic
# ---------------------------------------------------------------------------

def _blockRecorder(monkeypatch, cls, produced):
    """Replace readTimeSteps with a canned stream and record writeList calls."""
    written = []

    def fakeReadTimeSteps(self, datasourcenamedict, timelist, fieldnames, regularMesh):
        written.append(("read", dict(
            filters=sorted(datasourcenamedict), timelist=list(timelist),
            fieldnames=fieldnames, regularMesh=regularMesh,
        )))
        return iter(produced)

    def fakeWriteList(self, theList, blockID, filtersDict, regularMesh, fileExt):
        written.append(("write", blockID, list(theList), fileExt))

    monkeypatch.setattr(cls, "readTimeSteps", fakeReadTimeSteps)
    monkeypatch.setattr(cls, "writeList", fakeWriteList)
    return written


@pytest.mark.unit
class TestWriteTimeStepBlocks:
    def test_a_full_block_is_flushed_and_the_block_counter_advances(self, base, monkeypatch, pvmod):
        produced = [{"a": i} for i in range(4)]
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, produced)
        base._writeTimeStepBlocks({"a": "/o/a"}, [0, 1, 2, 3], None, False, "parquet", 2)
        blocks = [entry for entry in written if entry[0] == "write"]
        assert [entry[1] for entry in blocks] == [0, 1]
        assert [len(entry[2]) for entry in blocks] == [2, 2]

    def test_a_partial_trailing_block_is_still_flushed(self, base, monkeypatch, pvmod):
        produced = [{"a": i} for i in range(5)]
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, produced)
        base._writeTimeStepBlocks({"a": "/o/a"}, list(range(5)), None, False, "parquet", 2)
        blocks = [entry for entry in written if entry[0] == "write"]
        assert [entry[1] for entry in blocks] == [0, 1, 2]
        assert [len(entry[2]) for entry in blocks] == [2, 2, 1]

    def test_an_exactly_full_stream_does_not_flush_an_empty_extra_block(self, base, monkeypatch, pvmod):
        produced = [{"a": i} for i in range(4)]
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, produced)
        base._writeTimeStepBlocks({"a": "/o/a"}, list(range(4)), None, False, "parquet", 4)
        assert [entry[1] for entry in written if entry[0] == "write"] == [0]

    def test_no_timesteps_means_nothing_is_written_at_all(self, base, monkeypatch, pvmod):
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, [])
        base._writeTimeStepBlocks({"a": "/o/a"}, [], None, False, "parquet", 50)
        assert [entry for entry in written if entry[0] == "write"] == []

    def test_the_field_list_mesh_flag_and_time_list_reach_the_reader(self, base, monkeypatch, pvmod):
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, [])
        base._writeTimeStepBlocks({"a": "/o/a"}, [1.0, 2.0], ["T"], True, "zarr", 50)
        assert written[0] == ("read", dict(
            filters=["a"], timelist=[1.0, 2.0], fieldnames=["T"], regularMesh=True,
        ))

    def test_the_file_extension_is_passed_through_to_every_block_write(self, base, monkeypatch, pvmod):
        written = _blockRecorder(monkeypatch, pvmod.paraviewOpenFOAM, [{"a": 0}])
        base._writeTimeStepBlocks({"a": "/o/a"}, [0], None, True, "zarr", 50)
        assert [entry[3] for entry in written if entry[0] == "write"] == ["zarr"]


# ---------------------------------------------------------------------------
# writeList -- deep for parquet, shallow for zarr
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWriteListParquet:
    def test_it_writes_one_temporary_block_per_filter_next_to_its_output(self, base, tmp_path):
        theList = [{"slice1": pandas.DataFrame({"time": [0.0], "T": [1.0]})}]
        base.writeList(theList, 0, {"slice1": str(tmp_path / "slice1.parquet")}, False, "parquet")
        assert (tmp_path / "tmp_slice1_000000.parquet").exists()

    def test_the_block_number_is_zero_padded_to_six_digits(self, base, tmp_path):
        theList = [{"slice1": pandas.DataFrame({"time": [0.0]})}]
        base.writeList(theList, 42, {"slice1": str(tmp_path / "s.parquet")}, False, "parquet")
        assert (tmp_path / "tmp_slice1_000042.parquet").exists()

    def test_dots_in_the_filter_name_become_dashes_in_the_block_name(self, base, tmp_path):
        theList = [{"a.b": pandas.DataFrame({"time": [0.0]})}]
        base.writeList(theList, 0, {"a.b": str(tmp_path / "ab.parquet")}, False, "parquet")
        assert (tmp_path / "tmp_a-b_000000.parquet").exists()

    def test_every_timestep_in_the_block_ends_up_in_the_one_file(self, base, tmp_path):
        import dask.dataframe as dd

        theList = [
            {"s": pandas.DataFrame({"time": [0.0], "T": [1.0]})},
            {"s": pandas.DataFrame({"time": [1.0], "T": [2.0]})},
        ]
        base.writeList(theList, 0, {"s": str(tmp_path / "s.parquet")}, False, "parquet")
        result = dd.read_parquet(str(tmp_path / "tmp_s_000000.parquet")).compute()
        assert result.index.name == "time"
        assert sorted(result.index.tolist()) == [0.0, 1.0]
        assert sorted(result["T"].tolist()) == [1.0, 2.0]

    def test_a_filter_that_produced_several_blocks_per_timestep_is_flattened(self, base, tmp_path):
        import dask.dataframe as dd

        theList = [{"s": [
            pandas.DataFrame({"time": [0.0], "T": [1.0]}),
            pandas.DataFrame({"time": [0.0], "T": [2.0]}),
        ]}]
        base.writeList(theList, 0, {"s": str(tmp_path / "s.parquet")}, False, "parquet")
        result = dd.read_parquet(str(tmp_path / "tmp_s_000000.parquet")).compute()
        assert sorted(result["T"].tolist()) == [1.0, 2.0]

    def test_the_filters_to_write_are_taken_from_the_first_timestep(self, base, tmp_path):
        theList = [
            {"a": pandas.DataFrame({"time": [0.0]})},
            {"a": pandas.DataFrame({"time": [1.0]}), "b": pandas.DataFrame({"time": [1.0]})},
        ]
        base.writeList(theList, 0, {"a": str(tmp_path / "a.parquet")}, False, "parquet")
        assert [p.name for p in tmp_path.glob("tmp_*")] == ["tmp_a_000000.parquet"]


@pytest.mark.unit
class TestWriteListZarr:
    def test_a_regular_mesh_block_is_concatenated_along_time_and_written_as_zarr(
        self, base, pvmod, monkeypatch, tmp_path,
    ):
        """zarr is not installed; only the arguments writeList builds are checked."""
        concatenated = MagicMock()
        fake = fake_xarray(concat=MagicMock(return_value=concatenated))
        monkeypatch.setattr(pvmod, "xarray", fake)

        first = xarray.Dataset({"T": (("time",), [1.0])})
        second = xarray.Dataset({"T": (("time",), [2.0])})
        base.writeList([{"s": first}, {"s": second}], 3,
                       {"s": str(tmp_path / "s.zarr")}, True, "zarr")

        assert fake.concat.call_args.args[0] == [first, second]
        assert fake.concat.call_args.kwargs == {"dim": "time"}
        concatenated.to_zarr.assert_called_once_with(
            str(tmp_path / "tmp_s_000003.zarr"), mode="w",
        )


# ---------------------------------------------------------------------------
# _mergeParquet -- deep round trip;  _mergeZarr / _mergeToFinalOutput -- shallow
# ---------------------------------------------------------------------------

def _writeParquetBlock(path, frame):
    import dask.dataframe as dd

    dd.from_pandas(frame, npartitions=1).set_index("time").to_parquet(str(path))
    return str(path)


@pytest.mark.unit
class TestMergeParquet:
    def test_every_temporary_block_ends_up_in_the_final_output(self, base, tmp_path):
        import dask.dataframe as dd

        blocks = [
            _writeParquetBlock(tmp_path / "t0.parquet", pandas.DataFrame({"time": [0.0], "T": [1.0]})),
            _writeParquetBlock(tmp_path / "t1.parquet", pandas.DataFrame({"time": [1.0], "T": [2.0]})),
        ]
        output = str(tmp_path / "final.parquet")
        base._mergeParquet(output, blocks, append=False)
        result = dd.read_parquet(output + ".final").compute()
        assert result.index.name == "time"
        assert sorted(result["T"].tolist()) == [1.0, 2.0]

    def test_appending_keeps_the_rows_already_saved_in_the_output(self, base, tmp_path):
        import dask.dataframe as dd

        output = str(tmp_path / "final.parquet")
        _writeParquetBlock(output, pandas.DataFrame({"time": [9.0], "T": [99.0]}))
        blocks = [_writeParquetBlock(
            tmp_path / "t0.parquet", pandas.DataFrame({"time": [0.0], "T": [1.0]}),
        )]
        base._mergeParquet(output, blocks, append=True)
        result = dd.read_parquet(output + ".final").compute()
        assert sorted(result["T"].tolist()) == [1.0, 99.0]

    def test_not_appending_ignores_whatever_the_output_already_held(self, base, tmp_path):
        import dask.dataframe as dd

        output = str(tmp_path / "final.parquet")
        _writeParquetBlock(output, pandas.DataFrame({"time": [9.0], "T": [99.0]}))
        blocks = [_writeParquetBlock(
            tmp_path / "t0.parquet", pandas.DataFrame({"time": [0.0], "T": [1.0]}),
        )]
        base._mergeParquet(output, blocks, append=False)
        result = dd.read_parquet(output + ".final").compute()
        assert result["T"].tolist() == [1.0]

    def test_appending_with_no_previous_output_is_the_same_as_not_appending(self, base, tmp_path):
        import dask.dataframe as dd

        output = str(tmp_path / "final.parquet")
        blocks = [_writeParquetBlock(
            tmp_path / "t0.parquet", pandas.DataFrame({"time": [0.0], "T": [1.0]}),
        )]
        base._mergeParquet(output, blocks, append=True)
        assert dd.read_parquet(output + ".final").compute()["T"].tolist() == [1.0]


@pytest.mark.unit
class TestMergeZarr:
    def test_the_temporary_blocks_are_opened_lazily_with_the_zarr_engine(
        self, base, pvmod, monkeypatch, tmp_path,
    ):
        lazy = MagicMock()
        fake = fake_xarray(open_mfdataset=MagicMock(return_value=lazy))
        monkeypatch.setattr(pvmod, "xarray", fake)
        base._mergeZarr(str(tmp_path / "out.zarr"), ["b0.zarr", "b1.zarr"], append=False)
        fake.open_mfdataset.assert_called_once_with(
            ["b0.zarr", "b1.zarr"], chunks="auto", engine="zarr",
        )
        lazy.to_zarr.assert_called_once_with(str(tmp_path / "out.zarr") + ".final", mode="w")

    def test_appending_concatenates_the_old_data_and_sorts_by_time(
        self, base, pvmod, monkeypatch, tmp_path,
    ):
        output = tmp_path / "out.zarr"
        output.mkdir()
        lazy = MagicMock(name="new")
        old = MagicMock(name="old")
        fake = fake_xarray(
            open_mfdataset=MagicMock(side_effect=[lazy, old]),
            concat=MagicMock(),
        )
        monkeypatch.setattr(pvmod, "xarray", fake)
        base._mergeZarr(str(output), ["b0.zarr"], append=True)
        assert fake.concat.call_args.args[0] == [lazy, old]
        assert fake.concat.call_args.kwargs == {"dim": "time"}
        fake.concat.return_value.sortby.assert_called_once_with("time")
        fake.concat.return_value.sortby.return_value.to_zarr.assert_called_once_with(
            str(output) + ".final", mode="w",
        )

    def test_a_missing_output_is_not_treated_as_old_data_to_append(
        self, base, pvmod, monkeypatch, tmp_path,
    ):
        fake = fake_xarray(open_mfdataset=MagicMock())
        monkeypatch.setattr(pvmod, "xarray", fake)
        base._mergeZarr(str(tmp_path / "absent.zarr"), ["b0.zarr"], append=True)
        assert fake.open_mfdataset.call_count == 1
        assert fake.concat.call_count == 0

    def test_a_zarr_writer_that_refuses_the_chunking_is_retried_rechunked(
        self, base, pvmod, monkeypatch, tmp_path,
    ):
        lazy = MagicMock()
        lazy.to_zarr.side_effect = NotImplementedError("cannot write chunks")
        fake = fake_xarray(open_mfdataset=MagicMock(return_value=lazy))
        monkeypatch.setattr(pvmod, "xarray", fake)
        base._mergeZarr(str(tmp_path / "out.zarr"), ["b0.zarr"], append=False)
        lazy.chunk.assert_called_once_with("auto")
        lazy.chunk.return_value.to_zarr.assert_called_once_with(
            str(tmp_path / "out.zarr") + ".final", mode="w",
        )


@pytest.mark.unit
class TestMergeToFinalOutput:
    def test_a_regular_mesh_is_merged_as_zarr(self, base, pvmod, monkeypatch):
        calls = []
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "_mergeZarr",
                            lambda self, *a: calls.append(("zarr",) + a))
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "_mergeParquet",
                            lambda self, *a: calls.append(("parquet",) + a))
        base._mergeToFinalOutput("/o/out", ["b0"], True, False)
        assert calls == [("zarr", "/o/out", ["b0"], False)]

    def test_an_unstructured_mesh_is_merged_as_parquet(self, base, pvmod, monkeypatch):
        calls = []
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "_mergeZarr",
                            lambda self, *a: calls.append(("zarr",) + a))
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "_mergeParquet",
                            lambda self, *a: calls.append(("parquet",) + a))
        base._mergeToFinalOutput("/o/out", ["b0"], False, True)
        assert calls == [("parquet", "/o/out", ["b0"], True)]


# ---------------------------------------------------------------------------
# writeCase -- shallow: which helper, in which order, with which arguments
# ---------------------------------------------------------------------------

@pytest.fixture()
def writeCaseSpy(pvmod, monkeypatch):
    """Replace every writeCase helper with a recorder on the class."""
    calls = []
    cls = pvmod.paraviewOpenFOAM

    def recorder(name, result=None):
        def _call(self, *args, **kwargs):
            calls.append((name, args, kwargs))
            return result
        return _call

    monkeypatch.setattr(cls, "_removeOldOutputs", recorder("_removeOldOutputs"))
    monkeypatch.setattr(cls, "_ensureOutputDirs", recorder("_ensureOutputDirs"))
    monkeypatch.setattr(cls, "_resolveTimeList", recorder("_resolveTimeList", [0.0, 1.0]))
    monkeypatch.setattr(cls, "_writeTimeStepBlocks", recorder("_writeTimeStepBlocks"))
    monkeypatch.setattr(cls, "_collectTmpFiles", recorder("_collectTmpFiles", ["tmp0"]))
    monkeypatch.setattr(cls, "_mergeToFinalOutput", recorder("_mergeToFinalOutput"))
    monkeypatch.setattr(cls, "_atomicReplace", recorder("_atomicReplace"))
    monkeypatch.setattr(cls, "_cleanupTmpFiles", recorder("_cleanupTmpFiles"))
    return calls


FILTERS = {"slice1": "/out/slice1.parquet"}


@pytest.mark.unit
class TestWriteCase:
    def test_it_prepares_streams_then_merges_in_that_order(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False)
        assert [name for name, _, _ in writeCaseSpy] == [
            "_ensureOutputDirs", "_resolveTimeList", "_writeTimeStepBlocks",
            "_collectTmpFiles", "_mergeToFinalOutput", "_atomicReplace", "_cleanupTmpFiles",
        ]

    def test_overwriting_removes_the_old_outputs_first(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, overwrite=True)
        assert writeCaseSpy[0][0] == "_removeOldOutputs"
        assert writeCaseSpy[0][1] == (FILTERS,)

    def test_not_overwriting_leaves_the_old_outputs_in_place(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, overwrite=False)
        assert "_removeOldOutputs" not in [name for name, _, _ in writeCaseSpy]

    def test_a_regular_mesh_streams_zarr_blocks(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=True)
        stream = next(call for call in writeCaseSpy if call[0] == "_writeTimeStepBlocks")
        assert stream[1] == (FILTERS, [0.0, 1.0], None, True, "zarr", 50)

    def test_an_unstructured_mesh_streams_parquet_blocks(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, fieldnames=["T"], tsBlockNum=5)
        stream = next(call for call in writeCaseSpy if call[0] == "_writeTimeStepBlocks")
        assert stream[1] == (FILTERS, [0.0, 1.0], ["T"], False, "parquet", 5)

    def test_the_time_list_and_latest_flag_are_resolved_before_streaming(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, timeList=[3.0], latestTimestamp=True)
        resolve = next(call for call in writeCaseSpy if call[0] == "_resolveTimeList")
        assert resolve[1] == ([3.0], True)

    def test_appending_is_the_negation_of_overwriting(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, overwrite=False)
        merge = next(call for call in writeCaseSpy if call[0] == "_mergeToFinalOutput")
        assert merge[1] == ("/out/slice1.parquet", ["tmp0"], False, True)

    def test_overwriting_switches_the_merge_out_of_append_mode(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False, overwrite=True)
        merge = next(call for call in writeCaseSpy if call[0] == "_mergeToFinalOutput")
        assert merge[1][3] is False

    def test_each_filter_is_merged_replaced_and_cleaned_independently(self, base, writeCaseSpy):
        filters = {"a": "/out/a.parquet", "b": "/out/b.parquet"}
        base.writeCase(filters, regularMesh=False)
        tail = [name for name, _, _ in writeCaseSpy][2:]
        assert tail == [
            "_writeTimeStepBlocks",
            "_collectTmpFiles", "_mergeToFinalOutput", "_atomicReplace", "_cleanupTmpFiles",
            "_collectTmpFiles", "_mergeToFinalOutput", "_atomicReplace", "_cleanupTmpFiles",
        ]

    def test_the_temporary_blocks_that_were_merged_are_the_ones_cleaned_up(self, base, writeCaseSpy):
        base.writeCase(FILTERS, regularMesh=False)
        cleanup = next(call for call in writeCaseSpy if call[0] == "_cleanupTmpFiles")
        assert cleanup[1] == (["tmp0"],)


# ---------------------------------------------------------------------------
# writeCase end to end, real filesystem, unstructured only (zarr is absent)
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWriteCaseEndToEnd:
    def test_an_unstructured_case_lands_as_one_parquet_output_per_filter(
        self, pv, monkeypatch, tmp_path,
    ):
        import dask.dataframe as dd

        obj = pv.cls(casePath="/c")
        obj.reader = FakeReader(timesteps=[0.0, 1.0])
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(
            pv.cls, "_readTimeStep",
            lambda self, datasource, timeslice, fieldnames, regularMesh: pandas.DataFrame(
                {"time": [timeslice], "T": [timeslice * 10]},
            ),
        )
        output = tmp_path / "results" / "slice1.parquet"
        obj.writeCase({"slice1": str(output)}, regularMesh=False, overwrite=True)

        assert output.exists()
        assert list(tmp_path.glob("results/tmp_*")) == []
        result = dd.read_parquet(str(output)).compute()
        assert sorted(result.index.tolist()) == [0.0, 1.0]
        assert sorted(result["T"].tolist()) == [0.0, 10.0]

    def test_a_second_non_overwriting_run_appends_to_the_existing_output(
        self, pv, monkeypatch, tmp_path,
    ):
        import dask.dataframe as dd

        obj = pv.cls(casePath="/c")
        monkeypatch.setattr(pv.pvsimple, "FindSource", lambda name: FakeProxy())
        monkeypatch.setattr(
            pv.cls, "_readTimeStep",
            lambda self, datasource, timeslice, fieldnames, regularMesh: pandas.DataFrame(
                {"time": [timeslice], "T": [timeslice * 10]},
            ),
        )
        output = tmp_path / "results" / "slice1.parquet"
        obj.writeCase({"slice1": str(output)}, regularMesh=False, timeList=[0.0], overwrite=True)
        obj.writeCase({"slice1": str(output)}, regularMesh=False, timeList=[1.0], overwrite=False)
        result = dd.read_parquet(str(output)).compute()
        assert sorted(result.index.tolist()) == [0.0, 1.0]


# ---------------------------------------------------------------------------
# Deprecated aliases -- shallow: pure forwarding
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestDeprecatedWriteAliases:
    @pytest.fixture()
    def spy(self, pvmod, monkeypatch):
        recorded = []
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "writeCase",
                            lambda self, *args, **kwargs: recorded.append((args, kwargs)))
        return recorded

    def test_write_netcdf_asks_write_case_for_a_regular_mesh(self, base, spy):
        base.write_netcdf({"a": "/o/a"}, [1.0], ["T"], 7, True)
        assert spy == [((({"a": "/o/a"}), True, [1.0], ["T"], 7, True), {})]

    def test_write_parquet_asks_write_case_for_an_unstructured_mesh(self, base, spy):
        base.write_parquet({"a": "/o/a"}, [1.0], ["T"], 7, True)
        assert spy == [((({"a": "/o/a"}), False, [1.0], ["T"], 7, True), {})]

    def test_the_alias_append_flag_never_reaches_write_case(self, base, spy):
        """writeCase derives append from overwrite, so the old flag is dead."""
        base.write_parquet({"a": "/o/a"}, append=True, overwrite=False)
        args, kwargs = spy[0]
        assert True not in args[1:]
        assert kwargs == {}

    def test_write_parquets_filter_list_argument_is_ignored(self, base, spy):
        base.write_parquet({"a": "/o/a"}, filterList=["only-this-one"])
        args, kwargs = spy[0]
        assert ["only-this-one"] not in args
        assert kwargs == {}


@pytest.mark.unit
class TestDeprecatedReadAliases:
    @pytest.fixture()
    def spy(self, pvmod, monkeypatch):
        recorded = []
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "readTimeSteps",
                            lambda self, *args, **kwargs: recorded.append((args, kwargs)))
        return recorded

    @pytest.mark.parametrize("alias", ["to_pandas", "to_dataFrame"])
    def test_the_pandas_aliases_ask_for_an_unstructured_read(self, base, spy, alias):
        getattr(base, alias)({"a": "/o/a"}, [1.0], ["T"])
        assert spy == [((({"a": "/o/a"}), [1.0], ["T"]), {"regularMesh": False})]

    @pytest.mark.parametrize("alias", ["to_xarray", "to_dataArray"])
    def test_the_xarray_aliases_ask_for_a_regular_mesh_read(self, base, spy, alias):
        getattr(base, alias)({"a": "/o/a"}, [1.0], ["T"])
        assert spy == [((({"a": "/o/a"}), [1.0], ["T"]), {"regularMesh": True})]

    @pytest.mark.parametrize(
        "alias", ["to_pandas", "to_dataFrame", "to_xarray", "to_dataArray"],
    )
    def test_every_read_alias_returns_whatever_read_time_steps_returned(
        self, base, pvmod, monkeypatch, alias,
    ):
        monkeypatch.setattr(pvmod.paraviewOpenFOAM, "readTimeSteps",
                            lambda self, *args, **kwargs: "the-generator")
        assert getattr(base, alias)({"a": "/o/a"}, [1.0]) == "the-generator"

    @pytest.mark.parametrize(
        "alias", ["to_pandas", "to_dataFrame", "to_xarray", "to_dataArray"],
    )
    def test_every_read_alias_is_marked_deprecated(self, base, alias):
        with pytest.deprecated_call():
            getattr(base, alias)({"a": "/o/a"}, [])
