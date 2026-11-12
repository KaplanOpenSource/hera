"""utils/freeCAD.py: ``getObjFileBoundaries``.

The function loads a mesh file through FreeCAD's ``Mesh`` module, asks the
resulting document for its objects, and reduces their bounding boxes to a
single corner pair.  The documented behaviour that is actually testable
without FreeCAD installed is the arithmetic on top of those boxes:

* the maxima are the element-wise maximum over every object's box and the
  minima the element-wise minimum, so the result bounds the *union* of the
  meshes rather than any one of them;
* every value is divided by 1000, i.e. FreeCAD's millimetres become the
  metres the rest of hera works in.

``Mesh`` and the document are supplied here as fakes -- ``Mesh.open`` and
``FreeCAD.getDocument("Unnamed")`` are the module's only two calls into
FreeCAD, and ``findObjects()`` its only call on the document.

One defect surfaced:

* B253: the module's ``try: import FreeCAD; import Mesh except
  ImportError`` is downgraded to ``warnings.warn`` so that the rest of
  ``hera.utils`` stays importable.  Nothing then guards the use, so
  calling ``getObjFileBoundaries`` without FreeCAD raises ``NameError:
  name 'Mesh' is not defined`` -- which says nothing about FreeCAD, the
  package the warning at import time told the user to install.  In this
  test environment the shape is the same but sharper: the unit layer's
  conftest stubs ``FreeCAD`` and not ``Mesh``, so ``FreeCAD`` is bound
  while ``Mesh`` is not, and the same NameError is raised on the first
  line of the function body.
"""
import sys
import types

import pytest

import hera.utils.freeCAD as freeCADModule


class _BoundBox:
    """The six attributes the function reads off a Mesh.BoundBox."""

    def __init__(self, xMin, yMin, zMin, xMax, yMax, zMax):
        self.XMin, self.YMin, self.ZMin = xMin, yMin, zMin
        self.XMax, self.YMax, self.ZMax = xMax, yMax, zMax


class _Document:
    def __init__(self, boundBoxes):
        self._objects = [
            types.SimpleNamespace(Mesh=types.SimpleNamespace(BoundBox=box))
            for box in boundBoxes
        ]

    def findObjects(self):
        return self._objects


@pytest.fixture()
def freeCADWith(monkeypatch):
    """Bind a fake ``Mesh`` and a fake document, and record the file opened.

    The module-level name is what the function body reads, so both the
    sys.modules entry and the module attribute have to be set -- the latter
    with raising=False, because the failed import left no attribute to
    replace.
    """
    opened = []

    def _install(boundBoxes):
        fakeMesh = types.SimpleNamespace(open=lambda fileName: opened.append(fileName))
        monkeypatch.setitem(sys.modules, "Mesh", fakeMesh)
        monkeypatch.setattr(freeCADModule, "Mesh", fakeMesh, raising=False)
        monkeypatch.setattr(
            freeCADModule.FreeCAD,
            "getDocument",
            lambda name: _Document(boundBoxes),
            raising=False,
        )
        return opened

    return _install


ONE_BOX = [_BoundBox(-1000.0, 0.0, 500.0, 1000.0, 2000.0, 3000.0)]


@pytest.mark.unit
class TestObjFileBoundaries:
    def test_it_reports_all_six_corners(self, freeCADWith):
        freeCADWith(ONE_BOX)
        corners = freeCADModule.getObjFileBoundaries("model.obj")
        assert set(corners) == {"XMax", "YMax", "ZMax", "XMin", "YMin", "ZMin"}

    def test_millimetres_become_metres(self, freeCADWith):
        """FreeCAD works in mm; the documented /1000 makes the result metres."""
        freeCADWith(ONE_BOX)
        corners = freeCADModule.getObjFileBoundaries("model.obj")
        assert corners == pytest.approx(
            {"XMin": -1.0, "YMin": 0.0, "ZMin": 0.5,
             "XMax": 1.0, "YMax": 2.0, "ZMax": 3.0}
        )

    def test_the_requested_file_is_the_one_opened(self, freeCADWith):
        opened = freeCADWith(ONE_BOX)
        freeCADModule.getObjFileBoundaries("/tmp/terrain.stl")
        assert opened == ["/tmp/terrain.stl"]

    def test_it_reads_the_unnamed_document(self, freeCADWith, monkeypatch):
        """Mesh.open puts the mesh into a document called 'Unnamed', which is
        what the function then asks FreeCAD for."""
        freeCADWith(ONE_BOX)
        requested = []
        monkeypatch.setattr(
            freeCADModule.FreeCAD,
            "getDocument",
            lambda name: requested.append(name) or _Document(ONE_BOX),
            raising=False,
        )
        freeCADModule.getObjFileBoundaries("model.obj")
        assert requested == ["Unnamed"]

    def test_the_maxima_never_fall_below_the_minima(self, freeCADWith):
        freeCADWith(ONE_BOX)
        corners = freeCADModule.getObjFileBoundaries("model.obj")
        for axis in "XYZ":
            assert corners[f"{axis}Min"] <= corners[f"{axis}Max"]

    def test_several_objects_are_bounded_as_a_union(self, freeCADWith):
        """The maxima are an element-wise numpy.max over the objects and the
        minima an element-wise numpy.min, so the box encloses every mesh --
        which is the only sensible reading of "its boundaries" for a file
        holding more than one solid."""
        freeCADWith(
            [
                _BoundBox(0.0, 0.0, 0.0, 1000.0, 1000.0, 1000.0),
                _BoundBox(-5000.0, 2000.0, -1000.0, -1000.0, 4000.0, 500.0),
            ]
        )
        corners = freeCADModule.getObjFileBoundaries("two.obj")
        assert corners == pytest.approx(
            {"XMin": -5.0, "YMin": 0.0, "ZMin": -1.0,
             "XMax": 1.0, "YMax": 4.0, "ZMax": 1.0}
        )

    def test_a_contained_object_does_not_change_the_bounds(self, freeCADWith):
        """A union is idempotent over nested boxes."""
        outer = _BoundBox(-1000.0, -1000.0, -1000.0, 1000.0, 1000.0, 1000.0)
        inner = _BoundBox(-10.0, -10.0, -10.0, 10.0, 10.0, 10.0)
        freeCADWith([outer, inner])
        corners = freeCADModule.getObjFileBoundaries("nested.obj")
        assert corners == pytest.approx(
            {"XMin": -1.0, "YMin": -1.0, "ZMin": -1.0,
             "XMax": 1.0, "YMax": 1.0, "ZMax": 1.0}
        )

    def test_negative_coordinates_survive_the_conversion(self, freeCADWith):
        """Dividing by 1000 must not be confused with taking a magnitude:
        an ITM-shifted model sits entirely at negative coordinates."""
        freeCADWith([_BoundBox(-9000.0, -8000.0, -7000.0, -6000.0, -5000.0, -4000.0)])
        corners = freeCADModule.getObjFileBoundaries("shifted.obj")
        assert corners["XMin"] == pytest.approx(-9.0)
        assert corners["XMax"] == pytest.approx(-6.0)

    def test_the_extent_scales_with_the_model(self, freeCADWith):
        """A model ten times larger has ten times the extent, and the /1000
        is a pure scaling with no offset."""
        freeCADWith([_BoundBox(0.0, 0.0, 0.0, 1000.0, 1000.0, 1000.0)])
        small = freeCADModule.getObjFileBoundaries("small.obj")
        freeCADWith([_BoundBox(0.0, 0.0, 0.0, 10000.0, 10000.0, 10000.0)])
        large = freeCADModule.getObjFileBoundaries("large.obj")
        assert large["XMax"] / small["XMax"] == pytest.approx(10.0)

    def test_a_file_with_no_objects_cannot_be_bounded(self, freeCADWith):
        """numpy.max over an empty sequence raises; there is no box to
        report for an empty document."""
        freeCADWith([])
        with pytest.raises(ValueError):
            freeCADModule.getObjFileBoundaries("empty.obj")


@pytest.mark.unit
class TestMissingFreeCAD:
    """B253: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B253: the module downgrades `import FreeCAD; import Mesh` "
               "failures to warnings.warn('freecad is not installed. some "
               "features will not work.') and then uses both names unguarded. "
               "Calling getObjFileBoundaries without FreeCAD raises NameError: "
               "name 'Mesh' is not defined, which names neither FreeCAD nor the "
               "remedy the import-time warning offered. The unit layer's conftest "
               "stubs FreeCAD but not Mesh, so the same NameError appears on the "
               "function's first line. See the consolidated findings issue.",
    )
    def test_a_missing_freecad_is_reported_as_such(self):
        with pytest.raises(ImportError, match="freecad"):
            freeCADModule.getObjFileBoundaries("model.obj")

    def test_it_currently_raises_a_bare_nameerror(self):
        """Characterisation of B253."""
        if hasattr(freeCADModule, "Mesh"):
            pytest.skip("a real FreeCAD/Mesh install is present")
        with pytest.raises(NameError, match="Mesh"):
            freeCADModule.getObjFileBoundaries("model.obj")

    def test_the_import_failure_is_only_a_warning(self):
        """Characterisation of B253's mechanism: importing the module
        succeeds and says so, which is why the failure surfaces late."""
        import inspect

        source = inspect.getsource(freeCADModule)
        assert "except ImportError:" in source
        assert 'warnings.warn("freecad is not installed' in source
        # The raise the warning replaced is still there, commented out.
        assert "#raise ImportError(" in source

    def test_the_two_names_are_bound_independently(self):
        """Characterisation of B253's blast radius: `import FreeCAD` and
        `import Mesh` are separate statements inside one try, so a partial
        install binds the first and leaves the second missing -- exactly what
        the conftest stub reproduces."""
        assert hasattr(freeCADModule, "FreeCAD")
        assert not hasattr(freeCADModule, "Mesh")
