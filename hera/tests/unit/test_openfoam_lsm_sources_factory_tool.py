"""openFoam/lagrangian/LSM/sourcesFactoryTool.py: same geometry-source
pattern as absractStochasticLagrangianSolver_toolkitExtension (batch 23),
duplicated here as its own factory class."""
import numpy
import pytest

from hera.simulations.openFoam.lagrangian.LSM.sourcesFactoryTool import sourcesFactoryTool


@pytest.fixture()
def factory():
    return sourcesFactoryTool()


@pytest.mark.unit
class TestSourcesTypeList:
    def test_it_lists_every_makesource_suffix(self, factory):
        assert set(factory.sourcesTypeList) == {"Point", "Circle", "Sphere", "Cylinder", "Rectangle", "Cube"}


@pytest.mark.unit
class TestMakeSourceDispatch:
    def test_it_dispatches_to_the_matching_method(self, factory):
        df = factory.makeSource(0, 0, 0, 5, type="Point")
        assert len(df) == 5
        assert (df["x"] == 0).all()

    def test_an_unknown_type_raises(self, factory):
        with pytest.raises(ValueError, match="type must be"):
            factory.makeSource(0, 0, 0, 5, type="NoSuchShape")


@pytest.mark.unit
class TestMakeSourcePoint:
    def test_every_particle_sits_at_the_exact_location(self, factory):
        df = factory.makeSource_Point(1, 2, 3, 5)
        assert (df["x"] == 1).all() and (df["y"] == 2).all() and (df["z"] == 3).all()


@pytest.mark.unit
class TestMakeSourceCircle:
    def test_all_particles_stay_within_the_radius(self, factory):
        df = factory.makeSource_Circle(0, 0, 5, 200, radius=10)
        distances = numpy.hypot(df["x"], df["y"])
        assert (distances <= 10 + 1e-9).all()

    def test_z_is_constant(self, factory):
        df = factory.makeSource_Circle(0, 0, 5, 10, radius=10)
        assert (df["z"] == 5).all()


@pytest.mark.unit
class TestMakeSourceSphere:
    def test_all_particles_stay_within_the_radius(self, factory):
        df = factory.makeSource_Sphere(0, 0, 0, 200, radius=10)
        distances = numpy.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)
        assert (distances <= 10 + 1e-9).all()


@pytest.mark.unit
class TestMakeSourceRectangleAndCube:
    def test_rectangle_stays_within_its_bounds(self, factory):
        df = factory.makeSource_Rectangle(0, 0, 0, 200, lengthX=10, lengthY=4, rotateAngle=0)
        assert (df["x"].abs() <= 5 + 1e-9).all()
        assert (df["y"].abs() <= 2 + 1e-9).all()

    def test_cube_stays_within_its_bounds(self, factory):
        df = factory.makeSource_Cube(0, 0, 0, 200, lengthX=10, lengthY=4, lengthZ=6, rotateAngle=0)
        assert (df["x"].abs() <= 5 + 1e-9).all()
        assert (df["z"].abs() <= 3 + 1e-9).all()


@pytest.mark.unit
class TestMakeSourceCylinder:
    def test_radial_distance_stays_within_the_radius(self, factory):
        df = factory.makeSource_Cylinder(0, 0, 0, 200, radius=5, height=20)
        radial = numpy.hypot(df["x"], df["y"])
        assert (radial <= 5 + 1e-9).all()

    def test_height_stays_within_plus_minus_half_height(self, factory):
        df = factory.makeSource_Cylinder(0, 0, 100, 200, radius=5, height=20)
        assert (df["z"] >= 90 - 1e-9).all()
        assert (df["z"] <= 110 + 1e-9).all()
