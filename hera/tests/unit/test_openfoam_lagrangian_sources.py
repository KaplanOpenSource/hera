"""absractStochasticLagrangianSolver_toolkitExtension: the pure-geometry
particle-source generators, and the file writer built on top of them.

Randomised distributions are checked against their bounding geometry (a
circle's particles all fall within its radius, a cylinder's within its
radius and half-height, etc.) rather than exact values.
"""
import os

import numpy
import pandas
import pytest

from hera.simulations.openFoam.lagrangian.abstractLagrangianSolver import (
    absractStochasticLagrangianSolver_toolkitExtension,
)


@pytest.fixture()
def ext():
    return absractStochasticLagrangianSolver_toolkitExtension.__new__(
        absractStochasticLagrangianSolver_toolkitExtension
    )


@pytest.mark.unit
class TestSourcesTypeList:
    def test_it_lists_every_makesource_suffix(self, ext):
        assert set(ext.sourcesTypeList) == {"Point", "Circle", "Sphere", "Cylinder", "Rectangle", "Cube"}


@pytest.mark.unit
class TestMakeSourcePoint:
    def test_every_particle_sits_at_the_exact_location(self, ext):
        df = ext.makeSource_Point(1, 2, 3, nParticles=5)
        assert len(df) == 5
        assert (df["x"] == 1).all() and (df["y"] == 2).all() and (df["z"] == 3).all()


@pytest.mark.unit
class TestMakeSourceCircle:
    def test_all_particles_stay_within_the_radius_on_the_xy_plane(self, ext):
        df = ext.makeSource_Circle(0, 0, 5, nParticles=200, radius=10)
        distances = numpy.hypot(df["x"], df["y"])
        assert (distances <= 10 + 1e-9).all()

    def test_z_is_constant_at_the_source_height(self, ext):
        df = ext.makeSource_Circle(0, 0, 5, nParticles=10, radius=10)
        assert (df["z"] == 5).all()


@pytest.mark.unit
class TestMakeSourceSphere:
    def test_all_particles_stay_within_the_radius_in_3d(self, ext):
        df = ext.makeSource_Sphere(0, 0, 0, nParticles=200, radius=10)
        distances = numpy.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)
        assert (distances <= 10 + 1e-9).all()


@pytest.mark.unit
class TestMakeSourceCylinder:
    def test_radial_distance_stays_within_the_radius(self, ext):
        df = ext.makeSource_Cylinder(0, 0, 0, nParticles=200, radius=5, height=20)
        radial = numpy.hypot(df["x"], df["y"])
        assert (radial <= 5 + 1e-9).all()

    def test_height_stays_within_plus_minus_half_height(self, ext):
        df = ext.makeSource_Cylinder(0, 0, 100, nParticles=200, radius=5, height=20)
        assert (df["z"] >= 90 - 1e-9).all()
        assert (df["z"] <= 110 + 1e-9).all()


@pytest.mark.unit
class TestMakeSourceRectangle:
    def test_unrotated_particles_stay_within_the_rectangle_bounds(self, ext):
        df = ext.makeSource_Rectangle(0, 0, 0, nParticles=200, lengthX=10, lengthY=4, rotateAngle=0)
        assert (df["x"].abs() <= 5 + 1e-9).all()
        assert (df["y"].abs() <= 2 + 1e-9).all()

    def test_z_is_constant(self, ext):
        df = ext.makeSource_Rectangle(0, 0, 7, nParticles=10, lengthX=10, lengthY=4)
        assert (df["z"] == 7).all()


@pytest.mark.unit
class TestMakeSourceCube:
    def test_unrotated_particles_stay_within_the_cuboid_bounds(self, ext):
        df = ext.makeSource_Cube(0, 0, 0, nParticles=200, lengthX=10, lengthY=4, lengthZ=6, rotateAngle=0)
        assert (df["x"].abs() <= 5 + 1e-9).all()
        assert (df["y"].abs() <= 2 + 1e-9).all()
        assert (df["z"].abs() <= 3 + 1e-9).all()


@pytest.mark.unit
class TestWriteParticlePositionFile:
    def test_an_unknown_source_type_raises(self, ext, tmp_path):
        with pytest.raises(ValueError, match="type must be"):
            ext.writeParticlePositionFile(0, 0, 0, nParticles=5, dispersionName=str(tmp_path), type="NoSuchShape")

    def test_it_writes_one_coordinate_line_per_particle(self, ext, tmp_path):
        (tmp_path / "constant").mkdir()
        ext.writeParticlePositionFile(1, 2, 3, nParticles=3, dispersionName=str(tmp_path), type="Point")
        content = (tmp_path / "constant" / "kinematicCloudPositions").read_text()
        assert content.count("(1 2 3)") == 3
        assert content.rstrip().endswith(")")

    def test_the_header_declares_the_particle_count(self, ext, tmp_path):
        (tmp_path / "constant").mkdir()
        ext.writeParticlePositionFile(0, 0, 0, nParticles=4, dispersionName=str(tmp_path), type="Point")
        content = (tmp_path / "constant" / "kinematicCloudPositions").read_text()
        assert "\n4\n(\n" in content
