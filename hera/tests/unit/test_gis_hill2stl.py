"""hill2stl: generating an ASCII STL mesh from a height function.

Pure numpy geometry -- no GIS toolkit, no DB.
"""
import io

import numpy
import pytest

from hera.measurements.GIS.raster.hill2stl import compute_normal, function, generate_solid_stl, write_triangle


@pytest.mark.unit
class TestComputeNormal:
    def test_the_normal_of_the_xy_plane_points_along_z(self):
        normal = compute_normal((0, 0, 0), (1, 0, 0), (0, 1, 0))
        assert numpy.allclose(normal, [0, 0, 1])

    def test_the_normal_is_a_unit_vector(self):
        normal = compute_normal((0, 0, 0), (2, 0, 0), (0, 3, 0))
        assert numpy.linalg.norm(normal) == pytest.approx(1.0)

    def test_the_normal_is_perpendicular_to_both_triangle_edges(self):
        v1, v2, v3 = (0, 0, 0), (1, 0, 1), (0, 1, 2)
        normal = compute_normal(v1, v2, v3)
        edge1 = numpy.array(v2) - numpy.array(v1)
        edge2 = numpy.array(v3) - numpy.array(v1)
        assert numpy.dot(normal, edge1) == pytest.approx(0.0, abs=1e-9)
        assert numpy.dot(normal, edge2) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.unit
class TestWriteTriangle:
    def test_it_writes_a_well_formed_ascii_stl_facet(self):
        buf = io.StringIO()
        write_triangle(buf, (0, 0, 0), (1, 0, 0), (0, 1, 0))
        text = buf.getvalue()
        assert text.startswith("  facet normal ")
        assert "outer loop" in text
        assert text.count("vertex") == 3
        assert text.rstrip().endswith("endfacet")

    def test_the_facet_normal_line_matches_compute_normal(self):
        buf = io.StringIO()
        v1, v2, v3 = (0, 0, 0), (1, 0, 0), (0, 1, 0)
        write_triangle(buf, v1, v2, v3)
        normal = compute_normal(v1, v2, v3)
        expected = f"  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n"
        assert buf.getvalue().splitlines()[0] + "\n" == expected


@pytest.mark.unit
class TestFunction:
    def test_it_is_the_documented_cos_squared_surface(self):
        assert function(0.0, 0.0) == pytest.approx(4.0)

    def test_it_never_drops_below_three(self):
        xs = numpy.linspace(-10, 10, 50)
        assert numpy.all(function(xs, xs) >= 3.0)


@pytest.mark.unit
class TestGenerateSolidStl:
    def test_it_writes_a_valid_solid_with_matching_facet_counts(self, tmp_path):
        target = str(tmp_path / "out.stl")
        generate_solid_stl(function, x_range=(-1, 1), y_range=(-1, 1), resolution=4, filename=target)
        text = open(target).read()
        assert text.startswith("solid solid_surface\n")
        assert text.rstrip().endswith("endsolid solid_surface")
        assert text.count("facet normal") == text.count("endfacet")
        assert text.count("facet normal") > 0

    def test_the_flat_base_sits_below_the_minimum_surface_height(self, tmp_path):
        target = str(tmp_path / "out.stl")
        generate_solid_stl(function, x_range=(-1, 1), y_range=(-1, 1), resolution=4, filename=target)
        text = open(target).read()
        vertex_zs = [float(line.split()[-1]) for line in text.splitlines() if line.strip().startswith("vertex")]
        # function() >= 3 everywhere, and the base is z_min = min(Z) - 0.5
        assert min(vertex_zs) < 3.0
