"""Spreading Gaussian puffs onto a regular mesh.

The physics to hold onto is mass conservation: a normalised 2-D Gaussian
integrated over the mesh must return the released quantity Q.  The tests below
recover sigma and total mass from the produced field, so they check the maths
rather than snapshotting an array.
"""
import numpy as np
import pandas as pd
import pytest

from hera.simulations.gaussian.MeshUtils import GaussianIntegrationToMesh, GaussianToMesh
from hera.utils.unitHandler import ureg


def recovered_std(coordinates, weights):
    """Standard deviation of a discrete distribution, for checking sigma."""
    normalised = weights / np.sum(weights)
    mean = np.sum(coordinates * normalised)
    return np.sqrt(np.sum(normalised * (coordinates - mean) ** 2))


def one_puff(x=0.0, y=0.0, sigmaX=5.0, sigmaY=None, quantity=1.0):
    """A single puff as the mesher expects it: a Series with named columns."""
    fields = {"x": x, "y": y, "sigmaXCorrected": sigmaX, "Q": quantity}
    if sigmaY is not None:
        fields["sigmaYCorrected"] = sigmaY
    return pd.Series(fields)


@pytest.fixture()
def mesher():
    return GaussianToMesh(dxdy=1.0)


@pytest.mark.unit
class TestMassConservation:
    def test_a_normalised_puff_integrates_to_its_quantity(self, mesher):
        field = mesher.gaussianToMesh(one_puff(quantity=1.0))
        assert float(field.sum()) * 1.0 * 1.0 == pytest.approx(1.0, rel=1e-6)

    @pytest.mark.parametrize("quantity", [0.5, 2.0, 10.0])
    def test_the_integral_scales_linearly_with_the_quantity(self, mesher, quantity):
        field = mesher.gaussianToMesh(one_puff(quantity=quantity))
        assert float(field.sum()) * 1.0 * 1.0 == pytest.approx(quantity, rel=1e-6)

    def test_two_puffs_carry_twice_the_mass(self, mesher):
        frame = pd.DataFrame([one_puff(x=0.0), one_puff(x=20.0)]).reset_index(drop=True)
        field = mesher.gaussianToMesh(frame)
        assert float(field.sum()) * 1.0 * 1.0 == pytest.approx(2.0, rel=1e-3)

    def test_the_field_is_never_negative(self, mesher):
        field = mesher.gaussianToMesh(one_puff())
        assert float(field.min()) >= 0


@pytest.mark.unit
class TestShape:
    def test_the_peak_sits_at_the_puff_centre(self, mesher):
        field = mesher.gaussianToMesh(one_puff(x=10.0, y=-5.0))
        peak = field.where(field == field.max(), drop=True)
        assert float(peak.coords["x"][0]) == pytest.approx(10.0, abs=1.0)
        assert float(peak.coords["y"][0]) == pytest.approx(-5.0, abs=1.0)

    @pytest.mark.parametrize("sigmaX", [2.0, 5.0, 10.0])
    def test_the_recovered_alongwind_sigma_matches_the_input(self, mesher, sigmaX):
        field = mesher.gaussianToMesh(one_puff(sigmaX=sigmaX))
        marginal = field.sum(dim="y").values
        assert recovered_std(field.coords["x"].values, marginal) == pytest.approx(
            sigmaX, rel=1e-2
        )

    def test_the_domain_spans_six_sigma_each_way(self, mesher):
        """Documented: [min - 6 sigma, max + 6 sigma]."""
        assert mesher.sigmashifts == 6.0
        coordinates, _ = mesher._defineCoordinates(one_puff(x=0.0, sigmaX=5.0))
        assert coordinates.min() == pytest.approx(-30.0)
        assert coordinates.max() < 30.0

    def test_a_wider_puff_needs_a_wider_domain(self, mesher):
        narrow, _ = mesher._defineCoordinates(one_puff(sigmaX=2.0))
        wide, _ = mesher._defineCoordinates(one_puff(sigmaX=20.0))
        assert wide.max() - wide.min() > narrow.max() - narrow.min()

    def test_the_cell_size_defaults_to_a_fifth_of_the_smallest_sigma(self):
        """Documented: number of points is domain width / (min sigma / 5)."""
        automatic = GaussianToMesh()
        coordinates, _ = automatic._defineCoordinates(one_puff(sigmaX=5.0))
        assert coordinates[1] - coordinates[0] == pytest.approx(1.0)


@pytest.mark.unit
class TestIntegrationVariant:
    """GaussianIntegrationToMesh uses erf differences instead of point values.

    Its docstring calls itself "more accurate than estimating the gaussian",
    which invites swapping it in for the base class.  It cannot be swapped in
    freely -- see B29 at the bottom of this class.
    """

    def test_it_also_conserves_mass(self):
        mesher = GaussianIntegrationToMesh(dxdy=1.0)
        field = mesher.gaussianToMesh(one_puff(quantity=1.0 * ureg.kg))
        assert float(field.sum()) * 1.0 * 1.0 == pytest.approx(1.0, rel=1e-3)

    def test_it_agrees_with_the_point_evaluation_on_a_fine_mesh(self):
        """The two differ by discretisation error, which shrinks with the cell."""
        point = GaussianToMesh(dxdy=0.25).gaussianToMesh(one_puff())
        integrated = GaussianIntegrationToMesh(dxdy=0.25).gaussianToMesh(
            one_puff(quantity=1.0 * ureg.kg)
        )
        assert float(point.sum()) * 0.0625 == pytest.approx(
            float(integrated.sum()) * 0.0625, rel=1e-2
        )

    def test_it_is_a_subclass_that_only_replaces_the_kernel(self):
        assert issubclass(GaussianIntegrationToMesh, GaussianToMesh)
        assert (
            GaussianIntegrationToMesh._gaussianToMesh
            is not GaussianToMesh._gaussianToMesh
        )

    def test_the_base_class_accepts_both_forms_of_quantity(self):
        """Establishes the contract the subclass is expected to honour."""
        mesher = GaussianToMesh(dxdy=1.0)
        plain = mesher.gaussianToMesh(one_puff(quantity=1.0))
        with_units = mesher.gaussianToMesh(one_puff(quantity=1.0 * ureg.kg))
        assert float(plain.sum()) == pytest.approx(float(with_units.sum()))

    @pytest.mark.xfail(
        strict=True,
        reason="B29: the subclass narrows its parent's input. The base kernel "
               "multiplies by the plain number and converts a unit factor; the "
               "integration kernel calls unumToPint on the value itself, so a "
               "plain float becomes dimensionless and .m_as(kg) raises "
               "DimensionalityError. Code written against the base class breaks "
               "on the substitution its own docstring invites. "
               "See the consolidated findings issue.",
    )
    def test_the_subclass_accepts_a_plain_quantity_like_its_parent(self):
        mesher = GaussianIntegrationToMesh(dxdy=1.0)
        field = mesher.gaussianToMesh(one_puff(quantity=1.0))
        assert float(field.sum()) == pytest.approx(1.0, rel=1e-3)


@pytest.mark.unit
class TestSigmaColumnNames:
    def test_the_column_names_are_configurable(self, mesher):
        """Both have setters, so a caller can point them at their own columns."""
        mesher.sigmaXName = "myX"
        mesher.sigmaYName = "myY"
        assert mesher.sigmaXName == "myX"
        assert mesher.sigmaYName == "myY"

    def test_anisotropy_works_once_the_names_are_distinct(self, mesher):
        """Proves the kernel itself is fine -- only the default name is wrong.

        With sigmaYName pointed at its own column, the recovered crosswind
        sigma follows that column, so the machinery supports an elliptical
        plume.  Only the constructor's default prevents it.
        """
        mesher.sigmaYName = "sigmaYCorrected"
        field = mesher.gaussianToMesh(one_puff(sigmaX=5.0, sigmaY=15.0))

        alongwind = recovered_std(field.coords["x"].values, field.sum(dim="y").values)
        crosswind = recovered_std(field.coords["y"].values, field.sum(dim="x").values)

        assert alongwind == pytest.approx(5.0, rel=1e-2)
        assert crosswind == pytest.approx(15.0, rel=1e-2)

    @pytest.mark.xfail(
        strict=True,
        reason="B27: GaussianToMesh.__init__ sets sigmaYName to 'sigmaXCorrected', "
               "the X column, so the crosswind spread always reads the alongwind "
               "sigma and the plume is forced isotropic. Verified: an input "
               "sigmaYCorrected of 50 yields a recovered sigma_y of 5.000. "
               "See the consolidated findings issue.",
    )
    def test_the_default_crosswind_column_is_the_crosswind_one(self, mesher):
        assert mesher.sigmaYName == "sigmaYCorrected"

    @pytest.mark.xfail(
        strict=True,
        reason="B27: because sigmaYName defaults to the X column, a puff carrying "
               "a distinct sigmaYCorrected is silently spread isotropically. "
               "See the consolidated findings issue.",
    )
    def test_a_distinct_crosswind_sigma_is_honoured_by_default(self, mesher):
        """An atmospheric plume is elliptical; sigma_x and sigma_y differ by design."""
        field = mesher.gaussianToMesh(one_puff(sigmaX=5.0, sigmaY=50.0))
        crosswind = recovered_std(field.coords["y"].values, field.sum(dim="x").values)
        assert crosswind == pytest.approx(50.0, rel=1e-2)


@pytest.mark.unit
class TestDataFrameIndexing:
    def test_a_default_index_works(self, mesher):
        frame = pd.DataFrame([one_puff(x=0.0), one_puff(x=10.0)]).reset_index(drop=True)
        coordinates, _ = mesher._defineCoordinates(frame)
        assert len(coordinates) > 0

    @pytest.mark.xfail(
        strict=True,
        reason="B28: _defineCoordinates takes a LABEL from .index[0] and feeds it to "
               ".iloc, which wants a POSITION. Any frame whose index is not a "
               "0-based RangeIndex raises IndexError -- which covers anything that "
               "came out of a groupby, a filter or a concat. "
               "See the consolidated findings issue.",
    )
    def test_a_non_zero_based_index_works(self, mesher):
        """A filtered or grouped frame keeps its original labels."""
        frame = pd.DataFrame(
            [one_puff(x=0.0), one_puff(x=10.0)], index=[7, 8]
        )
        coordinates, _ = mesher._defineCoordinates(frame)
        assert len(coordinates) > 0

    @pytest.mark.xfail(
        strict=True,
        reason="B28: the same label-versus-position confusion, reached through the "
               "public entry point after a filter. "
               "See the consolidated findings issue.",
    )
    def test_a_filtered_frame_can_be_meshed(self, mesher):
        frame = pd.DataFrame(
            [one_puff(x=0.0), one_puff(x=10.0), one_puff(x=20.0)]
        ).reset_index(drop=True)
        filtered = frame[frame["x"] > 5.0]
        assert mesher.gaussianToMesh(filtered) is not None
