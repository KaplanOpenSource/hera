"""LSM/singleSimulation.py: SingleSimulation wraps an LSM run's final
xarray dataset and derives dosage/concentration from it.

B98: ``getConcentrationAtPoint`` -- called exactly as documented, with a
single scalar x/y/datetime (the "point" the name promises) -- always
raises ``IndexError``. ``.interp(x=x, y=y, datetime=datetime)`` on a
fully-specified single point collapses every matching dimension and
returns a 0-dimensional DataArray; ``.values`` on that is a 0-d numpy
array (shape ``()``), and indexing it with ``[0]`` raises
``IndexError: too many indices for array: array is 0-dimensional, but 1
were indexed``.
"""
import numpy
import pandas
import pytest
import xarray

from hera.simulations.LSM.singleSimulation import SingleSimulation
from hera.utils import ureg


def _dataset(datetime_coord):
    x = numpy.array([0.0, 1.0])
    y = numpy.array([0.0, 1.0])
    dosage = xarray.DataArray(
        numpy.arange(len(datetime_coord) * 2 * 2, dtype=float).reshape(len(datetime_coord), 2, 2),
        dims=["datetime", "x", "y"],
        coords={"datetime": datetime_coord, "x": x, "y": y},
    )
    return xarray.Dataset({"Dosage": dosage})


@pytest.fixture()
def datetime_times():
    return pandas.date_range("2020-01-01", periods=4, freq="1min")


@pytest.fixture()
def sim(datetime_times):
    return SingleSimulation(_dataset(datetime_times))


@pytest.mark.unit
class TestConstructionFromXarray:
    def test_a_dataset_is_used_directly(self, sim):
        assert isinstance(sim._finalxarray, xarray.Dataset)

    def test_a_document_delegates_to_get_data(self):
        ds = _dataset(pandas.date_range("2020-01-01", periods=2, freq="1min"))

        class _FakeDoc(dict):
            def getData(self):
                return ds

        doc = _FakeDoc({"desc": {"params": {"a": 1}, "version": [1, 0, 0]}})
        wrapped = SingleSimulation(doc)
        assert wrapped._finalxarray is ds

    def test_params_and_version_read_from_the_document_desc(self):
        class _FakeDoc(dict):
            def getData(self):
                return _dataset(pandas.date_range("2020-01-01", periods=2, freq="1min"))

        doc = _FakeDoc({"desc": {"params": {"a": 1}, "version": [1, 0, 0]}})
        wrapped = SingleSimulation(doc)
        assert wrapped.params == {"a": 1}
        assert wrapped.version == [1, 0, 0]


@pytest.mark.unit
class TestGetDosage:
    def test_dt_is_derived_from_the_datetime_spacing(self, sim):
        dosage = sim.getDosage()
        assert dosage.attrs["dt"].m_as(ureg.min) == pytest.approx(1.0)

    def test_it_also_works_with_a_numeric_datetime_coordinate(self):
        """The is_numeric_dtype branch: datetime as plain seconds, not datetime64."""
        numeric_times = numpy.array([0.0, 60.0, 120.0, 180.0])
        sim = SingleSimulation(_dataset(numeric_times))
        dosage = sim.getDosage()
        assert dosage.attrs["dt"].m_as(ureg.min) == pytest.approx(1.0)

    def test_dosage_is_scaled_by_the_q_factor(self, sim):
        default_q = sim.getDosage()
        double_q = sim.getDosage(Q=2 * ureg.kg)
        assert double_q["Dosage"].values[1, 0, 0] == pytest.approx(
            2 * default_q["Dosage"].values[1, 0, 0]
        )

    def test_q_and_c_attrs_use_the_requested_units(self, sim):
        dosage = sim.getDosage(q_units=ureg.g)
        assert dosage.attrs["Q"].units == ureg.g
        assert dosage.attrs["C"].units == ureg.g / ureg.m**3


@pytest.mark.unit
class TestGetConcentration:
    def test_it_has_one_fewer_datetime_step_than_dosage(self, sim):
        concentration = sim.getConcentration()
        assert len(concentration.datetime) == len(sim._finalxarray.datetime) - 1

    def test_c_is_the_dosage_derivative_over_dt(self, sim):
        concentration = sim.getConcentration()
        expected = concentration["dDosage"] / concentration.attrs["dt"].m_as(ureg.min)
        assert numpy.allclose(concentration["C"].values, expected.values)


@pytest.mark.unit
class TestGetConcentrationAtPointIsBroken:
    @pytest.mark.xfail(
        strict=True,
        reason="B98: interp(x=x, y=y, datetime=datetime) on a single fully "
               "-specified point returns a 0-dimensional DataArray; "
               ".values[0] on a 0-d array raises IndexError. "
               "See the consolidated findings issue.",
    )
    def test_it_should_return_the_concentration_at_a_point(self, sim, datetime_times):
        sim.getConcentrationAtPoint(x=0.5, y=0.5, datetime=datetime_times[1])

    def test_it_currently_raises_indexerror(self, sim, datetime_times):
        """Characterisation of B98."""
        with pytest.raises(IndexError, match="0-dimensional"):
            sim.getConcentrationAtPoint(x=0.5, y=0.5, datetime=datetime_times[1])
