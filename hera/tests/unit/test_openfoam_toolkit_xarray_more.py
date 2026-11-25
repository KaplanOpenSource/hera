"""OFToolkit.xarrayToSetFieldsDictDomain -- the last uncovered method of
openFoam/toolkit.py.

The only test that reached it lives in ``test_openfoam_dataset2of.py``, whose
module-level ``pytest.importorskip("foamlib")`` skips the whole file in this
venv, so the method was uncovered from line 781 down.  It needs no foamlib at
all: it builds ``boxToCell`` text out of an xarray dataset and takes nothing
from ``self`` but the logger, so a bare stand-in object is enough -- the
toolkit itself would need a project and a database.

What the method does, and what is asserted here: for every *cell* of the
dataset (one fewer than the number of coordinate values along each active
axis, because the coordinates are read as cell boundaries) it emits one
``boxToCell`` region whose ``box`` runs from coordinate ``i`` to coordinate
``i+1``, carrying one ``volScalarFieldValue`` / ``volVectorFieldValue`` /
``volTensorFieldValue`` entry per mapped field.  Values are read with
``isel`` at the *lower* index of the cell.  The assertions are derived from
that definition: the datasets below use coordinates 0, 1, 2 and a field
equal to a known function of the indices, so both the box corners and the
values are known in advance.

Passing ``None`` for an axis drops it from the product and from the box, so
the same method serves a 2-D dataset -- that is asserted too.

No bugs found: the four defects the docstring of the sibling file mentions
(``product`` handed integers, a string field falling into the vector branch,
the vector branch overwriting the dataset, the vector formatter joining
floats as strings) are all fixed in the current source.
"""
import numpy
import pytest
import xarray

from hera.simulations.openFoam.toolkit import OFToolkit

_EDGES = numpy.array([0.0, 1.0, 2.0])


class _LoggerHost:
    """Only `self` is used, and only to name the logger."""


@pytest.fixture()
def host():
    return _LoggerHost()


@pytest.fixture()
def dataset():
    """t[i, j, k] = 100i + 10j + k, so every isel value is readable by eye."""
    iGrid, jGrid, kGrid = numpy.meshgrid(numpy.arange(3.0), numpy.arange(3.0),
                                         numpy.arange(3.0), indexing="ij")
    return xarray.Dataset(
        data_vars=dict(t=(("x", "y", "z"), 100.0 * iGrid + 10.0 * jGrid + kGrid),
                       u=(("x", "y", "z"), numpy.ones((3, 3, 3))),
                       v=(("x", "y", "z"), numpy.full((3, 3, 3), 2.0))),
        coords=dict(x=_EDGES, y=_EDGES, z=_EDGES),
    )


def _boxLines(text):
    return [line.strip() for line in text.split("\n") if line.strip().startswith("box (")]


def _convert(host, dataset, **kwargs):
    arguments = dict(xarrayData=dataset, xColumnName="x", yColumnName="y", zColumnName="z")
    arguments.update(kwargs)
    return OFToolkit.xarrayToSetFieldsDictDomain(host, **arguments)


@pytest.mark.unit
class TestTheRegionsItEmits:
    def test_there_is_one_region_per_cell_not_per_coordinate_value(self, host, dataset):
        # 3 coordinate values along each of 3 axes -> 2 cells each -> 8 cells
        assert _convert(host, dataset, T="t").count("boxToCell") == 8

    def test_every_region_carries_a_box_and_a_fieldvalues_list(self, host, dataset):
        text = _convert(host, dataset, T="t")
        assert text.count("box (") == 8
        assert text.count("fieldValues") == 8

    def test_a_box_runs_from_one_coordinate_value_to_the_next(self, host, dataset):
        boxes = _boxLines(_convert(host, dataset, T="t"))
        assert "box (0.0 0.0 0.0 ) (1.0 1.0 1.0 );" in boxes
        assert "box (1.0 1.0 1.0 ) (2.0 2.0 2.0 );" in boxes

    def test_the_boxes_tile_the_domain_without_repeating_one(self, host, dataset):
        boxes = _boxLines(_convert(host, dataset, T="t"))
        assert len(set(boxes)) == 8

    def test_the_last_coordinate_value_is_only_ever_an_upper_corner(self, host, dataset):
        """The coordinates are cell boundaries, so the top one starts no box."""
        lowerCorners = [box.split(") (")[0] for box in _boxLines(_convert(host, dataset, T="t"))]
        assert all("2.0" not in corner for corner in lowerCorners)


@pytest.mark.unit
class TestTheFieldValues:
    def test_a_string_mapping_is_written_as_a_scalar(self, host, dataset):
        text = _convert(host, dataset, T="t")
        # the cell at indices (0, 0, 0) holds 100*0 + 10*0 + 0
        assert "volScalarFieldValue T 0.0" in text

    def test_the_value_is_read_at_the_lower_index_of_the_cell(self, host, dataset):
        text = _convert(host, dataset, T="t")
        # indices (1, 0, 1) -> 100 + 0 + 1
        assert "volScalarFieldValue T 101.0" in text
        # the top index 2 is never read, so 200/20/2 never appear as a value
        assert "volScalarFieldValue T 222.0" not in text

    def test_a_three_component_mapping_is_written_as_a_vector(self, host, dataset):
        assert "volVectorFieldValue U (1.0 2.0 0 )" in _convert(host, dataset, U=("u", "v", 0))

    def test_a_fixed_number_in_a_vector_is_written_through_unchanged(self, host, dataset):
        text = _convert(host, dataset, U=("u", "v", 7))
        assert "(1.0 2.0 7 )" in text

    def test_a_nine_component_mapping_is_written_as_a_tensor(self, host, dataset):
        nine = dataset.assign(**{f"c{index}": dataset["u"] for index in range(9)})
        text = _convert(host, nine, R=tuple(f"c{index}" for index in range(9)))
        assert "volTensorFieldValue R" in text
        assert "(1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 )" in text

    def test_several_fields_are_all_written_into_the_same_region(self, host, dataset):
        firstRegion = _convert(host, dataset, T="t", U=("u", "v", 0)).split("boxToCell")[1]
        assert "volScalarFieldValue T" in firstRegion
        assert "volVectorFieldValue U" in firstRegion

    def test_a_component_that_is_neither_a_name_nor_a_number_is_rejected(self, host, dataset):
        with pytest.raises(ValueError, match="not a string or number"):
            _convert(host, dataset, U=("u", "v", object()))

    @pytest.mark.parametrize("componentCount", [2, 4])
    def test_only_one_three_or_nine_components_make_a_field(self, host, dataset,
                                                            componentCount):
        mapping = tuple(["u"] * componentCount)
        with pytest.raises(ValueError, match="must be 1,3 or 9"):
            _convert(host, dataset, U=mapping)


@pytest.mark.unit
class TestTheAxesCanBeDroppedIndividually:
    def test_a_none_z_gives_two_dimensional_boxes(self, host, dataset):
        text = _convert(host, dataset.isel(z=0), zColumnName=None, T="t")
        assert text.count("boxToCell") == 4
        assert "box (0.0 0.0 ) (1.0 1.0 );" in _boxLines(text)

    def test_a_none_x_drops_the_first_component_of_the_box(self, host, dataset):
        text = _convert(host, dataset.isel(x=0), xColumnName=None, T="t")
        assert text.count("boxToCell") == 4
        # the corners are now (y, z) pairs, taken from the two remaining axes
        assert "box (0.0 0.0 ) (1.0 1.0 );" in _boxLines(text)
        # and the value read is t[0, j, k], i.e. below 100
        assert "volScalarFieldValue T 11.0" in text

    def test_a_single_axis_dataset_gives_one_region_per_interval(self, host, dataset):
        text = _convert(host, dataset.isel(y=0, z=0), yColumnName=None, zColumnName=None,
                        T="t")
        assert text.count("boxToCell") == 2


@pytest.mark.unit
class TestTheTimeSelection:
    @pytest.fixture()
    def timeDataset(self, dataset):
        withTime = dataset.expand_dims(time=[0.0, 10.0]).copy()
        withTime["t"] = withTime["t"] + xarray.DataArray([0.0, 1000.0], dims="time")
        return withTime

    def test_the_nearest_time_step_is_selected(self, host, timeDataset):
        text = _convert(host, timeDataset, time=9.0, T="t")
        assert "volScalarFieldValue T 1000.0" in text
        assert "volScalarFieldValue T 0.0" not in text

    def test_the_earlier_step_is_chosen_when_it_is_the_nearer_one(self, host, timeDataset):
        text = _convert(host, timeDataset, time=1.0, T="t")
        assert "volScalarFieldValue T 0.0" in text
        assert "volScalarFieldValue T 1000.0" not in text

    def test_the_name_of_the_time_coordinate_is_an_argument(self, host, dataset):
        renamed = dataset.expand_dims(Times=[0.0, 10.0]).copy()
        renamed["t"] = renamed["t"] + xarray.DataArray([0.0, 1000.0], dims="Times")
        text = _convert(host, renamed, time=9.0, timeColumn="Times", T="t")
        assert "volScalarFieldValue T 1000.0" in text

    def test_without_a_time_the_index_is_ignored_and_every_step_is_kept(self, host,
                                                                       timeDataset):
        """isel over the spatial axes only, so the time dimension survives and
        .item() would fail on more than one element -- which is what makes
        `time=None` on a dataset that still has a time axis unusable."""
        with pytest.raises(ValueError):
            _convert(host, timeDataset, T="t")


@pytest.mark.unit
class TestTheReturnValue:
    def test_it_is_a_single_string_of_the_regions_only(self, host, dataset):
        text = _convert(host, dataset, T="t")
        assert isinstance(text, str)
        # the regions text only: no FoamFile header, no defaultFieldValues
        assert "FoamFile" not in text
        assert "defaultFieldValues" not in text

    def test_an_empty_field_map_still_emits_one_region_per_cell(self, host, dataset):
        text = _convert(host, dataset)
        assert text.count("boxToCell") == 8
        assert "FieldValue" not in text
