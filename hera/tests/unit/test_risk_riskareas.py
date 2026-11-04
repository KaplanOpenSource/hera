"""riskAreas: risk-area estimator dispatch and the sweep algorithm's settings.

``riskAreaAlgorithm_Sweep.calculate()`` drives a full geopandas/shapely
pipeline over a demographic layer and effect isopleths, forked across
``multiprocessing`` workers -- out of scope for a hermetic unit test. What is
covered here is the pure part: constructing the estimator, its property
round-trips, and the ``pydoc.locate`` dispatch in ``getRiskAreaAlgorithm``.

B60, resolved independently: this pass found ``getRiskAreaAlgorithm``
calling ``pydoc.locate("pyriskassessment.datalayer.riskAreas....")``, a
package that does not exist in this codebase or its dependencies (this
module actually lives at ``hera.riskassessment.analysis.riskAreas``), so
dispatch always raised -- even for ``"Sweep"``, the one algorithm actually
implemented here. That was a real, accurate finding at the time; commit
193249e6 ("fix: repair stale module path and nonexistent method call in
risk analysis") fixed the lookup to the real module path, unrelated to
this test-expansion effort. Dispatch by name now works correctly.
"""
import multiprocessing

import pytest

from hera.riskassessment.analysis.riskAreas import getRiskAreaAlgorithm, riskAreaAlgorithm_Sweep


@pytest.mark.unit
class TestConstruction:
    def test_the_documented_defaults_are_applied(self):
        sweep = riskAreaAlgorithm_Sweep()
        assert sweep.dxdy == pytest.approx(150.0)
        assert sweep.outlayers == 3
        assert sweep.parallel is True

    def test_worker_count_defaults_to_the_cpu_count(self):
        assert riskAreaAlgorithm_Sweep().workerCount == multiprocessing.cpu_count()

    def test_constructor_arguments_override_the_defaults(self):
        sweep = riskAreaAlgorithm_Sweep(dxdy=50, outlayers=1, parallel=False)
        assert (sweep.dxdy, sweep.outlayers, sweep.parallel) == (50.0, 1, False)


@pytest.mark.unit
class TestPropertySetters:
    def test_dxdy_setter_coerces_to_float(self):
        sweep = riskAreaAlgorithm_Sweep()
        sweep.dxdy = "75"
        assert sweep.dxdy == pytest.approx(75.0)
        assert isinstance(sweep.dxdy, float)

    def test_worker_count_setter_coerces_to_int(self):
        sweep = riskAreaAlgorithm_Sweep()
        sweep.workerCount = "4"
        assert sweep.workerCount == 4
        assert isinstance(sweep.workerCount, int)

    def test_parallel_setter_coerces_to_bool(self):
        sweep = riskAreaAlgorithm_Sweep()
        sweep.parallel = 0
        assert sweep.parallel is False
        sweep.parallel = 1
        assert sweep.parallel is True

    def test_outlayers_has_no_setter(self):
        """Unlike dxdy/parallel/workerCount, outlayers is read-only."""
        with pytest.raises(AttributeError):
            riskAreaAlgorithm_Sweep().outlayers = 5


@pytest.mark.unit
class TestGetRiskAreaAlgorithm:
    def test_an_unknown_algorithm_name_raises(self):
        with pytest.raises(ValueError, match="not found"):
            getRiskAreaAlgorithm("NoSuchAlgorithm")

    def test_the_error_message_lists_the_available_estimator(self):
        """B60, resolved: the real module now resolves, so the error
        message correctly lists 'Sweep' as available."""
        with pytest.raises(ValueError, match=r"Available estimators are: Sweep"):
            getRiskAreaAlgorithm("NoSuchAlgorithm")

    def test_sweep_is_the_only_implemented_algorithm(self):
        sweep = getRiskAreaAlgorithm("Sweep")
        assert isinstance(sweep, riskAreaAlgorithm_Sweep)

    def test_kwargs_are_forwarded_to_the_constructor(self):
        sweep = getRiskAreaAlgorithm("Sweep", dxdy=42, parallel=False)
        assert sweep.dxdy == pytest.approx(42.0)
        assert sweep.parallel is False
