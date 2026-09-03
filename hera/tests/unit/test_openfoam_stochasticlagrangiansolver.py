"""StochasticLagrangianSolver_toolkitExtension.createDispersionFlowField --
the only member of that module left uncovered.

The module is four lines of real code: a constructor that forwards to the
abstract Lagrangian extension, and an override of
``createDispersionFlowField`` whose entire job is to guarantee that the two
fields the solver cannot run without -- ``ustar`` and
``distanceFromWalls`` -- are present in ``flowData['dispersionFields']``,
defaulting to 0, before delegating to the base implementation.

The base implementation copies an OpenFOAM case directory tree, links or
copies its mesh, writes dispersion fields with PyFoam and registers a
document in MongoDB; it is covered elsewhere.  Here it is monkeypatched on
the *base class* (never on an instance) and used as a seam: what the
override put into ``flowData``, which keyword arguments crossed, and what
came back.

The extension is built with ``__new__`` rather than its constructor, because
``absractStochasticLagrangianSolver_toolkitExtension.__init__`` opens a
``dask.distributed.Client()``, which the unit layer's socket guard refuses
(rightly -- it would start a scheduler process).  Nothing the override does
reads any attribute the constructor sets, so ``__new__`` is sufficient and
keeps the test hermetic.

Bug pinned here (with a strict xfail for the intended behaviour and a
passing characterisation of what happens today):

* B293 the override calls ``super().createDispersionFlowField(...)``
  without ``return``, so the document the base implementation registers and
  returns -- its own docstring's "Returns" is the dispersion-flow document --
  is thrown away and every caller of the concrete solver gets ``None``.  The
  base method's last statement is ``return ret``; every other override in
  this hierarchy that wraps a returning method returns it.
"""
import pytest

from hera.simulations.openFoam import FLOWTYPE_COMPRESSIBLE, FLOWTYPE_INCOMPRESSIBLE
from hera.simulations.openFoam.lagrangian.abstractLagrangianSolver import (
    absractStochasticLagrangianSolver_toolkitExtension,
)
from hera.simulations.openFoam.lagrangian.StochasticLagrangianSolver import (
    StochasticLagrangianSolver_toolkitExtension,
)

_SENTINEL_DOCUMENT = object()


@pytest.fixture()
def extension():
    """The concrete extension without its dask-Client-opening constructor."""
    return StochasticLagrangianSolver_toolkitExtension.__new__(
        StochasticLagrangianSolver_toolkitExtension)


@pytest.fixture()
def baseSeam(monkeypatch):
    """Replace the base implementation on the CLASS; return the call list."""
    calls = []

    def recorder(self, **kwargs):
        calls.append(kwargs)
        return _SENTINEL_DOCUMENT

    monkeypatch.setattr(absractStochasticLagrangianSolver_toolkitExtension,
                        "createDispersionFlowField", recorder)
    return calls


def _flowData(**dispersionFields):
    return dict(originalFlow=dict(time=dict(temporalType="steadyState", timestep=400)),
                dispersionFields=dict(dispersionFields))


@pytest.mark.unit
class TestTheRequiredDispersionFields:
    def test_ustar_is_added_when_the_caller_did_not_ask_for_it(self, extension, baseSeam):
        flowData = _flowData(Hmix=1000)
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["dispersionFields"]["ustar"] == 0

    def test_distancefromwalls_is_added_too(self, extension, baseSeam):
        flowData = _flowData(Hmix=1000)
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["dispersionFields"]["distanceFromWalls"] == 0

    def test_a_value_the_caller_supplied_is_not_overwritten(self, extension, baseSeam):
        flowData = _flowData(ustar=0.25, distanceFromWalls=3.5)
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["dispersionFields"]["ustar"] == 0.25
        assert flowData["dispersionFields"]["distanceFromWalls"] == 3.5

    def test_a_falsy_value_the_caller_supplied_is_still_not_overwritten(self, extension,
                                                                       baseSeam):
        """setdefault, not `or`: an explicit 0.0 must survive."""
        flowData = _flowData(ustar=0.0)
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["dispersionFields"]["ustar"] == 0.0

    def test_the_other_dispersion_fields_are_left_alone(self, extension, baseSeam):
        flowData = _flowData(Hmix=1000)
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["dispersionFields"]["Hmix"] == 1000

    def test_the_caller_s_own_dictionary_is_the_one_that_is_updated(self, extension, baseSeam):
        """The defaults must be visible to the caller, not only to super()."""
        flowData = _flowData()
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert baseSeam[0]["flowData"] is flowData

    def test_the_rest_of_the_flow_data_is_untouched(self, extension, baseSeam):
        flowData = _flowData()
        original = dict(flowData["originalFlow"])
        extension.createDispersionFlowField(flowName="f", flowData=flowData,
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert flowData["originalFlow"] == original

    def test_a_flow_data_without_a_dispersionfields_key_is_not_invented(self, extension,
                                                                       baseSeam):
        """Characterisation: the override indexes ['dispersionFields'] directly."""
        with pytest.raises(KeyError, match="dispersionFields"):
            extension.createDispersionFlowField(flowName="f", flowData=dict(),
                                                OriginalFlowField="orig",
                                                dispersionDuration=100)
        assert baseSeam == []


@pytest.mark.unit
class TestWhatIsForwardedToTheBaseImplementation:
    def test_every_argument_crosses_by_keyword(self, extension, baseSeam):
        extension.createDispersionFlowField(flowName="myFlow", flowData=_flowData(),
                                            OriginalFlowField="origFlow",
                                            dispersionDuration=3600)
        assert baseSeam[0]["flowName"] == "myFlow"
        assert baseSeam[0]["OriginalFlowField"] == "origFlow"
        assert baseSeam[0]["dispersionDuration"] == 3600

    def test_the_flow_type_defaults_to_incompressible(self, extension, baseSeam):
        extension.createDispersionFlowField(flowName="f", flowData=_flowData(),
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert baseSeam[0]["flowType"] == FLOWTYPE_INCOMPRESSIBLE

    def test_an_explicit_flow_type_is_forwarded(self, extension, baseSeam):
        extension.createDispersionFlowField(flowName="f", flowData=_flowData(),
                                            OriginalFlowField="orig", dispersionDuration=100,
                                            flowType=FLOWTYPE_COMPRESSIBLE)
        assert baseSeam[0]["flowType"] == FLOWTYPE_COMPRESSIBLE

    def test_overwrite_defaults_to_false_and_db_support_to_true(self, extension, baseSeam):
        extension.createDispersionFlowField(flowName="f", flowData=_flowData(),
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert baseSeam[0]["overwrite"] is False
        assert baseSeam[0]["useDBSupport"] is True

    @pytest.mark.parametrize("overwrite, useDBSupport", [(True, False), (False, True)])
    def test_both_flags_are_forwarded_as_given(self, extension, baseSeam, overwrite,
                                               useDBSupport):
        extension.createDispersionFlowField(flowName="f", flowData=_flowData(),
                                            OriginalFlowField="orig", dispersionDuration=100,
                                            overwrite=overwrite, useDBSupport=useDBSupport)
        assert baseSeam[0]["overwrite"] is overwrite
        assert baseSeam[0]["useDBSupport"] is useDBSupport

    def test_the_base_implementation_is_called_exactly_once(self, extension, baseSeam):
        extension.createDispersionFlowField(flowName="f", flowData=_flowData(),
                                            OriginalFlowField="orig", dispersionDuration=100)
        assert len(baseSeam) == 1


@pytest.mark.unit
class TestTheReturnValue:
    @pytest.mark.xfail(
        strict=True,
        reason="B293: the override ends with a bare "
               "super().createDispersionFlowField(...) call and no `return`, "
               "so the dispersion-flow document the base implementation "
               "registers and returns (its last statement is `return ret`) is "
               "discarded and every caller of the concrete stochastic "
               "Lagrangian solver receives None.  See the consolidated "
               "findings issue.",
    )
    def test_the_document_the_base_returns_should_reach_the_caller(self, extension, baseSeam):
        assert extension.createDispersionFlowField(
            flowName="f", flowData=_flowData(), OriginalFlowField="orig",
            dispersionDuration=100) is _SENTINEL_DOCUMENT

    def test_it_currently_returns_none_however_the_base_answered(self, extension, baseSeam):
        """Characterisation of B293."""
        assert extension.createDispersionFlowField(
            flowName="f", flowData=_flowData(), OriginalFlowField="orig",
            dispersionDuration=100) is None
        # the base really was called and really did answer
        assert len(baseSeam) == 1

    def test_the_base_method_itself_does_return_a_value(self):
        """Characterisation of B293: there is something to lose."""
        import inspect

        source = inspect.getsource(
            absractStochasticLagrangianSolver_toolkitExtension.createDispersionFlowField)
        assert source.rstrip().endswith("return ret")
