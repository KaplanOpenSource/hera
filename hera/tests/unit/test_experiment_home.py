"""experimentHome: the factory/registry side of the experiment toolkit.

``experimentSetupWithData`` and everything under it need a real argos
experiment zip file to construct, so this file only covers the
``experimentHome`` factory itself, backed by mongomock like every other
toolkit.
"""
import pytest

from hera import toolkitHome


@pytest.fixture()
def home(unit_toolkit_factory):
    return unit_toolkit_factory(toolkitHome.EXPERIMENT)


@pytest.mark.unit
class TestExperimentMap:
    def test_a_fresh_project_has_no_experiments(self, home):
        assert home.experimentMap == {}

    def test_experiments_table_is_the_same_as_get_experiments_table(self, home):
        assert home.experimentsTable.equals(home.getExperimentsTable())

    def test_keys_lists_no_experiments_on_a_fresh_project(self, home):
        assert home.keys() == []


@pytest.mark.unit
class TestExperimentDataType:
    def test_it_defaults_to_none(self, home):
        assert home.experimentDataType() is None


@pytest.mark.unit
class TestGetExperiment:
    def test_an_unknown_experiment_name_raises(self, home):
        with pytest.raises(ValueError, match="load the experiment"):
            home.getExperiment("NoSuchExperiment")

    def test_dict_style_access_delegates_to_get_experiment(self, home):
        with pytest.raises(ValueError, match="load the experiment"):
            home["NoSuchExperiment"]
