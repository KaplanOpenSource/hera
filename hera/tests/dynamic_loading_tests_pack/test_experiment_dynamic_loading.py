import os
import sys
import pytest


@pytest.mark.integration
def test_experiment_home_sees_dummy_experiment(load_dummy_experiment_to_project, hera_repo_root):
    """
    After we load the dummy experiment into the project,
    verify experimentHome sees it and can instantiate it.
    """
    proj = load_dummy_experiment_to_project["project"]
    exp_name = load_dummy_experiment_to_project["experiment"]

    sys.path.insert(0, str(hera_repo_root))

    from hera.measurements.experiment.experiment import experimentHome

    home = experimentHome(projectName=proj)
    keys = list(home.getExperimentsMap().keys())
    assert exp_name in keys

    exp = home.getExperiment(exp_name)
    assert exp is not None
    assert exp.__class__.__name__ == exp_name
