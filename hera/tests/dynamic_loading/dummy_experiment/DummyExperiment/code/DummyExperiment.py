from hera.measurements.experiment.experiment import experimentSetupWithData


class DummyExperiment(experimentSetupWithData):
    """
    Dummy experiment used only for testing dynamic experiment loading.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def hello(self):
        return "hello from DummyExperiment"
