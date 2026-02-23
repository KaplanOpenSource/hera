
from hera.measurements.experiment.experiment import experimentSetupWithData
from hera.measurements.experiment.analysis import experimentAnalysis
from hera.measurements.experiment.presentation import experimentPresentation


class DummyExperiment(experimentSetupWithData):
    def __init__(self,projectName,pathToExperiment,filesDirectory):
        super().__init__(projectName=projectName,pathToExperiment=pathToExperiment,filesDirectory=filesDirectory)
        self._analysis = DummyExperimentAnalysis(self)
        self._presentation = DummyExperimentPresentation(self,self.analysis)

class DummyExperimentAnalysis(experimentSetupWithData):
    datalayer = None 
    def __init__(self,datalayer):
        self.datalayer = datalayer


class DummyExperimentPresentation(experimentSetupWithData):
    datalayer = None 
    analysis  = None 
    
    def __init__(self,datalayer,analysis):
        self.datalayer = datalayer
        self.analysis  = analysis

