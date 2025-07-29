from hera.utils.logging import with_logger, get_classMethod_logger
from hera.toolkit import abstractToolkit

class hermesWorkflowToolkit(abstractToolkit):
    """
        The class handles machine/deep learning models.

        It helps saving hyper parameters and it provide simple
        tools (like batch/train splitting).

        Notes:
            * Torch models it requires pytorch installed.
            * SkiLearn requires scikitlearn installed.

    """

    def __init__(self, projectName: str, filesDirectory: str = None):
        """
            Initializes the machineLearning/deepLearning toolkit.

        Parameters
        ----------
        projectName: str
            The project where the models are stored.

        filesDirectory : str
            The directory to write all the Workflow and the outputs. default is current directory.
        """
        super().__init__(projectName=projectName,
                         filesDirectory=filesDirectory,
                         toolkitName= "machineLearningDeepLearningToolkit")

