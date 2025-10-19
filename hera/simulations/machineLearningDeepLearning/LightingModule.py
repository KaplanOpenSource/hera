
class LightningModuleHera:
    """
        This class adds the features of saving the container.
        THis will allow the user to reload the datasets every x epochs.
    """
    modelContainer = None

    def setModelContainer(self,modelConteiner):
        self.modelContainer = modelConteiner

