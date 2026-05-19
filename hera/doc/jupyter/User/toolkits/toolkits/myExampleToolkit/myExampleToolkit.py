
class myExampleToolkit:
    def __init__(self, projectName=None, filesDirectory=None, **kwargs):
        self.projectName = projectName
        self.filesDirectory = filesDirectory

    def hello(self):
        return f"Hello from {self.__class__.__name__}"
