from jinja2 import Environment, FileSystemLoader

class InputForModelsCreator(object):
    _paramsMap = None
    _templateName = None
    _templatesDir = None
    _renderedTemplate = None

    def __init__(self, templatesDir):
        self._templatesDir = templatesDir

    @property
    def paramsMap(self):
        return self._paramsMap

    @property
    def templateName(self):
        return self._templateName

    @property
    def renderedTemplate(self):
        return self._renderedTemplate

    def setParamsMap(self, paramsMap):
        self._paramsMap = paramsMap

    def setTemplate(self, templateName):
        self._templateName = templateName

    def render(self, savePath=None):
        if self._templateName is None or self._paramsMap is None:
            print("templateName and paramsMap are not set yet")
        else:
            file_loader = FileSystemLoader(self._templatesDir)
            env = Environment(loader=file_loader)
            template = env.get_template(self._templateName)
            renderedTemplate = template.render(self._paramsMap)
            if savePath is not None:
                with open(savePath, 'w') as file:
                    file.write(renderedTemplate)
            self._renderedTemplate = renderedTemplate
        return renderedTemplate