


import functools
import os
import logging as _stdlib_logging

_initialized = False
ToolkitHome = toolkitHome = loadJSON = None
_get_logger = None  # hera's get_logger, deferred

def _setup():
    global _initialized, ToolkitHome, toolkitHome, loadJSON, _get_logger
    if _initialized:
        return
    _initialized = True
    from hera import ToolkitHome as _TH, toolkitHome as _th
    ToolkitHome, toolkitHome = _TH, _th
    from hera.utils.jsonutils import loadJSON as _lj
    loadJSON = _lj
    from hera.utils.logging import get_logger as _gl
    _get_logger = _gl

def _lazy_setup(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        _setup()
        return func(*args, **kwargs)
    return wrapper

# Compatibility shim so existing 'logging.get_logger(...)' calls still work
class _LoggingShim:
    def get_logger(self, *a, **kw):
        if _get_logger is not None:
            return _get_logger(*a, **kw)
        return _stdlib_logging.getLogger(*a, **kw)
    def __getattr__(self, item):
        return getattr(_stdlib_logging, item)

logging = _LoggingShim()

@_lazy_setup
def _confirm_project_name(arguments, logger):
    if arguments.projectName is None:
        logger.debug(
            f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        arguments.projectName = caseConfiguration['projectName']

@_lazy_setup
def list_templates(arguments):
    # for template in os.listdir("templates"):
    logger = logging.get_logger("hera.bin.hera_lsm.load_template")
    _confirm_project_name(arguments, logger)
    lsm = toolkitHome.getToolkit(toolkitName=ToolkitHome.LSM, projectName=arguments.projectName)
    templates = lsm.getTemplates()
    if len(templates)==0:
        print("There are no templates")
        return
    print("---====[ List of template]====---")    
    for template in templates:
        print(f"* {template.templateName}:")
        print(f"\t version: {template.version}")
        print(f"\t path: {template.dirPath}")
        print(f"\t model folder: {template.modelFolder}")
    
@_lazy_setup
def setup_template(arguments):
    # for template in os.listdir("templates"):
    logger = logging.get_logger("hera.bin.hera_lsm.load_template")
    _confirm_project_name(arguments, logger)
    lsm = toolkitHome.getToolkit(toolkitName=ToolkitHome.LSM, projectName=arguments.projectName)
    template = lsm.getTemplateByName(arguments.templateName,templateVersion=arguments.templateVersion)
    if template is None:
        raise RuntimeError(f"The template {arguments.templateName} is not in the project. You must load the appropriate repository that contains the template and then do `hera-project project updateRepositories`")
    
    logger.info(f"Found template {arguments.templateName}")


    codeDir = os.path.join(arguments.codeDir, arguments.templateName.split("-")[0])
    outDir = template.modelFolder
    metDir = os.path.join(outDir,"tozaot","Meteorology")

    try:
        os.unlink(metDir)
        logger.info(f"unlinked {metDir} from template")
    except FileNotFoundError:
        logger.error(f"file {metDir} not found")
        pass

    try:
        os.unlink(os.path.join(outDir,"a.out"))
        logger.info(f"unlinked {os.path.join(outDir,'a.out')} from template")
    except FileNotFoundError:
        logger.error(f"file {os.path.join(outDir,'a.out')} not found")
        pass

    os.makedirs(os.path.join(outDir,'tozaot'),exist_ok=True)
    os.makedirs(os.path.join(outDir,'tozaot','machsan'),exist_ok=True)
    logger.info(f"setup directories successfully")


    os.system(f"ln -s {os.path.join(codeDir,'a.out')} {outDir}")
    logger.info(f"linked {os.path.join(outDir,'a.out')} from template")
    os.system(f"ln -s {os.path.join(codeDir,'tozaot/Meteorology')} {metDir}")
    logger.info(f"linked {os.path.join(codeDir,'tozaot/Meteorology')} from template")
