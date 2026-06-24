


import os

from hera import ToolkitHome, toolkitHome
from hera.utils import logging
from hera.utils.jsonutils import loadJSON

def _confirm_project_name(arguments, logger):
    if arguments.projectName is None:
        logger.debug(
            f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        arguments.projectName = caseConfiguration['projectName']

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


    src_aout = os.path.join(codeDir, 'a.out')
    dst_aout = os.path.join(outDir, 'a.out')
    if not os.path.exists(dst_aout):
        os.symlink(src_aout, dst_aout)
    logger.info(f"linked {dst_aout} from template")
    src_met = os.path.join(codeDir, 'tozaot', 'Meteorology')
    if not os.path.exists(metDir):
        os.symlink(src_met, metDir)
    logger.info(f"linked {src_met} from template")
