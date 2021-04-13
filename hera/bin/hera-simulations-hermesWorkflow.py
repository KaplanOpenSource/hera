#! /usr/bin/env python

import argparse
from hermes import expandWorkflow
from hermes import hermesWorkflow
import json
import os
import pathlib
from hera.datalayer import Project
from pathlib import Path


class argsHandler(Project):

    templateDocType = "HermesOpenFOAM"

    def __init__(self,projectName=None):
        projectName = "OpenFoamRuns" if projectName is None else projectName
        super().__init__(projectName,loggerName="bin")
        self.logger.info("Initialize")

    def _expand_and_load(self,templatePath,newTemplatePath,loadToDB=True):

        """
        parameters
        ----------
        templatePath: string. the fileName/path to workflow json file
        newTemplatePath: string. the fileName/path for resulted expanded workflow json file
        loadToDB: boolean. load/not the workflow to DB. determined by the -noDB flag

        """

        expander = expandWorkflow()
        newTemplate = expander.expand(templatePath)

        if loadToDB:
            self.logger.info("Saving template to the DB")
            doc=self.addSimulationsDocument(resource=newTemplate['workflow']['nodes']['Parameters']['GUI']["WebGui"]["formData"]['caseDirectory'],
                                       dataFormat='string',
                                       type=self.templateDocType,
                                       desc=dict(OF_workflow=newTemplate)) #desc=dict(OF_workflow=newTemplate
            print(doc.id)
            newPath=os.path.join(newTemplate['workflow']['nodes']['Parameters']['GUI']["WebGui"]["formData"]['caseDirectory'],str(doc.id))
            #newPath=os.path.join(newTemplate['CaseDirectory'],str(doc.id))
            #newTemplate['workflow']['nodes']['copyDirectory']["GUI"]["Properties"]["Target"]["current_val"]=newPath
            newTemplate['workflow']['nodes']['Parameters']['GUI']["WebGui"]["formData"]['caseDirectory']=newPath

            if 'copyDirectory' in newTemplate['workflow']['nodeList']:
                sourceFolder = newTemplate['workflow']['nodes']['Parameters']['GUI']["WebGui"]["formData"]['OFtemplateDirectory']
                json_files = [f for f in os.listdir(sourceFolder) if f.endswith('.json')]

                if len(json_files) > 1:
                    raise self.logger.warning(f"found more then single json in {sourceFolder}. only parameter json file should exists")
                try:
                    with open(os.path.join(sourceFolder, json_files[0])) as f:
                        templateParams = json.load(f)

                    newTemplate['templateParams'] = templateParams
                except:
                    self.logger.warning(f"no template parameter json file in {sourceFolder}")
            self.logger.info("updating template to the DB")
            doc.resource=newPath
            doc.desc['OF_workflow']=newTemplate
            doc.save()
            self.logger.info("workflow parameters updated and saved to DB")


        with open(newTemplatePath, 'w') as fp:
            json.dump(newTemplate, fp)

        self.logger.info("Done")


    def _build(self,templatePath,WDPath,builder,pythonPath):

        """

        parameters
        ----------
        templatePath: string. the fileName/path to the expanded workflow json file
        WDPath:
        builder:
        pythonPath: string. the fileName/path for resulted python file

        """


        flow = hermesWorkflow(templatePath, WDPath,"")
        build = flow.build(builder)
        with open(pythonPath, "w") as file:
            file.write(build)

        self.logger.info("Done")


    def _executeLuigi(self,pythonPath):

        """

        parameters
        ----------
        pythonPath: string. the fileName/path of the python file

        """

        cwd = pathlib.Path().absolute()
        moduleParent = pathlib.Path(pythonPath).parent.absolute()
        os.chdir(moduleParent)
        os.system(f"python3 -m luigi --module {os.path.basename(pythonPath)} finalnode_xx_0 --local-scheduler")
        os.chdir(cwd)


    def expand_handler(self,args):
        """
        parameters
        ----------
        args: argparse object' resulted from CLI inputs

        """
        arguments=args.args
        templatePath = arguments[0]
        newTemplatePath = arguments[1]
        loadToDB=False if args.noDB else True

        self._expand_and_load(templatePath, newTemplatePath, loadToDB)


    def buildPython_handler(self,args):
        """
        parameters
        ----------
        args: argparse object' resulted from CLI inputs

        """

        arguments=args.args
        templatePath = arguments[0]
        pythonPath = arguments[1]
        WDPath = arguments[2] if len(arguments) > 2 else str(pathlib.Path(pythonPath).parent.absolute())
        builder = arguments[3] if len(arguments) > 3 else "luigi"

        self._build(templatePath,WDPath,builder,pythonPath)


    def executeLuigi_handler(self,args):
        """
        parameters
        ----------
        args: argparse object' resulted from CLI inputs

        """

        arguments=args.args
        pythonPath = arguments[0]

        self._executeLuigi(pythonPath)


    def runAll_handler(self,args):
        """
        parameters
        ----------
        args: argparse object' resulted from CLI inputs

        """

        arguments=args.args

        with open(arguments[0]) as f:
            argDict = json.load(f)

        templatePath=argDict["templatePath"]
        newTemplatePath=argDict["newTemplatePath"]
        loadToDB=False if args.noDB else True

        pythonPath=argDict.get('pythonPath')
        WDPath=argDict.get('WDPath',str(pathlib.Path(pythonPath).parent.absolute()))
        builder = argDict.get('builder', "luigi")

        self._expand_and_load(templatePath,newTemplatePath,loadToDB)
        self._build(newTemplatePath,WDPath,builder,pythonPath)
        self._executeLuigi(Path(pythonPath).stem)


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('command', nargs=1, type=str)
    parser.add_argument('args', nargs='*', type=str)
    parser.add_argument('-noDB', action='store_true')
    args = parser.parse_args()
    funcName = args.command[0]
    projectName = args.args[-1] if not args.noDB and funcName in ['runAll','expand'] else None

    handler = argsHandler(projectName)
    function = getattr(handler,f"{funcName}_handler")
    function(args)
