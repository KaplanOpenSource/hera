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
        super().__init__(projectName)


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
        with open(newTemplatePath, 'w') as fp:
            json.dump(newTemplate, fp)

        if loadToDB:
            self.logger.info("Saving template to the DB")
            self.addSimulationsDocument(resource=newTemplate['CaseDirectory'],
                                       dataFormat='string',
                                       type=self.templateDocType,
                                       desc=dict(OF_Workflow=newTemplate)) #desc=dict(OF_Workflow=newTemplate
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
    projectName = args.args[-1] if not args.noDB and funcName=='expand' else None

    handler = argsHandler(projectName)
    function = getattr(handler,f"{funcName}_handler")
    function(args)
