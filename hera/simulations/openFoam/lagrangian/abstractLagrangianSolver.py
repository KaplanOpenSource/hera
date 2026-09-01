import io
import logging
import subprocess

import dask.dataframe
import numpy
import pandas
import os
import glob
import shutil
import json
import xarray
from itertools import product

from dask.delayed import delayed
import dask.dataframe as dask_dataframe
from pandas.core.array_algos.transforms import shift
from dask.distributed import Client

from hera.utils.unitHandler import *

from hera.datalayer import datatypes
from hera.datalayer.document.metadataDocument import MetadataFrame
from hera.utils import dictToMongoQuery,JSONToConfiguration,get_classMethod_logger,ConfigurationToJSON,tonumber
from hera.simulations.openFoam  import FLOWTYPE_INCOMPRESSIBLE
from hera.simulations.openFoam.OFWorkflow import workflow_StochasticLagrangianSolver

class absractStochasticLagrangianSolver_toolkitExtension:
    """
        A base class that handles all the stochastic lagrangian things (datalauer
    """
    toolkit = None

    analysis = None  # link to the analysis of the StochasticLagrnagian
    presentation = None  # link to the presentation of the stochasticLagrangian

    DOCTYPE_LAGRANGIAN_CACHE = "lagrangianCacheDocType"
    DOCTYPE_CONCENTRATIONEULERIAN_CACHE = "EulerianConcentrationCacheDocType"

    def __init__(self, toolkit):
        """Initialize with a reference to the parent toolkit."""
        self.toolkit = toolkit
        self.daskClient = Client()
        self.analysis = analysis(self)

    def createDispersionFlowField(self, flowName, flowData, OriginalFlowField, dispersionDuration,
                                  flowType=FLOWTYPE_INCOMPRESSIBLE, overwrite: bool = False, useDBSupport: bool = True):
        """
        Prepare the case directory of the flow for dispersion and assign a local name.

        Orchestrates five steps: locate original flow → detect parallelism and
        timesteps → build time mapping → check/register in DB → copy/link
        directories and write dispersion fields.

        Parameters
        ----------
        flowName : str
            Suffix appended to the group name for this dispersion flow.
        flowData : dict
            Configuration dict with ``originalFlow`` and ``dispersionFields``.
        OriginalFlowField : str
            Name or path of the base flow (DB name, directory, or resource).
        dispersionDuration : float
            Maximum timestep of the dispersion simulation.
        flowType : str
            ``"incompressible"`` or ``"compressible"``.
        overwrite : bool
            If True, replace existing dispersion flow.
        useDBSupport : bool
            If False, skip all DB queries and registration.

        Returns
        -------
        str or None
            Path to the dispersion flow directory, or None if no DB support.
        """
        logger = get_classMethod_logger(self, "createDispersionFlowField")
        logger.info(f"Creating dispersion flow field : {OriginalFlowField} ")
        originalFlow = dict(flowData['originalFlow'])
        originalFlow['source'] = OriginalFlowField
        dispersionDuration = self._toTimeFormat(dispersionDuration)

        # Step 1: Locate the original flow case directory.
        originalFlowCaseDir, workflowGroup = self._locateOriginalFlowCase(originalFlow)
        workflowGroup = f"{workflowGroup}_DFF"
        dispersionFlowFieldName = f"{workflowGroup}_{flowName}"
        logger.info(f"Found original flow: {originalFlowCaseDir}. Group: {workflowGroup}, Name: {dispersionFlowFieldName}")

        # Step 2: Detect parallel vs serial and discover timesteps.
        parallelOriginal, TS = self._detectParallelAndTimesteps(originalFlowCaseDir)

        # Step 3: Build the time mapping (original_time → dispersion_time).
        dynamicType = originalFlow['time']['temporalType']
        timeStep = originalFlow.get("timeStep", None)
        timeList, uts = self._buildTimeMapping(dynamicType, TS, timeStep, dispersionDuration)
        logger.info(f"Simulation type: {dynamicType}. Time mapping (orig,dest): {timeList}")

        linkMeshSymbolically = originalFlow.get("linkMeshSymbolically", True)
        linkDataSymbolically = originalFlow.get("linkDataSymbolically", True)

        # Step 4: Check DB and handle existing records.
        originalFlow['baseFlowDirectory'] = originalFlowCaseDir
        originalFlow['timeStep'] = uts
        doc = self._checkAndRegisterDispersionInDB(
            useDBSupport, dispersionFlowFieldName, workflowGroup,
            originalFlow, flowData, dispersionDuration, overwrite
        )

        # Handle overwrite of existing directory.
        if doc is not None:
            if overwrite:
                resource = doc.resource
                logger.info(f"Overwriting existing dispersion field at {resource}")
                if os.path.exists(resource):
                    shutil.rmtree(resource)
            else:
                raise FileExistsError("Dispersion flow field already exists. Use overwrite=True.")

        # Step 5: Create the dispersion case directory structure.
        dispersionFlowFieldDirectory = os.path.abspath(
            os.path.join(self.toolkit.FilesDirectory, dispersionFlowFieldName))
        self._prepareDispersionDirectory(
            dispersionFlowFieldDirectory, originalFlowCaseDir,
            parallelOriginal, timeList, linkMeshSymbolically,
            linkDataSymbolically, flowData, flowType
        )

        # Step 6: Register in DB.
        ret = self._finalizeDispersionInDB(
            useDBSupport, doc, dispersionFlowFieldName, workflowGroup,
            dispersionFlowFieldDirectory, originalFlow, flowData, dispersionDuration
        )

        logger.info("... Done")
        return ret

    # ------------------------------------------------------------------
    # Private helpers for createDispersionFlowField
    # ------------------------------------------------------------------

    @staticmethod
    def _toTimeFormat(timeDirectory):
        """Convert a time directory name to int if possible, else float."""
        int_version = int(float(timeDirectory))
        float_version = float(timeDirectory)
        return int_version if float_version == int_version else float_version

    def _locateOriginalFlowCase(self, originalFlow):
        """Find the original flow case directory from DB or filesystem.

        Returns
        -------
        tuple of (str, str)
            (caseDirectory, workflowGroup)
        """
        logger = get_classMethod_logger(self, "_locateOriginalFlowCase")
        source = originalFlow['source']
        logger.debug(f"Trying to find flow {source} in the DB")
        docList = self.toolkit.getWorkflowListDocumentFromDB(source)

        if len(docList) == 0:
            logger.debug(f"Not in DB. Trying as a directory.")
            if not os.path.isdir(source):
                raise FileNotFoundError(
                    f"The original flow {source} is not a directory and does not exist in the DB.")
            if not (os.path.isdir(os.path.join(source, 'system')) and
                    os.path.isdir(os.path.join(source, 'constant'))):
                raise ValueError(
                    f"The directory {source} is not a case directory (no system/ or constant/).")
            return source, os.path.basename(source)

        if len(docList) > 1:
            raise ValueError(f"The name {source} matches more than one simulation.")

        flowCaseDir= docList[0].resource
        workflowName = docList[0].desc['workflowName']
        if docList[0].type == "hermesWorkflow":
            flowCaseDir = os.path.join(os.path.dirname(flowCaseDir), workflowName)

        return flowCaseDir, workflowName

    def _detectParallelAndTimesteps(self, caseDir):
        """Detect parallel case structure and discover available timesteps.

        Returns
        -------
        tuple of (bool, list of float)
            (isParallel, sortedTimesteps)
        """
        logger = get_classMethod_logger(self, "_detectParallelAndTimesteps")
        if os.path.exists(os.path.join(caseDir, "processor0")):
            logger.debug("Found processor0 — parallel case")
            ptPath = ["processor0", "*"]
            isParallel = True
        else:
            logger.debug("No processor0 — single processor case")
            ptPath = ["*"]
            isParallel = False

        # Time directories have names that parse as numbers.
        TS = [float(os.path.basename(ts))
              for ts in glob.glob(os.path.join(caseDir, *ptPath))
              if os.path.basename(ts).replace(".", "").isdigit()]
        TS.sort()
        logger.info(f"Found timesteps: {TS}")
        return isParallel, TS

    def _buildTimeMapping(self, dynamicType, TS, timeStep, dispersionDuration):
        """Build the (original_time, dispersion_time) mapping list.

        Returns
        -------
        tuple of (list of tuple, float)
            (timeList, selectedTimestep)
        """
        logger = get_classMethod_logger(self, "_buildTimeMapping")
        fmt = self._toTimeFormat

        if dynamicType == self.toolkit.TIME_STEADYSTATE:
            # Steady state: freeze flow at one timestep.
            if timeStep is None:
                uts = TS[-1]
            else:
                uts = TS[min(range(len(TS)), key=lambda i: abs(TS[i] - timeStep))]
            logger.debug(f"Steady state: using timestep {uts}")
            timeList = [
                (str(fmt(uts)), str(fmt(0))),
                (str(fmt(uts)), str(fmt(dispersionDuration))),
            ]

        elif dynamicType == self.toolkit.TIME_DYNAMIC:
            # Dynamic: shift times so dispersion starts at t=0.
            if timeStep is None:
                uts = TS[0]
            else:
                uts = TS[min(range(len(TS)), key=lambda i: abs(TS[i] - timeStep))]
            logger.debug(f"Dynamic: first timestep {uts}")
            timeList = [(str(fmt(x)), str(fmt(x - uts))) for x in TS if x >= uts]
            # Extend if flow ends before dispersion duration.
            if TS[-1] < dispersionDuration:
                timeList.append((str(fmt(TS[-1])), str(fmt(dispersionDuration))))

        else:
            raise ValueError(
                f"Invalid temporal type '{dynamicType}'. "
                f"Must be '{self.toolkit.TIME_STEADYSTATE}' or '{self.toolkit.TIME_DYNAMIC}'.")

        return timeList, uts

    def _checkAndRegisterDispersionInDB(self, useDBSupport, dispersionFlowFieldName,
                                         workflowGroup, originalFlow, flowData,
                                         dispersionDuration, overwrite):
        """Query DB for existing dispersion flow. Returns existing doc or None."""
        if not useDBSupport:
            return None

        logger = get_classMethod_logger(self, "_checkAndRegisterDispersionInDB")
        querydict = dict(
            groupName=workflowGroup,
            flowParameters=dict(
                flowFields=flowData.get("dispersionFields", {}),
                dispersionDuration=dispersionDuration,
                originalFlow=originalFlow
            )
        )
        logger.debug(f"Querying DB for: {json.dumps(querydict, indent=4)}")
        docList = self.toolkit.getSimulationsDocuments(
            type=self.toolkit.DOCTYPE_OF_FLOWDISPERSION,
            **dictToMongoQuery(querydict),
            workflowName=dispersionFlowFieldName
        )

        if len(docList) > 1:
            raise ValueError(
                f"Found {len(docList)} documents for {dispersionFlowFieldName}. Fix manually.")

        if len(docList) == 0 or overwrite:
            if len(docList) == 1:
                logger.info("Overwriting — deleting old document")
                docList[0].delete()
            else:
                # Check if name exists with different parameters.
                docList = self.toolkit.getSimulationsDocuments(
                    type=self.toolkit.DOCTYPE_OF_FLOWDISPERSION,
                    workflowName=dispersionFlowFieldName
                )
                if len(docList) > 0:
                    return docList[0]
            return None

        return docList[0]

    def _prepareDispersionDirectory(self, destDir, origCaseDir, isParallel,
                                     timeList, linkMesh, linkData, flowData, flowType):
        """Copy/link directories and write dispersion-specific fields."""
        logger = get_classMethod_logger(self, "_prepareDispersionDirectory")
        os.makedirs(destDir, exist_ok=True)

        # Copy constant, system, 0 from original case.
        for subdir in ["constant", "system", "0"]:
            src = os.path.join(origCaseDir, subdir)
            dst = os.path.join(destDir, subdir)
            if os.path.exists(dst):
                shutil.rmtree(dst)
            logger.debug(f"Copying {src} → {dst}")
            shutil.copytree(src, dst)

        origDirsList = glob.glob(os.path.join(origCaseDir, "processor*")) if isParallel else [origCaseDir]

        for orig_time, dest_time in timeList:
            # Create root-level time directory (stochastic solver workaround).
            os.makedirs(os.path.join(destDir, str(dest_time)), exist_ok=True)
            logger.info(f"Mapping {orig_time} → {dest_time}")

            for origDir in origDirsList:
                self._copyOrLinkTimestep(
                    origDir, destDir, orig_time, dest_time,
                    isParallel, linkData, linkMesh
                )

            # Write dispersion-specific fields (e.g. ustar, Hmix).
            for fieldName, value in flowData.get("dispersionFields", {}).items():
                logger.info(f"Writing field {fieldName} = {value} at t={dest_time}")
                field = self.toolkit.OFObjectHome.getEmptyFieldFromCase(
                    fieldName=fieldName, flowType=flowType,
                    internalValue=value, caseDirectory=destDir
                )
                field.writeToCase(caseDirectory=destDir, timeOrLocation=str(dest_time))

    def _copyOrLinkTimestep(self, origDir, destDir, origTime, destTime,
                             isParallel, linkData, linkMesh):
        """Copy or symlink a single timestep directory, plus mesh handling."""
        logger = get_classMethod_logger(self, "_copyOrLinkTimestep")
        orig_ts = os.path.join(origDir, origTime)
        dest_ts = (os.path.join(destDir, os.path.basename(origDir), destTime)
                   if isParallel else os.path.join(destDir, destTime))

        if os.path.exists(dest_ts):
            shutil.rmtree(dest_ts)

        if linkData:
            os.makedirs(dest_ts, exist_ok=True)
            for f in glob.glob(os.path.join(orig_ts, "*")):
                os.symlink(os.path.abspath(f), os.path.join(dest_ts, os.path.basename(f)))
        else:
            shutil.copytree(orig_ts, dest_ts)

        # Handle mesh: copy or symlink polyMesh.
        proc_prefix = [os.path.basename(origDir)] if isParallel else []
        if not linkMesh:
            src = os.path.join(origDir, "constant")
            dst = os.path.join(destDir, *proc_prefix, "constant")
            if not os.path.exists(dst):
                shutil.copytree(src, dst)
        else:
            src = os.path.abspath(os.path.join(origDir, "constant", "polyMesh"))
            dst = os.path.join(destDir, *proc_prefix, "constant", "polyMesh")
            if not os.path.exists(dst):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                os.symlink(src, dst)

    def _finalizeDispersionInDB(self, useDBSupport, doc, name, group,
                                 directory, originalFlow, flowData, duration):
        """Register or update dispersion flow in the database."""
        if not useDBSupport:
            return None

        logger = get_classMethod_logger(self, "_finalizeDispersionInDB")
        querydict = dict(
            groupName=group,
            workflowName=name,
            flowParameters=dict(
                flowFields=flowData.get("dispersionFields", {}),
                dispersionDuration=duration,
                originalFlow=originalFlow
            )
        )

        if doc is None:
            logger.info("Adding new dispersion flow to database")
            self.toolkit.addSimulationsDocument(
                resource=directory,
                type=self.toolkit.DOCTYPE_OF_FLOWDISPERSION,
                dataFormat=datatypes.STRING,
                desc=querydict
            )
            return directory
        else:
            logger.info(f"Updating existing document. Returning {doc.resource}")
            doc.desc = querydict
            doc.save()
            return doc.resource

    def createDispersionCaseDirectory(self, hermes_dispersionWorkflow, updateDB=True, exportFromDB=False,
                                      allowDuplicate=False, rewrite=False):
        """
        Create a dispersion case directory and synchronise with the database.

        Ensures consistency between the input workflow, DB records, and disk.
        Handles name conflicts, parameter mismatches, and directory conflicts.

        Parameters
        ----------
        hermes_dispersionWorkflow : workflow_StochasticLagrangianSolver
            The workflow instance to create a case for.
        updateDB : bool
            On inconsistency, update DB from the input workflow.
        exportFromDB : bool
            On inconsistency, update the input workflow from DB.
        allowDuplicate : bool
            Allow adding even if same workflow exists under a different name.
        rewrite : bool
            Remove and recreate the directory if it exists and is inconsistent.

        Returns
        -------
        None
        """
        logger = get_classMethod_logger(self, "createDispersionCaseDirectory")

        # Validate flags — updateDB and exportFromDB are mutually exclusive.
        if updateDB and exportFromDB:
            raise ValueError("Cannot use both updateDB and exportFromDB.")

        if hermes_dispersionWorkflow.name is None:
            raise ValueError("Must set the name property in the dispersionWorkflow object.")

        dispersionDirectoryName = os.path.abspath(
            os.path.join(self.toolkit.filesDirectory, hermes_dispersionWorkflow.name))

        # Step 1: Resolve the dispersion flow field from DB.
        dispersionFlowFieldName, dispersionFlowFieldDirectory = self._resolveDispersionFlowField(
            hermes_dispersionWorkflow
        )

        # Step 2: Check DB consistency — does a workflow with this name or
        # these parameters already exist? Determine what action to take.
        hermes_dispersionWorkflow, updateWorkflow, needDBAdd, foundInconsistency = \
            self._checkDBConsistency(
                hermes_dispersionWorkflow, updateDB, exportFromDB, allowDuplicate
            )

        # Step 3: Handle existing directory on disk.
        self._resolveDirectoryConflict(
            dispersionDirectoryName, foundInconsistency, rewrite
        )

        # Step 4: Create the dispersion case and link to the flow field.
        logger.info(f"Creating dispersion case {dispersionDirectoryName}")
        self.toolkit.stochasticLagrangian.createAndLinkDispersionCaseDirectory(
            dispersionDirectoryName,
            dispersionFlowDirectory=dispersionFlowFieldDirectory
        )

        # Step 5: Add or update DB record.
        self._syncDispersionToDB(
            hermes_dispersionWorkflow, dispersionDirectoryName,
            dispersionFlowFieldName, needDBAdd, updateWorkflow
        )

    def _resolveDispersionFlowField(self, hermes_dispersionWorkflow):
        """Look up the dispersion flow field in DB, updating the workflow's reference if needed.

        Returns
        -------
        tuple of (str, str)
            (flowFieldName, flowFieldDirectory)
        """
        logger = get_classMethod_logger(self, "_resolveDispersionFlowField")
        flowFieldRef = hermes_dispersionWorkflow.dispersionFlowFieldName

        # Try the name first, then the absolute path as fallback.
        docList = self.toolkit.getCaseListDocumentFromDB(flowFieldRef)
        if len(docList) == 0:
            docList = self.toolkit.getCaseListDocumentFromDB(os.path.abspath(flowFieldRef))
            if len(docList) == 0:
                raise ValueError(
                    f"Dispersion flow field '{flowFieldRef}' not found in DB. Add it first."
                )

        # If the workflow used a directory path, update to the canonical DB name.
        dbName = docList[0].desc['workflowName']
        if flowFieldRef != dbName:
            logger.debug(f"Updating flow field reference: {flowFieldRef} → {dbName}")
            hermes_dispersionWorkflow.dispersionFlowFieldName = dbName

        return dbName, docList[0].resource

    def _checkDBConsistency(self, workflow, updateDB, exportFromDB, allowDuplicate):
        """Check DB for existing workflows with the same name or parameters.

        Returns
        -------
        tuple of (workflow, updateWorkflow, needDBAdd, foundInconsistency)
        """
        logger = get_classMethod_logger(self, "_checkDBConsistency")
        updateWorkflow = False
        foundInconsistency = None
        needDBAdd = False

        # Check 1: Does a workflow with the same NAME exist in DB?
        doc_nameQueryList = self.toolkit.getCaseListDocumentFromDB(workflow.name)

        if len(doc_nameQueryList) > 0:
            # Same name found — compare parameters to detect inconsistency.
            doc_similarWorkflow = self.toolkit.getHemresWorkflowFromDocument(doc_nameQueryList[0])
            res = self.toolkit.compareWorkflowsObj([doc_similarWorkflow, workflow])

            if len(res.columns) == 1:
                # Identical — no action needed.
                logger.info("Workflow in DB is identical to input.")
            else:
                # Parameters differ — resolve based on flags.
                logger.debug(f"DB workflow differs from input:\n{res}")
                if exportFromDB:
                    workflow = doc_similarWorkflow
                elif updateDB:
                    updateWorkflow = True
                    foundInconsistency = res
                else:
                    raise ValueError(
                        f"Workflow '{workflow.name}' exists in DB with different parameters. "
                        f"Use updateDB=True or exportFromDB=True."
                    )
        else:
            # Check 2: Does a workflow with the same PARAMETERS exist under a different name?
            otherWorkflows = self.toolkit.getHermesWorkflowFromDB(workflow)
            if otherWorkflows and len(otherWorkflows) > 0:
                foundNames = [x.desc['workflowName'] for x in otherWorkflows]
                if not allowDuplicate:
                    raise FileExistsError(
                        f"Workflow exists under names {foundNames}. "
                        f"Use allowDuplicate=True to add anyway."
                    )
                needDBAdd = True

        return workflow, updateWorkflow, needDBAdd, foundInconsistency

    def _resolveDirectoryConflict(self, directoryPath, foundInconsistency, rewrite):
        """Handle existing directory: remove if rewrite=True and inconsistent."""
        logger = get_classMethod_logger(self, "_resolveDirectoryConflict")

        if not os.path.isdir(directoryPath):
            return  # No conflict.

        if foundInconsistency is not None:
            # Directory exists but is inconsistent with the chosen workflow.
            if rewrite:
                logger.info(f"Rewrite flag: removing {directoryPath}")
                shutil.rmtree(directoryPath)
            else:
                raise ValueError(
                    f"Directory {directoryPath} exists with inconsistent workflow. "
                    f"Use rewrite=True to recreate."
                )

    def _syncDispersionToDB(self, workflow, directoryName, flowFieldName,
                             needDBAdd, updateWorkflow):
        """Add a new DB document or update an existing one."""
        logger = get_classMethod_logger(self, "_syncDispersionToDB")

        if needDBAdd:
            # New workflow — extract group/ID from naming convention.
            groupName = workflow.name.split("_")[0]
            groupID = workflow.name.split("_")[1]
            logger.info(f"Adding new document: group={groupName}, id={groupID}")
            self.toolkit.addSimulationsDocument(
                resource=directoryName,
                dataFormat=datatypes.STRING,
                type=self.toolkit.DOCTYPE_WORKFLOW,
                desc=dict(
                    groupName=groupName,
                    groupID=groupID,
                    workflowName=flowFieldName,
                    workflowType=workflow.workflowType,
                    workflow=workflow.json,
                    parameters=workflow.parametersJSON,
                ),
            )
        elif updateWorkflow:
            # Existing workflow — update parameters in DB.
            docList = self.toolkit.getCaseListDocumentFromDB(workflow.name)
            if len(docList) > 0:
                doc = docList[0]
                doc.desc['workflow'] = workflow.json
                doc.desc['parameters'] = workflow.parametersJSON
                logger.info("Updating workflow in DB")
                doc.save()

    def createAndLinkDispersionCaseDirectory(self, dispersionDirectory, dispersionFlowDirectory):
        """Create a dispersion case directory and link it to the flow field directory."""
        logger = get_classMethod_logger(self, "createAndLinkDispersionCaseDirectory")
        dispersionDirectory = os.path.abspath(dispersionDirectory)
        dispersionFlowDirectory = os.path.abspath(dispersionFlowDirectory)
        logger.info(f"Creating dispersion at {dispersionDirectory} with base flow {dispersionFlowDirectory}")

        if not os.path.isdir(dispersionFlowDirectory):
            err = f"The {dispersionFlowDirectory} is not not a directory."
            logger.error(err)
            raise ValueError(err)

        constantDir = os.path.join(dispersionDirectory, "constant")
        systemDir = os.path.join(dispersionDirectory, "system")

        logger.debug(f"Copying constant and system from the base flow")
        if os.path.exists(constantDir):
            shutil.rmtree(constantDir)

        if os.path.exists(systemDir):
            shutil.rmtree(systemDir)

        shutil.copytree(os.path.join(dispersionFlowDirectory, "constant"), constantDir)
        shutil.copytree(os.path.join(dispersionFlowDirectory, "system"), systemDir)

        for proc in glob.glob(os.path.join(dispersionFlowDirectory, "processor*")):
            logger.debug(f"Working on processor {proc}")
            fullpath = os.path.abspath(os.path.join(proc, "constant", "polyMesh"))
            destination = os.path.join(dispersionDirectory, os.path.basename(proc), "constant", "polyMesh")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            logger.debug(f"\t Linking: ln -s {fullpath} {destination}")
            logger.debug(
                f"\t Linking root case : ln -s {os.path.abspath(proc)} {os.path.join(dispersionDirectory, os.path.basename(proc))}/rootCase")
            if not os.path.exists(destination):
                os.symlink(fullpath, destination)
            rootcase_link = os.path.join(dispersionDirectory, os.path.basename(proc), "rootCase")
            if not os.path.exists(rootcase_link):
                os.symlink(os.path.abspath(proc), rootcase_link)

            # create the 0 directory in all processors.
            os.makedirs(os.path.join(dispersionDirectory, os.path.basename(proc), '0'), exist_ok=True)
            logger.debug(f"Creating the 0 time in {os.path.join(dispersionDirectory, os.path.basename(proc), '0')} ")

        # linking the decpomposePar dict from the root.
        # root_decomposePar = os.path.abspath(os.path.join(dispersionFlow,"system","decomposeParDict"))
        # decompose_dest    = os.path.abspath(os.path.join(dispersionDirectory,"system"))
        # os.system(f"ln -s {root_decomposePar} {decompose_dest}")

        # linking the rootCase in the root directory of the dispersion dispersionDirectory.
        logger.debug(
            f"Linking the root case: ln -s {dispersionFlowDirectory} {os.path.join(dispersionDirectory, 'rootCase')}")
        rootcase_root = os.path.join(dispersionDirectory, 'rootCase')
        if not os.path.exists(rootcase_root):
            os.symlink(dispersionFlowDirectory, rootcase_root)

        # create the 0 directory in the root.
        logger.debug(f"Making the 0 in {os.path.join(dispersionDirectory, '0')}")
        os.makedirs(os.path.join(dispersionDirectory, '0'), exist_ok=True)

        logger.info("... Done")

    @property
    def sourcesTypeList(self):
        """List of available source geometry type names."""
        return [x.split("_")[1] for x in dir(self) if "makeSource" in x and "_" in x]

    def writeParticlePositionFile(self, x, y, z, nParticles, dispersionName, type="Point",
                                  fileName="kinematicCloudPositions", **kwargs):
        """
            Writes a source position files.
            Saves the file to the constant directory of the case.
        Parameters:
        -----------
            x: float
                The x coordinate of the source.
            y: float
                The y coordinate of the source.
            z: float
             The z coordinate of the source.
            nParticles: int
                The number of particles.
            type: str
                The type of the source. Must be one of the sources.sourcesTypeList
            fileName: str
                The output filename (without the directory).
            kwargs:
                additional parameters for the source creation.

        Returns
        -------
            None


        """
        string = "/*--------------------------------*- C++ -*----------------------------------*\n " \
                 "| =========                 |                                                 |\n" \
                 "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
                 "|  \\    /   O peration     | Version:  dev                                   |\n" \
                 "|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n" \
                 "|    \\/     M anipulation  |                                                 |\n" \
                 "\\*---------------------------------------------------------------------------*/\n" \
                 "FoamFile\n{    version     2.0;\n    format      ascii;\n    class       vectorField;\n" \
                 "    object      kinematicCloudPositions;\n}\n" \
                 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n" \
                 f"{nParticles}\n(\n"

        slist = self.sourcesTypeList
        if type not in slist:
            raise ValueError(f"The type must be [{','.join(slist)}]. Got {type} instead")

        source = getattr(self, f"makeSource_{type}")(x=x, y=y, z=z, nParticles=nParticles, **kwargs)
        for i in range(nParticles):
            string += f"({source.loc[i].x} {source.loc[i].y} {source.loc[i].z})\n"
        string += ")\n"
        with open(os.path.join(dispersionName, "constant", fileName), "w") as writeFile:
            writeFile.write(string)

    def makeSource_Point(self, x, y, z, nParticles, **kwargs):
        """Generate particles at a single point location."""
        return pandas.DataFrame(
            {"x": [x for i in range(nParticles)], "y": [y for i in range(nParticles)],
             "z": [z for i in range(nParticles)]})

    def makeSource_Circle(self, x, y, z, nParticles, radius, distribution="uniform", **kwargs):
        """Generate particles distributed in a circle on the XY plane."""
        Rs = getattr(numpy.random, distribution)(0, radius, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        for R, theta in zip(Rs, thetas):
            xs.append(x + R * numpy.cos(theta))
            ys.append(y + R * numpy.sin(theta))
        return pandas.DataFrame({"x": xs, "y": ys, "z": [z for i in range(nParticles)]})

    def makeSource_Sphere(self, x, y, z, nParticles, radius, distribution="uniform", **kwargs):
        """Generate particles distributed within a sphere."""
        Rs = getattr(numpy.random, distribution)(0, radius, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        phis = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta, phi in zip(Rs, thetas, phis):
            xs.append(x + R * numpy.sin(theta) * numpy.cos(phi))
            ys.append(y + R * numpy.sin(theta) * numpy.sin(phi))
            zs.append(z + R * numpy.cos(theta))
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})

    def makeSource_Cylinder(self, x, y, z, nParticles, radius, height, horizontalDistribution="uniform",
                            verticalDistribution="uniform", **kwargs):
        """Generate particles distributed within a cylinder."""
        Rs = getattr(numpy.random, horizontalDistribution)(0, radius, nParticles)
        Hs = getattr(numpy.random, verticalDistribution)(-height / 2, height / 2, nParticles)
        thetas = numpy.random.uniform(0, 2 * numpy.pi, nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta, H in zip(Rs, thetas, Hs):
            xs.append(x + R * numpy.cos(theta))
            ys.append(y + R * numpy.sin(theta))
            zs.append(z + H)
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})

    def makeSource_Rectangle(self, x, y, z, nParticles, lengthX, lengthY, rotateAngle=0, **kwargs):
        """Generate particles distributed within a rotated rectangle on the XY plane."""
        xdist = numpy.random.uniform(-lengthX / 2, lengthX / 2, nParticles)
        ydist = numpy.random.uniform(-lengthY / 2, lengthY / 2, nParticles)
        xs = []
        ys = []
        for i in range(nParticles):
            xs.append(x + xdist[i] * numpy.cos(rotateAngle) + ydist[i] * numpy.sin(rotateAngle))
            ys.append(y - xdist[i] * numpy.sin(rotateAngle) + ydist[i] * numpy.cos(rotateAngle))
        return pandas.DataFrame({"x": xs, "y": ys, "z": [z for i in range(nParticles)]})

    def makeSource_Cube(self, x, y, z, nParticles, lengthX, lengthY, lengthZ, rotateAngle=0, **kwargs):
        """Generate particles distributed within a rotated cuboid volume."""
        xdist = numpy.random.uniform(-lengthX / 2, lengthX / 2, nParticles)
        ydist = numpy.random.uniform(-lengthY / 2, lengthY / 2, nParticles)
        zdist = numpy.random.uniform(-lengthZ / 2, lengthZ / 2, nParticles)
        xs = []
        ys = []
        zs = []
        for i in range(nParticles):
            xs.append(x + xdist[i] * numpy.cos(rotateAngle) + ydist[i] * numpy.sin(rotateAngle))
            ys.append(y - xdist[i] * numpy.sin(rotateAngle) + ydist[i] * numpy.cos(rotateAngle))
            zs.append(z + zdist[i])
        return pandas.DataFrame({"x": xs, "y": ys, "z": zs})

    def getMassFromLog(self, logFile, solver="StochasticLagrangianSolver"):
        """Parse mass and parcel fate information from a solver log file.

        The log file has a repeating structure per timestep:
        1. "Time = <t>" line marks a new timestep
        2. Cloud name line marks the start of cloud output
        3. For each source: release count + mass, then parcel fate (escape/stick)
        4. "ExecutionTime" marks end of timestep output
        """
        count = 0
        times = []
        names = []
        actions = []
        mass = []
        parcels = []

        def addToLists(time, name, action, m, parcel):
            """Append a parsed record to the accumulator lists."""
            times.append(time)
            names.append(name)
            actions.append(action)
            mass.append(m)
            parcels.append(parcel)

        with open(logFile, "r") as readFile:
            Lines = readFile.readlines()

        # Phase 1: Skip header lines until we find the first solver execution line.
        # The format is "Exec <solver>" — skip 2 more lines after that.
        for line in Lines:
            count += 1
            if "Exec" in line and solver in line:
                count += 2
                break

        # Phase 2: Parse timesteps until "End" marker.
        # State machine: outer loop finds "Time" lines, inner loop parses
        # cloud output until "ExecutionTime".
        while "End" not in Lines[count]:
            if "Time" in Lines[count]:
                time = float(Lines[count].split()[-1])
                # Scan forward to the cloud name line.
                while f"{self._datalayer.cloudName}\n" not in Lines[count]:
                    count += 1
                count += 1
                # Parse cloud output: each source has release info + parcel fate.
                while "ExecutionTime" not in Lines[count]:
                    if "Parcel fate" in Lines[count]:
                        # Parcel fate block: next 2 lines are escape and stick counts.
                        name = Lines[count].split()[-1]
                        count += 1
                        addToLists(time=time, name=name, action="escape", parcel=float(Lines[count].split()[-2][:-1]),
                                   m=float(Lines[count].split()[-1]))
                        count += 1
                        addToLists(time=time, name=name, action="stick", parcel=float(Lines[count].split()[-2][:-1]),
                                   m=float(Lines[count].split()[-1]))
                    elif ":\n" in Lines[count]:
                        # Release block: source name, then parcel count, then mass.
                        name = Lines[count].split()[0][:-1]
                        count += 1
                        parcel = float(Lines[count].split()[-1])
                        count += 1
                        addToLists(time=time, name=name, action="release", parcel=parcel,
                                   m=float(Lines[count].split()[-1]))
                    count += 1
            count += 1
        return pandas.DataFrame(
            {"time": times, "cloudName": self._datalayer.cloudName, "name": names, "action": actions,
             "parcels": parcels,
             "mass": mass})


    def getOriginalFlowFieldMesh(self,nameOrWorkflowFileOrJSONOrResource,readParallel=True, time=0):
        """
            Returns the mesh of the original flow field. name from the workflow

        Parameters
        ----------
        nameOrWorkflowFileOrJSONOrResource : string or dict
        The name/dict that defines the item

        readParallel: bool
                If parallel case exists, read it .

        time : float
            The time to read the mesh from. (relevant for mesh moving cases).

        Returns
        -------

        """
        logger = get_classMethod_logger(self,"getMeshFromLagrangianName")
        logger.info(f"Getting the mesh for {nameOrWorkflowFileOrJSONOrResource}")
        logger.debug(f"Getting the original flow field")
        originalFlowField = self.getOriginalFlowDocument(nameOrWorkflowFileOrJSONOrResource)
        logger.debug(f"Getting the mesh from the original flow field")
        return self.toolkit.getMesh(originalFlowField.getData())


    def getCaseResults(self, caseDescriptor, timeList=None, withVelocity=True, withReleaseTimes=False, withMass=True,
                       cloudName="kinematicCloud", forceSingleProcessor=False, cache=True, overwrite=False):
        """
        Read Lagrangian cloud data from disk, with caching support.

        Checks for cached results first; if not cached (or overwrite=True),
        reads particle data from the OpenFOAM case (parallel or serial),
        optionally caches to parquet.

        Parameters
        ----------
        caseDescriptor : str or MetadataFrame
            Case name, directory path, or DB document.
        timeList : list, optional
            Specific timesteps to read. If None, reads all.
        withVelocity : bool
            Include velocity fields.
        withReleaseTimes : bool
            Include particle age.
        withMass : bool
            Include particle mass.
        cloudName : str
            OpenFOAM cloud name (default: ``"kinematicCloud"``).
        forceSingleProcessor : bool
            Force serial read even if parallel directories exist.
        cache : bool
            Save results to DB cache.
        overwrite : bool
            Force recomputation even if cache exists.

        Returns
        -------
        dask.dataframe.DataFrame
        """
        logger = get_classMethod_logger(self, "getCaseResults")

        # Step 1: Resolve the case descriptor to a name string.
        caseDescriptorName = self._resolveCaseDescriptorName(caseDescriptor)

        # Step 2: Check cache — return early if cached and not overwriting.
        cacheDoc, cachedData = self._checkCache(
            caseDescriptorName, self.DOCTYPE_LAGRANGIAN_CACHE, overwrite
        )
        if cachedData is not None:
            return cachedData

        # Step 3: Locate the case directory (DB or filesystem).
        finalCasePath, workflowName = self._locateCaseDirectory(caseDescriptorName, caseDescriptor)

        # Step 4: Build the record loader for Lagrangian particle data.
        loader = lambda timeName: readLagrangianRecord(
            timeName, casePath=finalCasePath,
            withVelocity=withVelocity, withReleaseTimes=withReleaseTimes,
            cloudName=cloudName, withMass=withMass
        )

        # Step 5: Discover timesteps and load data via Dask.
        
        ret = self._loadCaseDataViaDask(
            finalCasePath, loader, timeList, forceSingleProcessor
        )

        # Step 6: Cache results if requested.
        if cache:
            ret = self._saveToCacheParquet(
                ret, cacheDoc, caseDescriptor, workflowName,
                cloudName, self.DOCTYPE_LAGRANGIAN_CACHE, datatypes.PARQUET
            )

        return ret

    # ------------------------------------------------------------------
    # Shared helpers for getCaseResults / getCaseConcentrationsEulerian
    # ------------------------------------------------------------------

    def _resolveCaseDescriptorName(self, caseDescriptor):
        """Extract a string name from a case descriptor (str or MetadataFrame)."""
        if isinstance(caseDescriptor, str):
            return caseDescriptor
        elif isinstance(caseDescriptor, MetadataFrame):
            return caseDescriptor.desc['workflowName']
        else:
            raise ValueError("caseDescriptor must be a string or MetadataFrame.")

    def _checkCache(self, caseDescriptorName, doctype, overwrite):
        """Check if cached results exist in the DB.

        Returns
        -------
        tuple of (cacheDoc or None, cachedData or None)
            If cache hit and not overwriting, returns (doc, data).
            If cache miss or overwriting, returns (doc_or_None, None).
        """
        logger = get_classMethod_logger(self, "_checkCache")
        cachedDocumentList = self.toolkit.getWorkflowDocumentFromDB(
            caseDescriptorName, doctype=doctype, dockind=self.toolkit.DOCKIND_CACHE
        )

        if len(cachedDocumentList) == 0:
            return None, None

        if len(cachedDocumentList) > 1:
            raise ValueError(
                f"Multiple cache entries for {caseDescriptorName}. Remove duplicates."
            )

        cacheDoc = cachedDocumentList[0]

        if overwrite:
            # Delete the cached file on disk so it will be recomputed.
            logger.info("Overwrite requested — removing cached file.")
            if os.path.isdir(cacheDoc.resource):
                shutil.rmtree(cacheDoc.resource)
            elif os.path.isfile(cacheDoc.resource):
                os.remove(cacheDoc.resource)
            return cacheDoc, None

        # Try to load from cache.
        try:
            data = cacheDoc.getData()
            logger.info("Returning cached data.")
            return cacheDoc, data
        except FileNotFoundError:
            # Cache document exists but file is missing — clean up.
            logger.error("Cached file not found. Removing stale cache document.")
            for doc in cachedDocumentList:
                doc.delete()
            return None, None

    def _locateCaseDirectory(self, caseDescriptorName, caseDescriptor):
        """Find the OpenFOAM case directory from DB or filesystem.

        Returns
        -------
        tuple of (str, str)
            (casePath, workflowName)
        """
        logger = get_classMethod_logger(self, "_locateCaseDirectory")
        docList = self.toolkit.getWorkflowDocumentFromDB(
            caseDescriptorName,
            doctype=self.toolkit.DOCTYPE_WORKFLOW,
            dockind=self.toolkit.DOCKIND_SIMULATIONS
        )

        if len(docList) == 0:
            # Not in DB — try as a filesystem path.
            casePath = os.path.abspath(str(caseDescriptor))
            if not os.path.isdir(casePath):
                raise FileNotFoundError(
                    f"{casePath} is not a directory and not in the DB."
                )
            logger.info(f"Found as directory: {casePath}")
            return casePath, str(caseDescriptor)

        logger.info(f"Found in DB: {docList[0].resource}")
        flowCaseDir= docList[0].resource
        workflowName = docList[0].desc['workflowName']
        if docList[0].type == "hermesWorkflow":
            flowCaseDir = os.path.join(os.path.dirname(flowCaseDir), workflowName)
        return flowCaseDir, workflowName

    def _discoverTimesteps(self, casePath, processorDir=None):
        """Discover numeric time directories, filtering out system dirs.

        Parameters
        ----------
        casePath : str
            Base case directory.
        processorDir : str, optional
            Processor subdirectory (e.g. "processor0"). If None, scan casePath directly.

        Returns
        -------
        list of str
            Sorted time directory names.
        """
        scanDir = os.path.join(casePath, processorDir) if processorDir else casePath
        excluded = {"constant", "system", "rootCase", "VTK"}
        return sorted(
            [x for x in os.listdir(scanDir)
             if (os.path.isdir(os.path.join(scanDir, x))
                 and x.isdigit()
                 and not x.startswith("processor")
                 and x not in excluded)],
            key=lambda x: int(x)
        )

    def _loadCaseDataViaDask(self, casePath, loader, timeList, forceSingleProcessor):
        """Load OpenFOAM case data using Dask distributed map.

        Detects parallel vs serial case, discovers timesteps if needed,
        builds the loader argument list, and returns a Dask DataFrame.

        Parameters
        ----------
        casePath : str
            Path to the OpenFOAM case.
        loader : callable
            Function that reads a single timestep (takes timeName string).
        timeList : list or None
            Specific timesteps. If None, discover all.
        forceSingleProcessor : bool
            Force serial read.

        Returns
        -------
        dask.dataframe.DataFrame
        """
        logger = get_classMethod_logger(self, "_loadCaseDataViaDask")

        isParallel = (os.path.exists(os.path.join(casePath, "processor0"))
                      and not forceSingleProcessor)

        if isParallel:
            logger.info("Processing as parallel case")
            processorList = [os.path.basename(p)
                             for p in glob.glob(os.path.join(casePath, "processor*"))]
            if timeList is None:
                timeList = self._discoverTimesteps(casePath, processorList[0])
            # Cartesian product: every (processor, timestep) combination.
            loaderList = [os.path.join(proc, t)
                          for proc, t in product(processorList, timeList)]
        else:
            logger.info("Processing as single-processor case")
            if timeList is None:
                timeList = self._discoverTimesteps(casePath)
            loaderList = list(timeList)

        logger.debug(f"Loading {len(loaderList)} items")
        return dask_dataframe.from_delayed(self.daskClient.map(loader, loaderList))

    def _saveToCacheParquet(self, data, cacheDoc, caseDescriptor,
                             workflowName, cloudName, doctype, dataFormat):
        """Save Dask DataFrame to parquet cache and register in DB."""
        logger = get_classMethod_logger(self, "_saveToCacheParquet")

        if cacheDoc is not None:
            fullname = cacheDoc.resource
        else:
            targetDir = os.path.join(
                self.toolkit.filesDirectory, "cachedLagrangianData", workflowName
            )
            os.makedirs(targetDir, exist_ok=True)
            fullname = os.path.join(targetDir, f"{cloudName}.parquet")
            self.toolkit.addCacheDocument(
                dataFormat=dataFormat, type=doctype,
                resource=fullname,
                desc={"workflowName": caseDescriptor, "cloudName": cloudName}
            )

        logger.info(f"Writing to parquet: {fullname}")
        data.set_index("datetime").repartition(partition_size="100MB").to_parquet(fullname)
        return dask.dataframe.read_parquet(fullname, engine='pyarrow')

    def _saveToCacheNetCDF(self, data, cacheDoc, caseDescriptor,
                            workflowName, cloudName, doctype):
        """Save Dask DataFrame as xarray netCDF cache and register in DB."""
        logger = get_classMethod_logger(self, "_saveToCacheNetCDF")

        if cacheDoc is not None:
            fullname = cacheDoc.resource
        else:
            targetDir = os.path.join(
                self.toolkit.filesDirectory, "cachedLagrangianData", workflowName
            )
            os.makedirs(targetDir, exist_ok=True)
            fullname = os.path.join(targetDir, f"{cloudName}ConcentrationEulerian.nc")
            self.toolkit.addCacheDocument(
                dataFormat=datatypes.NETCDF_XARRAY, type=doctype,
                resource=fullname,
                desc={"workflowName": caseDescriptor, "cloudName": cloudName}
            )

        logger.info(f"Writing to netCDF: {fullname}")
        # Group by grid cell, sum concentrations from all processors, convert to xarray.
        ret = data.groupby(["datetime", "x", "y", "z"]).sum().compute().to_xarray().fillna(0)
        ret.to_netcdf(fullname)
        return ret

    def getDispersionDocument(self, nameOrDispersionWorkflow):
        """Return the DB document for a dispersion simulation.

        Parameters
        ----------
        nameOrDispersionWorkflow : str or workflow_StochasticLagrangianSolver
            The name of the workflow or a workflow instance.

        Returns
        -------
        document or None
            The DB document, or None if not found.
        """
        # TODO: This method was incomplete in the original code (bare reference
        # to getWorkflowDocumentFromDB without calling it). Implemented as a
        # simple delegation to the toolkit's workflow document lookup.
        return self.toolkit.getWorkflowDocumentFromDB(nameOrDispersionWorkflow)

    def getDispersionFlowDocument(self,nameOrDispersionWorkflow):
        """
            Returns the DB document of the dispersion workflow.
            We assume that it is a name or a nameOrDispersionWorkflow.

        Parameters
        ----------
        nameOrWorkflow : str, workflow_StochasticLagrangianSolver
                The name of the workflow or an instance of the hermes workflow of the StochasticLagrangianSolver).

        Returns
        -------
            DB document.
        """
        logger = get_classMethod_logger(self,"getDispersionFlowDocument")
        if isinstance(nameOrDispersionWorkflow,str):
            wf = self.toolkit.getHermesWorkflowFromDB(nameOrDispersionWorkflow)
            dffname = wf.dispersionFlowFieldName
        elif isinstance(nameOrDispersionWorkflow,workflow_StochasticLagrangianSolver):
            dffname = nameOrDispersionWorkflow.dispersionFlowFieldName
        else:
            err = f"{nameOrDispersionWorkflow} must be of type str or workflow_StochasticLagrangianSolver, got {type(nameOrDispersionWorkflow)}"

        logger.info(f"Trying to retireve the document for {dffname}")
        ret = self.toolkit.getWorkflowDocumentFromDB(dffname, doctype=self.toolkit.DOCTYPE_OF_FLOWDISPERSION)
        return ret[0] if len(ret) > 0 else None

    def getOriginalFlowDocument(self,nameOrDispersionWorkflow):
        """
            Returns the flow document of the original workflow.

        Parameters
        ----------
        nameOrWorkflow : str, workflow_StochasticLagrangianSolver
                The name of the workflow or an instance of the hermes workflow of the StochasticLagrangianSolver).

        Returns
        -------
            DB document.
        """
        logger = get_classMethod_logger(self,"getDispersionFlowDocument")
        if isinstance(nameOrDispersionWorkflow,str):
            wf = self.toolkit.getHermesWorkflowFromDB(nameOrDispersionWorkflow)
            dffname = wf.originalFlowFieldName

        elif isinstance(nameOrDispersionWorkflow,workflow_StochasticLagrangianSolver):
            dffname = nameOrDispersionWorkflow.originalFlowFieldName
        else:
            err = f"{nameOrDispersionWorkflow} must be of type str or workflow_StochasticLagrangianSolver, got {type(nameOrDispersionWorkflow)}"

        ret = self.toolkit.getWorkflowDocumentFromDB(dffname, doctype=self.toolkit.DOCTYPE_WORKFLOW)
        return ret[0] if len(ret) > 0 else None

    def getOriginalFlowFieldExtent(self,nameOrDispersionWorkflow):
        """
            Returns the extends of the mesh of the original flow field.

        Parameters
        ----------
        nameOrDispersionWorkflow : str, workflow_StochasticLagrangianSolver
                The name of the workflow or an instance of the hermes workflow of the StochasticLagrangianSolver).


        Returns
        -------

        """
        logger = get_classMethod_logger(self,"getOriginalFlowFieldExtent")
        logger.info(f"Getting the mesh for {nameOrDispersionWorkflow}")
        logger.debug(f"Getting the original flow field")
        originalFlowField = self.getOriginalFlowDocument(nameOrDispersionWorkflow)
        logger.debug(f"Getting the mesh from the original flow field")
        mesh = self.toolkit.getMesh(originalFlowField.getData())
        logger.debug(f"Computing the extents. ")
        return mesh.getDataFrame()[['Cx','Cy','Cz']].agg(["min","max"])

    def getOriginalFlowFieldExtentAsDict(self,nameOrDispersionWorkflow):
        """
            Returns the extents as dict with the keys xmin,xmax,ymin,ymax,zmin,zmax

        Parameters
        ----------
        nameOrDispersionWorkflow : str
            The name of the nameOrDispersionWorkflow

        Returns
        -------

        """
        lims = self.getOriginalFlowFieldExtent(nameOrDispersionWorkflow=nameOrDispersionWorkflow)
        return dict(
            xmin = lims.Cx.loc['min']*ureg.m,
            xmax = lims.Cx.loc['max']*ureg.m,
            ymin = lims.Cy.loc['min']*ureg.m,
            ymax = lims.Cy.loc['max']*ureg.m,
            zmin = lims.Cz.loc['min']*ureg.m,
            zmax = lims.Cy.loc['max']*ureg.m
        )



    def getCaseConcentrationsEulerian(self,
                                      caseDescriptor,
                                      timeList=None,
                                      overwrite=False,
                                      cloudName="kinematicCloud",
                                      cache=True,
                                      forceSingleProcessor=False):
        """
        Read Eulerian concentration fields from disk, with caching support.

        Reads the CSV output of the concentrationField cloud function.
        Uses the same cache/load/compute pattern as ``getCaseResults``
        but stores results as netCDF (xarray) instead of parquet.

        Parameters
        ----------
        caseDescriptor : str or MetadataFrame
            Case name, directory path, or DB document.
        timeList : list, optional
            Specific timesteps. If None, reads all.
        overwrite : bool
            Force recomputation.
        cloudName : str
            OpenFOAM cloud name.
        cache : bool
            Save results to DB cache.
        forceSingleProcessor : bool
            Force serial read.

        Returns
        -------
        xarray.Dataset
            Concentration field indexed by (datetime, x, y, z).
        """
        logger = get_classMethod_logger(self, "getCaseConcentrationsEulerian")

        # Step 1: Resolve case descriptor.
        caseDescriptorName = self._resolveCaseDescriptorName(caseDescriptor)

        # Step 2: Check cache.
        cacheDoc, cachedData = self._checkCache(
            caseDescriptorName, self.DOCTYPE_CONCENTRATIONEULERIAN_CACHE, overwrite
        )
        if cachedData is not None:
            return cachedData

        # Step 3: Locate case directory.
        finalCasePath, workflowName = self._locateCaseDirectory(caseDescriptorName, caseDescriptor)

        # Step 4: Build Eulerian concentration loader.
        loader = lambda timeName: readEulerianConcentration(
            timeName, casePath=finalCasePath, cloudName=cloudName
        )

        # Step 5: Load data via Dask.
        ret = self._loadCaseDataViaDask(
            finalCasePath, loader, timeList, forceSingleProcessor, self.daskClient
        )

        # Step 6: Cache as netCDF (group by grid cell, sum concentrations).
        if cache:
            ret = self._saveToCacheNetCDF(
                ret, cacheDoc, caseDescriptor, workflowName,
                cloudName, self.DOCTYPE_CONCENTRATIONEULERIAN_CACHE
            )

        return ret



class analysis:
    """Analysis utilities for computing Lagrangian concentration fields."""

    DOCTYPE_CONCENTRATION = "xarray_concentration"
    DOCTYPE_CONCENTRATION_POINTWISE = "dask_concentration"

    datalayer = None

    def __init__(self, datalayer):
        """Initialize with a reference to the parent toolkit extension."""
        self.datalayer = datalayer
        self.datalayer.toolkit.initConfig(analysisFullMeshCounter=0,
                                          analysisPointWiseCounter=0)

        self.datalayer.toolkit.defineCounter("cartesianMeshCounter")


    def calcConcentrationPointWise(self, data, dxdydz, xfield="x", yfield="y", zfield="z"):
        """
            Calculates the concentration in cells from point data.

            The mass of each particle is given in the mass fields.

            Note that it is a multi index dask that cannot be saved as is to parquet.

        :param data: parquet/dask
            The particle location data in time.

        :param dxdydz: float
                The size of the cell

        :param xfield: str
                The name of the X coordinate

        :param yfield: str
                The name of the Y coordinate

        :param zfield: str
                The name of the Z coordinate
        :return: dask.


        """

        dH = dxdydz ** 3
        return data.assign(xI=dxdydz * (data[xfield] // dxdydz),
                           yI=dxdydz * (data[yfield] // dxdydz),
                           zI=dxdydz * (data[zfield] // dxdydz)).groupby(["xI", "yI", "zI", "time"])[
                   'mass'].sum().to_frame("C") / dH

    def calcDocumentConcentrationPointWise(self, dataDocument, dxdydz, xfield="x", yfield="y",
                                           zfield="z", overwrite=False, saveAsDask=False, simulationID=None,
                                           **metadata):
        """
           Calculates the concentration from the cells where particles exists.

        :param dataDocument: hera.MetadataDocument
                The document with the particles.

        :param dxdydz: float
                The size of the cell

        :param xfield: str
                The name of the X coordinate

        :param yfield: str
                The name of the Y coordinate

        :param zfield: str
                The name of the Z coordinate

        :param overwrite: bool
                If true, recalculates.

        :param saveAsDask: bool
                If true, reset the index because dask cannot read multi index (the result ofthe concentration is multi index with time and x,y,z)

        :param simulationID: str
                If not None, use this str instead of the docmentID. used in cases where the
                cache was transferred from another DB.

        :param metadata:
                Any additional parameters to add to the cache.
        :return:
                Document
        """
        simID = str(dataDocument.id) if simulationID is None else simulationID

        docList = self.datalayer.getCacheDocuments(simID=simID, dxdydz=dxdydz, **metadata)
        if len(docList) == 0 or overwrite:

            data = dataDocument.getData()
            C = self.calcConcentrationPointWise(data, dxdydz=dxdydz, xfield=xfield, yfield=yfield, zfield=zfield)

            if saveAsDask:
                C = C.compute().reset_index(level=['xI', 'yI', 'zI'])

            if len(docList) == 0:
                desc = dict(simID=simID, dxdydz=dxdydz)
                desc.update(metadata)

                doc = self.datalayer.addCacheDocument(dataFormat=datatypes.PARQUET,
                                                      type=self.DOCTYPE_CONCENTRATION_POINTWISE,
                                                      resource="",
                                                      desc=desc)

                outputDirectory = os.path.join(self.datalayer.toolkit.filesDirectory, str(doc.id))
                os.makedirs(outputDirectory, exist_ok=True)

                finalFileName = os.path.join(outputDirectory, "concentration.parquet")
                doc.resource = finalFileName
                doc.save()

            else:
                doc = docList[0]
                finalFileName = doc.resource

            C.to_parquet(finalFileName, compression="GZIP")
            ret = doc

    def calcConcentrationTimeStepFullMesh(self, timeData, extents, dxdydz, xfield="x", yfield="y",
                                          zfield="z"):
        """
            Converts a xyz particle data (with mass field) to a concentration field in the requested domain (defined by extent).

            Embed in a larget timestep.

        Parameters
        -----------

        datimeDatata: pandas dataframe.
            The data to convert. Takes only one time step at a time (for multiple time steps use iterConcentationField)

        extents: dict
            The domain in which the concentation will be calculated.
            has keys: xmin,xmax,ymin,ymax,zmin,zmax of the entire domain.

        dxdydz: float
                    The mesh steps.

        dxdydz: float
                    The size of a mesh unit (that will be created for the concentration).

        xfield: str
                The column name of the x coordinates.

        yfield: str
            The column name of the y coordinates.

        zfield: str
            The column name of the z coordinates.


        Returns
        --------

            Xarray.dataframe
        """
        dxdydz = tonumber(dxdydz,ureg.m)
        x_full = numpy.arange(tonumber(extents['xmin'],ureg.m), tonumber(extents['xmax'],ureg.m)+1, dxdydz)
        y_full = numpy.arange(tonumber(extents['ymin'],ureg.m), tonumber(extents['ymax'],ureg.m)+1, dxdydz)
        z_full = numpy.arange(tonumber(extents['zmin'],ureg.m), tonumber(extents['zmax'],ureg.m)+1, dxdydz)

        dH = dxdydz ** 3

        fulldata = xarray.DataArray(coords=dict(xI=x_full, yI=y_full, zI=z_full), dims=['xI', 'yI', 'zI']).fillna(0)
        #fulldata.filterType = "C"
        if not isinstance(timeData,int):
            C = timeData.assign(xI=dxdydz * (timeData[xfield] // dxdydz), yI=dxdydz * (timeData[yfield] // dxdydz),
                                zI=dxdydz * (timeData[zfield] // dxdydz)).groupby(["xI", "yI", "zI", "datetime"])[
                    'mass'].sum().to_xarray().squeeze().fillna(0) / dH

            # assign the timestep into the large mesh
            fulldata.loc[dict(xI=C.xI, yI=C.yI, zI=C.zI)] = C
            fulldata.attrs['field'] = "1*kg/m**3"

            timeList = [timeData.datetime.unique()[0]]
        else:
            timeList = [timeData]

        return fulldata.expand_dims(dict(datetime=timeList), axis=-1)

    def calcConcentrationFieldFullMesh(self, caseDescriptor, dxdydz, extents=None, xfield="x", yfield="y",
                                       zfield="z", overwrite=False, reReadResults=None, **metadata):
        """
        Compute Eulerian concentration on a Cartesian mesh for each timestep.

        Bins Lagrangian particle positions onto a regular grid, stores results
        as partitioned netCDF files. Pads with zero-concentration timesteps
        up to the dispersion duration for moving-window averaging.

        Parameters
        ----------
        caseDescriptor : str or MetadataFrame
            Case name, path, or DB document.
        dxdydz : float
            Grid cell size in meters.
        extents : dict, optional
            Domain bounds ``{xmin, xmax, ymin, ymax, zmin, zmax}``.
            If None, derived from the original flow field mesh.
        xfield, yfield, zfield : str
            Column names for coordinates.
        overwrite : bool
            Force recomputation even if cache exists.
        reReadResults : bool, optional
            Force re-read of Lagrangian results. Defaults to *overwrite*.
        **metadata
            Additional metadata for the cache document.

        Returns
        -------
        xarray.Dataset
            Concentration field with ``C`` variable and ``dt`` attribute.
        """
        logger = get_classMethod_logger(self, "calcConcentrationFieldFullMesh")

        # Step 1: Resolve case name and check for existing cache.
        caseDescriptorName = self._resolveCaseDescriptorName(caseDescriptor)

        mdata = dict(extents=extents, dxdydz=dxdydz, caseDescriptorName=caseDescriptorName)
        mdata.update(**metadata)
        mdata = ConfigurationToJSON(mdata)
        docList = self.datalayer.toolkit.getCacheDocuments(type=self.DOCTYPE_CONCENTRATION, **mdata)

        if len(docList) == 0 or overwrite:
            # Step 2: Clean up old cache if overwriting.
            if len(docList) > 0:
                self._removeConcentrationCache(docList[0])

            # Step 3: Create new cache document with a counter-based path.
            xryDoc = self._createConcentrationCacheDoc(caseDescriptorName, mdata)
            path_to_data = os.path.dirname(xryDoc.resource)
            os.makedirs(path_to_data, exist_ok=True)

            # Step 4: Resolve extents from flow field if not provided.
            if extents is None:
                extents = self.datalayer.getOriginalFlowFieldExtentAsDict(caseDescriptorName)
                logger.info(f"Extents from flow field: {json.dumps(ConfigurationToJSON(extents), indent=4)}")

            # Step 5: Load Lagrangian particle data.
            reReadResults = overwrite if reReadResults is None else reReadResults
            workflow = self.datalayer.toolkit.getHermesWorkflowFromDB(caseDescriptorName)
            data = self.datalayer.getCaseResults(
                caseDescriptorName, overwrite=reReadResults,
                withVelocity=True, withReleaseTimes=True, withMass=True
            )

            # Step 6: Process each Dask partition — bin particles onto mesh per timestep.
            timeName = 0
            partitionID = 0
            for partitionID, partition in enumerate(data.partitions):
                timeName = self._processConcentrationPartition(
                    partition, partitionID, path_to_data, extents, dxdydz,
                    xfield, yfield, zfield
                )

            # Step 7: Pad remaining timesteps with zero concentrations.
            # This ensures the time series extends to the full dispersion
            # duration, needed for moving-window risk averaging (assumes dt=1).
            self._padRemainingTimesteps(
                timeName, workflow.dispersionDuration, partitionID + 1,
                path_to_data, extents, dxdydz, xfield, yfield, zfield
            )
        else:
            xryDoc = docList[0]

        # Step 8: Load the result and attach units metadata.
        ret = xryDoc.getData()
        ret.attrs['field'] = dict(C=1*ureg.kg/ureg.m**3)
        ret.attrs['dt'] = f"{(ret.datetime[-1]-ret.datetime[-2]).item()}s"
        return ret

    def _removeConcentrationCache(self, cacheDoc):
        """Remove an existing concentration cache document and its files."""
        logger = get_classMethod_logger(self, "_removeConcentrationCache")
        logger.info(f"Removing old cache at {cacheDoc.resource}")
        for f in glob.glob(cacheDoc.resource):
            try:
                os.remove(f)
            except OSError as e:
                logger.error(f"Error removing {f}: {e.strerror}")
        cacheDoc.delete()
        cache_dir = os.path.dirname(cacheDoc.resource)
        if os.path.exists(cache_dir):
            shutil.rmtree(cache_dir)

    def _createConcentrationCacheDoc(self, caseDescriptorName, mdata):
        """Create a new cache document for concentration data."""
        newID = self.datalayer.toolkit.addCounter("cartesianMeshCounter")
        resourcePath = os.path.join(
            self.datalayer.toolkit.filesDirectory, "cachedLagrangianData",
            f"{caseDescriptorName}_fullMeshCache_{newID}", "Concentrations*.nc"
        )
        return self.datalayer.toolkit.addCacheDocument(
            dataFormat=datatypes.NETCDF_XARRAY,
            resource=resourcePath,
            type=self.DOCTYPE_CONCENTRATION,
            desc=mdata
        )

    def _processConcentrationPartition(self, partition, partitionID, outputDir,
                                        extents, dxdydz, xfield, yfield, zfield):
        """Process a single Dask partition: bin particles per timestep, write netCDF.

        Returns the last timeName processed (for padding calculation).
        """
        logger = get_classMethod_logger(self, "_processConcentrationPartition")
        logger.info(f"Processing partition {partitionID}")
        L = []
        timeName = 0

        for timeName, timeData in partition.compute().reset_index().groupby("datetime"):
            xry = self.calcConcentrationTimeStepFullMesh(
                timeData, extents=extents, dxdydz=dxdydz,
                xfield=xfield, yfield=yfield, zfield=zfield
            )
            L.append(xry)

        if len(L) > 0:
            self._writeConcentrationPartition(L, partitionID, outputDir)

        return timeName

    def _padRemainingTimesteps(self, lastTime, totalDuration, nextPartitionID,
                                outputDir, extents, dxdydz, xfield, yfield, zfield):
        """Fill remaining timesteps with zero-concentration fields."""
        logger = get_classMethod_logger(self, "_padRemainingTimesteps")
        L = []
        for t in range(int(lastTime + 1), totalDuration):
            xry = self.calcConcentrationTimeStepFullMesh(
                timeData=t, extents=extents, dxdydz=dxdydz,
                xfield=xfield, yfield=yfield, zfield=zfield
            )
            L.append(xry)

        if len(L) > 0:
            self._writeConcentrationPartition(L, nextPartitionID, outputDir)

    @staticmethod
    def _writeConcentrationPartition(xarray_list, partitionID, outputDir):
        """Concatenate timestep xarrays and write to a partitioned netCDF file."""
        pxry = xarray.concat(xarray_list, dim="datetime")
        outFile = os.path.join(outputDir, f"Concentrations{partitionID:04}.nc")
        # Rename grid indices to standard names and transpose to (y, x, z, time)
        # order expected by downstream consumers.
        pxry.rename(dict(xI="x", yI="y", zI="z")).transpose(
            "y", "x", "z", "datetime"
        ).to_dataset(name="C").to_netcdf(outFile)

    def getConcentrationField(self, dataDocument, returnFirst=True, **metadata):
        """
            Returns the concentration field that is related to this simulation (with the metadata that is provided).

            Return None if simulation is not found.

        Parameters
        ----------
        dataDocument: datalayer.metadatadocument

        returnFirst: bool
            If true, returns only the first simulation and None if none was found.
            If false returns a list (might be empty).

        metadata: any metadata that is needed

        Returns
        -------
            The document of the concentrations (or none).

        """
        dataID = str(dataDocument.id)
        docList = self.datalayer.getCacheDocuments(dataID=dataID, type=self.DOCTYPE_CONCENTRATION, **metadata)

        if returnFirst:
            if len(docList) > 0:
                data = docList[0]
            else:
                data = None

        else:
            data = docList

        return data


##########################
####   Utils
##########################

def robustOpenFOAMFileValuesParser(path, columnNames):
    # The sed command explained:
    # '1,18d': Deletes the header, comments, and the count/bracket lines.
    # 's/[()]//g': Globally removes all parentheses.
    # 's/^[[:space:]]*//; s/[[:space:]]*$//': Trims leading/trailing whitespace.
    # 's/[[:space:]]\+/,/g': Replaces any whitespace chunk with a comma.
    # '/^$/d': Deletes empty lines.
    # '/\//d': Deletes lines containing '//' (foam footers).
    
    if not os.path.exists(path):
        return pandas.DataFrame(columns=columnNames).astype(float)

    START_OF_FILE_VALUES = 18
    FILE_SAMPLE_SIZE = 2048

    sed_command = [
        "sed",
        "-e", f"1,{START_OF_FILE_VALUES}d",
        "-e", r"s/[()]//g",
        "-e", r"s/^[[:space:]]*//",
        "-e", r"s/[[:space:]]*$//",
        "-e", r"s/[[:space:]]\+/,/g",
        "-e", r"/^$/d",
        "-e", r"/\/\//d",
        path,
    ]
    with open(path, 'r') as f:
        content = f.read(FILE_SAMPLE_SIZE)

    uniform_match = re.search(r'(\d+)\{(.*)\}', content)
    if uniform_match:
        # print("hit", uniform_match.group(1), uniform_match.group(2))
        count = int(uniform_match.group(1))
        # Clean up the value string (remove brackets if vector)
        val_str = uniform_match.group(2).replace('(', '').replace(')', '')
        single_val = numpy.array([float(x) for x in val_str.split()])

        data = numpy.tile(single_val, (count, 1))
        return pandas.DataFrame(data, columns=columnNames).astype(float)

    # failing in the process should give an unexpected error meaning we missed some case
    proc = subprocess.run(sed_command, capture_output=True, text=True, check=True)
    
    if not proc.stdout.strip():
        return pandas.DataFrame(columns=columnNames).astype(float)

    return pandas.read_csv(
        io.StringIO(proc.stdout),
        sep=',',
        header=None,
        names=columnNames,
        engine='c',
        dtype=float
    ).astype(float)




def readLagrangianRecord(timeName, casePath, withVelocity=False, withReleaseTimes=False, withMass=True,
                         cloudName="kinematicCloud"):
    """Read a single Lagrangian time step record from the case directory."""

    columnsDict = dict(x=[], y=[], z=[], id=[], procId=[], globalID=[],datetime=[])
    if withVelocity:
        columnsDict['U_x'] = []
        columnsDict['U_y'] = []
        columnsDict['U_z'] = []
    if withReleaseTimes:
        columnsDict["releaseTime"] = []
        columnsDict['age'] = []
    if withMass:
        columnsDict['mass'] = []

    newData = pandas.DataFrame(columnsDict, dtype=numpy.float64)

    newData = robustOpenFOAMFileValuesParser(
        os.path.join(casePath, timeName, "lagrangian", cloudName, "globalPositions"),
        ['x', 'y', 'z'])
    for fld in ['x', 'y', 'z']:
        newData[fld] = newData[fld].astype(numpy.float64)

    dataID = robustOpenFOAMFileValuesParser(os.path.join(casePath, timeName, "lagrangian", cloudName, "origId"),
                            ['id'])
    newData['id'] = dataID['id'].astype(numpy.float64)

    dataprocID = robustOpenFOAMFileValuesParser(
        os.path.join(casePath, timeName, "lagrangian", cloudName, "origProcId"), ['procId'])
    newData['procId'] = dataprocID['procId'].astype(numpy.float64)

    newData = newData.ffill().assign(globalID=1000000000 * newData.procId + newData.id)

    # dataGlobal = extractFile(
    #     os.path.join(casePath, timeName, "lagrangian", cloudName, "globalPositions"),
    #     ['globalX', 'globalY', 'globalZ'])
    #
    # for col in ['globalX', 'globalY', 'globalZ']:
    #     newData[col] = dataGlobal[col].astype(numpy.float64)

    theTime = os.path.split(timeName)[-1]
    newData['datetime'] = float(theTime)


    if withVelocity:
        dataU = robustOpenFOAMFileValuesParser(os.path.join(casePath, timeName, "lagrangian", cloudName, "U"),
                            ['U_x', 'U_y', 'U_z'])
        for col in ['U_x', 'U_y', 'U_z']:
            newData[col] = dataU[col]

    if withReleaseTimes:
        dataM = robustOpenFOAMFileValuesParser(os.path.join(casePath, timeName, "lagrangian", cloudName, "age"),
                            ['age'])
        newData["releaseTime"] = newData["datetime"] - dataM["age"]
        newData["age"] = dataM["age"]
    if withMass:
        dataM = robustOpenFOAMFileValuesParser(os.path.join(casePath, timeName, "lagrangian", cloudName, "mass"),
                            ['mass'])
        try:
            newData["mass"] = dataM["mass"]
        except:
            newData = newData.compute()
            newData["mass"] = dataM["mass"]

    return newData.dropna()

def readEulerianConcentration(timeName, casePath,cloudName="kinematicCloud"):
    """
        Reads the concentration file that was outputed by the cloud function ConcentrationField
    Parameters
    ----------
    timeName : string
        The time to read
    casePath  :  string
        The case path.

    cloudName : string
        The name of the cloud.

    Returns
    -------

    """
    fileName = os.path.join(casePath, timeName, f"{cloudName}EulerConcentrations")
    if os.path.exists(fileName):
        try:
            return pandas.read_csv(fileName).assign(datetime=float(timeName.split("/")[-1]))
        except pandas.errors.EmptyDataError:
            return pandas.DataFrame(dict(x=[], y=[], z=[],  C=[],datetime=[]))
    else:
        return pandas.DataFrame(dict(x=[],y=[],z=[],C=[],datetime=[]))

