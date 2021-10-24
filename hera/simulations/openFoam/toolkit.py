from ...toolkit import abstractToolkit
from .readDaskParallel import loadEulerianDataParallel
import pandas
import xarray
import numpy
import os
import glob


class OFToolkit(abstractToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """

    def __init__(self, projectName, filesDirectory=None, toolkitName="OFToolkit"):
        super().__init__(toolkitName=toolkitName, projectName=projectName, filesDirectory=filesDirectory)


    def getMesh(self,caseDirectory,parallel=True,time=0):
        """
            Reads the mesh from the mesh directory.

            Reads the decomposed case if it exists and parallel is true,
            otherwise, reads just the single case.

            Unfortunately, we have to use the OF postProcess utility in order to interpolate the
            mesh points to their centers.

        Parameters
        ----------
            caseDirectory: str
                    The path to the case

            parallel: bool
                    If parallel case exists, read it .
        Returns
        -------
            pandas dataframe with the points in the columns             x,y,z
            the index column (don't mix up with the index of pandas)  is the sequential number of the point.

            If the case is decomposed, return processorNumber and index columns.
            The index is the internal order in the processor.
        """

        # 1. Run the postProcess utility to set the cell centers
        self.logger.info(f"Getting mesh for case {caseDirectory}.")


        useParallel= False
        cmd = f"postProcess -func writeCellCentres -case {caseDirectory}"
        if parallel:
            self.logger.debug(f"Attempt to load parallel case")
            # Check if the case is decomposed, if it is, run it.
            proc0dir = os.path.join(caseDirectory,"processor0")

            if os.path.exists(proc0dir):
                self.logger.debug(f"Found parallel case, using decomposed case")
                cmd = f"foamJob -parallel postProcess -func writeCellCentres -case {caseDirectory}"
                useParallel = True
            else:
                self.logger.debug(f"parallel case NOT found. Using composed case")

        os.system(cmd)

        return loadEulerianDataParallel(caseDirectory,fieldName="C",columnNames=['x','y','z'],times=time,parallelCase=useParallel)



    def updateInternalDataFromPandas(self,caseDirectory,time,fileName,data,parallel=True,cols=None,isField=True):
        """
            Updates the field internal data to the input.
            Should **not** update the boundaries data.

            use cols to write the data.

            If parallel is true, and there is a decomposed case, write the data to the different processors.
            In that case, the series must have processorNumber and index columns. Otherwise, only index columns is required.

            The new data is written to the disk

        Paramaters
        ----------
        caseDirectory: str
                The path of the case
        time: str
                The time to update.

        fileName: str
                The field name
        data: pandas
                The data to update. Must include the column index (for composed cases) and processorName/index for decomposed cases.

        cols : list
                The list of columns to write (and their order).

        parallel: bool
                If true, check if there is a decomposed case

        isField: bool

        Determine whether the fileName we update is a field (i.e with boundary conditions) or a list (i.e without).

                If true:
                    check if file exists and if it does use its boundary conditions.
                    if file does not exist, add the general boundary conditions with zerGradient for all possible boundary conditions.

        Returns
        -------
            String with the new data. if parallel, a dict with the processor number as key and the
            value as the data.

        """


        if parallel:
            procPaths = [proc for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]




        else:
            fullFileName = os.path.join(caseDirectory, time, fileName)


        if os.path.exists(fullFileName):




        else:
            # create a new file.
            newFileStr = self._getHeader()






    def runWorkflow(self,workflowNameOrJSON,saveMode):
        """
            Runs a workflow. Checks if it is in the DB first for the appropriate saveMode.


        :param workflowNameOrJSON:
        :return:
        """
        pass


    ##############################




    def _getHeader(self):

        return """
/*--------------------------------*- C++ -*----------------------------------*\
 =========                 |                                                 
 \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           
  \    /   O peration     | Version:  dev                                   
   \  /    A nd           | Web:      www.OpenFOAM.org                      
    \/     M anipulation  |                                                 
\*---------------------------------------------------------------------------*/
FoamFile
{    
    version     2.0;
    format      ascii;
    class       vectorField;
    object      kinematicCloudPositions;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

    def _pandasToFoamFormat(self,data,columns=None):

        D = data if columns is None else data[columns]

        newStr = f"{str(data.shape[0])}\n"
        newStr += "(\n"
        newStr += "\n".join([f"({x})" for x in D.to_csv(sep=' ', header=False, index=False).split("\n")[:-1]])
        newStr += ")"
        return newStr