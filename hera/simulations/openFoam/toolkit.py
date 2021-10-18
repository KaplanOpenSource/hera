from ...toolkit import abstractToolkit
import pandas
import xarray
import numpy


class OFToolkit(abstractToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """

    def __init__(self,projectName,FilesDirectory=None,toolkitName="OFToolkit"):
        super().__init__(toolkitName=toolkitName,projectName=projectName,FilesDirectory=FilesDirectory)


    def getMesh(self,caseDirectory,parallel=True):
        """
            Reads the mesh from the mesh directory.

            Reads the decomposed case if it exists and parallel is true,
            otherwise, reads just the single case.

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
        pass




    def updateInternalFieldFromPandas(self,caseDirectory,time,fieldName,data,parallel=True,cols=None):
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

        fieldName: str
                The field name
        data: pandas
                The data to update. Must include the column index (for composed cases) and processorName/index for decomposed cases.

        parallel: bool
                If true, check if there is a decomposed case

        Returns
        -------
            String with the new data. if parallel, a dict with the processor number as key and the
            value as the data.

        """
        pass


    def runWorkflow(self,workflowNameOrJSON,saveMode):
        """
            Runs a workflow. Checks if it is in the DB first for the appropriate saveMode.


        :param workflowNameOrJSON:
        :return:
        """
        pass

