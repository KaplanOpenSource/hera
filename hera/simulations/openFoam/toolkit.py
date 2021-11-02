from ...toolkit import abstractToolkit
import os
import glob
from .OFObjects import OFField


class OFToolkit(abstractToolkit):
    """
        The goal of this toolkit is to provide the functions that are required to run workflows.
        and to mange the workflows in the DB.

        This toolkit might relay on the hermes project in order to manipulate the nodes
        of the workflow. (TBD).

    """

    def __init__(self, projectName, filesDirectory=None, toolkitName="OFToolkit"):
        super().__init__(toolkitName=toolkitName, projectName=projectName, filesDirectory=filesDirectory)


    def getFieldDimensions(self,fieldName):
        """
            Return the dimensions of the field.

            For now, we assume that the flow is incompressible.

        Parameters
        ----------
        fieldName : str
            The field name

        Returns
        -------
            str

        """
        dimensions = dict(U="[ 0 1 -1 0 0 0 0 ]",
                          p="[ 0 2 -2 0 0 0 0 ]",
                          epsilon="[ 0 2 -3 0 0 0 0 ]",
                          f="[0 0 -1 0 0 0 0]",
                          k="[0 2 -2 0 0 0 0]",
                          nut="[0 2 -1 0 0 0 0]")

        return dimensions[fieldName]



    def processorList(self,caseDirectory):
        """
            Returns the list of processors directories in the case
        Parameters
        ----------
        caseDirectory : str
            Path to the directory.

        Returns
        -------

        """
        return [os.path.basename(proc) for proc in glob.glob(os.path.join(caseDirectory, "processor*"))]


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
        self.logger.info(f"Start. case {caseDirectory}.")


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

        self.logger.info(f"End")

        cellCenters = OFField(name="C",dimensions="",componentNames=['x','y','z'])

        return cellCenters.load(caseDirectory,times=time,parallelCase=useParallel)


    def runWorkflow(self,workflowNameOrJSON,saveMode):
        """
            Runs a workflow. Checks if it is in the DB first for the appropriate saveMode.


        :param workflowNameOrJSON:
        :return:
        """
        pass



    def writeEmptyFieldFile(self,caseDirectory,time,fieldName,parallel=False,dimensions=None):
        """
            Writes an empty field file to the target case directory/time. If parallel, then write it in
            every time directory.

        Parameters
        ----------
        caseDirectory : str
            The case directory
        time  : str,int
            The time to write to.

        fieldName : str
            The field name

        parallel:  bool
            If true, use parallel case. (write to every processor*/time directory)

        dimensions : str
            If None, take from the default of the field.

        Returns
        -------

        """

        fileStr = """
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 1 -1 0 0 0 0 ];

internalField   uniform {value};

boundaryField
{
    "proc.*"
    {
        type            processor;
    }

}


// ************************************************************************* //
"""



