from .... import datalayer
from .... import GIS
from ... import WRF
import os
import pandas
import numpy
import geopandas

class process():

    _projectName = None
    _GISdatalayer = None

    def __init__(self, projectName):

        self._projectName = projectName
        self._GISdatalayer = GIS.GIS_datalayer(projectName=self._projectName, FilesDirectory="")

    def setTemplate(self, path):

        datalayer.Measurements.addDocument(projectName=self._projectName, resource=path, dataFormat="string", type="Template", desc={"Template":"OpenFoam"})

    def runCase(self, data, wrfZ=500, Time=None, wrfPath=None, velocityVector=[1,0,0], dirName="OFsimulation", newPath=None, dxdy=50, zMax=1000, xCells=50, yCells=50, zCells=50):
        """
        Runs an OF simulation. May use wrf results or a user's velocity vctor as inlet velocity values.

        Parameters
        ----------
        data: The data that should be converted to stl. May be a dataframe or a name of a saved polygon in the database.
        wrfZ: Used as a reference height.
        Time: A datetime.time object with time in hours and minutes, for wrf results.
        wrfPath: The path of the wrf results.
        velocityVector: The velocity components, if wrfPath is None.
        dirName: A name for the directory of the OF simulation, if None use "OFsimulation".
        newPath: The path for the new directory, if None uses current directory.
        dxdy: the dimention of each cell in the mesh in meters, the default is 50.
        zMax: The maximum height of the mesh, default is 1000.
        xCells: number of xCells, default is 50.
        yCells: number of yCells, default is 50.
        zCells: number of zCells, default is 50.
        -------

        """

        # Makes the stl
        newPath = newPath if newPath is not None else os.path.abspath(os.getcwd())
        os.system("cp -avr %s %s" % (datalayer.Measurements.getDocuments(projectName=self._projectName,Template="OpenFoam")[0].asDict()["resource"], newPath))
        os.system("mv %s/template %s/%s" % (newPath, newPath, dirName))
        GIScon = GIS.convert(projectName=self._projectName, FilesDirectory="%s/%s/constant/triSurface" % (newPath, dirName))
        stlstr, data = GIScon.toSTL(data=data, NewFileName="topo", dxdy=dxdy)
        data = data.reset_index()

        # Edit the blockMesh
        my_file = open('%s/%s/system/blockMeshDict' % (newPath, dirName))
        string_list = my_file.readlines()
        my_file.close()
        for key, value in zip(["xMin", "xMax", "yMin", "yMax", "zMin", "zMax", "xCells", "yCells", "zCells"],
                              [data.gridxMin[0],data.gridxMax[0],data.gridyMin[0],data.gridyMax[0],(data.gridzMin[0]-10),zMax,xCells,yCells,zCells]):
            for i in range(len(string_list)):
                if key in string_list[i]:
                    string_list[i] = "    %s   %s;\n" % (key, value)
                    break
        new_file_contents = "".join(string_list)
        my_file = open('%s/%s/system/blockMeshDict' % (newPath, dirName), "w")
        my_file.write(new_file_contents)
        my_file.close()

        # Edit the snappyHexMesh
        my_file = open('%s/%s/system/snappyHexMeshDict' % (newPath, dirName))
        string_list = my_file.readlines()
        my_file.close()
        for i in range(len(string_list)):
            if "locationInMesh" in string_list[i]:
                string_list[i] = '    locationInMesh   (%s %s %s);\n' % ((data.gridxMin[0] + data.gridxMax[0]) / 2, (data.gridyMin[0] + data.gridyMax[0]) / 2,(data.gridzMin[0] + 1000) / 2)
        new_file_contents = "".join(string_list)
        my_file = open('%s/%s/system/snappyHexMeshDict' % (newPath, dirName), "w")
        my_file.write(new_file_contents)
        my_file.close()

        #Find velocity components
        if wrfPath is not None:
            wrfdatalayer = WRF.wrfDatalayer()
            datalat = wrfdatalayer.getPandas(datapath=wrfPath, Time=Time, lon=data.gridxMin[0], heightLimit=wrfZ).query("LAT>=%s and LAT<=%s" % (data.gridyMin[0],data.gridyMax[0]))
            datalong = wrfdatalayer.getPandas(datapath=wrfPath, Time=Time, lat=data.gridyMax[0], heightLimit=wrfZ).query("LONG>=%s and LONG<=%s" % (data.gridxMin[0],data.gridxMax[0]))
            newData = pandas.concat([datalat.iloc[(datalat['height']-wrfZ).abs().argsort()[:2]],
                                     datalong.iloc[(datalong['height']-wrfZ).abs().argsort()[:2]]])
            U = newData.U.mean()
            V = newData.V.mean()
            W = newData.W.mean()
            height = newData.height.mean()
        else:
            U = velocityVector[0]
            V = velocityVector[1]
            W = velocityVector[2]
            height = wrfZ

        vel = numpy.sqrt(numpy.square(U) + numpy.square(V) + numpy.square(W))

        #Edit k, U and epsilon
        for p in ["k", "U", "epsilon"]:
            my_file = open('%s/%s/0/%s' % (newPath, dirName, p))
            string_list = my_file.readlines()
            my_file.close()
            for i in range(len(string_list)):
                if "flowDir" in string_list[i]:
                    string_list[i] = '        flowDir   (%s %s %s);\n' % (U, V, W)
                if "Zref" in string_list[i]:
                    string_list[i] = '        Zref   %s;\n' % height
                if "Uref" in string_list[i]:
                    string_list[i] = '        Uref   %s;\n' % vel
            new_file_contents = "".join(string_list)
            my_file = open('%s/%s/0/%s' % (newPath, dirName, p), "w")
            my_file.write(new_file_contents)
            my_file.close()

        # Run the simulation
        os.chdir("%s/%s" % (newPath, dirName))
        os.system("blockMesh")
        os.system("surfaceFeatureExtract")
        os.system(("snappyHexMesh -overwrite"))
        os.system("simpleFoam")


