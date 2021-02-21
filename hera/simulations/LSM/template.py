import os
# from ...datalayer.document.metadataDocument import Simulations as SimulationsDoc
from ...datalayer.document import getDBObject
from ..LSM import LagrangianReader
from .DataLayer import SingleSimulation
from ..utils.inputForModelsCreation import InputForModelsCreator
from hera.datalayer import Project
import shutil
import xarray

class LSMTemplate(Project):
    _document = None

    LSM_RUN_TYPE = 'LSM_run'

    def __init__(self, document,projectName,to_xarray=True,to_database=False,forceKeep=False):
        """
        Initializes the template object.

    Parameters
    ----------

        document: DB document (or JSON)
            contains the defaults parameters and the data under document['desc']

        projectName: str
            The project name of the current simulations.

        to_xarray: bool
            Save the simulation results into xarray or not

        to_database: bool
            Save the simulation run in the database or not

        forceKeep: bool
            If to_xarray is true, determine wehter to keep the original files.
            if False, removes the Lagrnagian files.



        """
        super().__init__(projectName=projectName)
        self._document = document

        self.to_xarray =to_xarray
        self.to_database = to_database
        self.forceKeep = forceKeep
        self.datetimeFormat = "timestamp"
        self.xColumn = "x"
        self.yColumn = "y"
        self.uColumn = "u"
        self.directionColumn = "direction"
        self.timeColumn = "datetime"


    @property
    def dirPath(self):
        return self._document['resource']

    @property
    def params(self):
        return self._document['desc']['params']

    @property
    def version(self):
        return self._document['desc']['version']

    @property
    def modelName(self):
        return self._document['type'].split('_')[0]

    @property
    def modelFolder(self):
        return self._document['desc']['modelFolder']

    def run(self, saveDir,topography=None, stations=None,overwrite=False,params=dict(),**descriptor):
        """
        Execute the LSM simulation

        Parameters
        ----------
        projectName: str
            The project name

        saveDir: str
            Path of the directory to put in the model run

        overwrite: bool
            False: execute the simulation event if it is in DB (the to_database is True).
            True : retrieve simulation from DB (if to_database is True)

            The retrieval will work only if the simulations were converted to xarray (to_xarray was true when they were calculated).

        params: dict
            overweriting the parameters of the simulation

        descriptors : a list of key/value that describe that simulations.

        Return:
            the xarray (if it is in the DB, or the simulation was converted to xarray)
            other wise None.

        """
        fileDict = {".true.":"OUTD3d03_3_",".TRUE.":"OUTD3d03_3_",".false.":"OUTD2d03_3_",".FALSE.":"OUTD2d03_3_"}
        saveDir = os.path.abspath(saveDir)

        # create the input file.
        # paramsMap['wind_dir'] = self.paramsMap['wind_dir_meteorological']
        self._document['desc']['params'].update(params)
        self._document['desc'].update(descriptor)
        if topography is None:
            self._document['desc']['params'].update(homogeneousWind=".TRUE.")
            print("setting homogeneous wind")
        else:
            # if stations is None:
            #     raise KeyError("When using topography stations must be delivered.")
            # else:
            self._document['desc']['params'].update(TopoFile="'TOPO'")
        if stations is not None:
            self._document['desc']['params'].update(homogeneousWind=".FALSE.",StationsFile="'STATIONS'")
        xshift = (self._document['desc']['params']["TopoXmax"] - self._document['desc']['params']["TopoXmin"]) * \
                 self._document['desc']['params']["sourceRatioX"]
        yshift = (self._document['desc']['params']["TopoYmax"] - self._document['desc']['params']["TopoYmin"]) * \
                 self._document['desc']['params']["sourceRatioY"]

        ifmc = InputForModelsCreator(self.dirPath) # was os.path.dirname(__file__)
        ifmc.setParamsMap(self._document['desc']['params'])
        ifmc.setTemplate('%s_%s' % (self.modelName, self.version))

        if self.to_database:
            docList = self.getSimulationsDocuments(type=self.LSM_RUN_TYPE,params=self._document['desc']['params'],version=self.version)
            if len(docList)==0:
                doc = self.addSimulationsDocument(
                                            type=self.LSM_RUN_TYPE,
                                            resource='None',
                                            dataFormat='None',
                                            desc=dict(params=self._document['desc']['params'],
                                                      version=self.version,
                                                      datetimeFormat=self.datetimeFormat
                                                      )
                                            )

                saveDir = os.path.join(saveDir, str(doc.id))
                if self.to_xarray:
                    doc['resource'] = os.path.join(saveDir, 'netcdf', '*')
                    doc['dataFormat'] = 'netcdf_xarray'
                else:
                    doc['resource'] = os.path.join(saveDir)
                    doc['dataFormat'] = 'string'
                doc.save()
            elif not overwrite:
                    return docList[0].getData()

        ## If overwrite, or document does not exist in DB, or running without DB.
        os.makedirs(saveDir, exist_ok=True)

        os.system('cp -rf %s %s' % (os.path.join(self.modelFolder, '*'), saveDir))
        # write to file.
        ifmc.render(os.path.join(saveDir, 'INPUT'))

        os.chdir(saveDir)
        if topography is not None:
            with open("TOPO","w") as topofile:
                topofile.write(topography)
            # make stations files

        if stations is not None:
            stations = stations.rename(columns={self.timeColumn:"datetime"})
            onlyStations = stations.drop(columns=["datetime",self.uColumn,self.directionColumn])
            stations[self.stationColumn] = stations[self.xColumn].astype(str) + stations[self.yColumn].astype(str)
            stationsXarray = stations.set_index([self.xColumn, self.yColumn, "datetime"]).drop(columns=self.stationColumn).to_xarray()
            resampled = stationsXarray.resample(datetime="5Min").interpolate()
            stations = resampled.to_dataframe().reset_index()
            stations = stations.set_index([self.xColumn,self.yColumn]).join(onlyStations.set_index([self.xColumn,self.yColumn])).reset_index().dropna()

            allStationsFile = f"{len(stations[self.stationColumn].drop_duplicates())}\n"
            i = 0
            j = 0
            for station in stations[self.stationColumn].drop_duplicates():
                stationName = f"{chr(65+i)}{chr(65+j)}"
                allStationsFile += f"{stationName}  {list(stations.loc[stations[self.stationColumn]==station][self.xColumn])[0]}"
                for k in range(9-len(str(list(stations.loc[stations[self.stationColumn]==station][self.xColumn])[0]))):
                    allStationsFile += " "
                allStationsFile += f"{list(stations.loc[stations[self.stationColumn]==station][self.yColumn])[0]}\n"
                stationFile = ""
                for k in range(2901):
                    stationFile += "0 0\n"
                for k, line in enumerate(stations.loc[stations[self.stationColumn]==station].iterrows()):
                    stationFile += f"{line[1][self.uColumn]} {line[1][self.directionColumn]}\n"
                with open(os.path.join("tozaot","Meteorology",f"{stationName}_st.txt"),"w") as newStationFile:
                    newStationFile.write(stationFile)
                if j == 25:
                    j = 0
                    i += 1
                else:
                    j += 1
            with open("STATIONS","w") as newStationFile:
                newStationFile.write(allStationsFile)

        # run the model.
        os.system('./a.out')
        if self.to_xarray:
            results_full_path = os.path.join(saveDir, "tozaot", "machsan", fileDict[self._document['desc']['params']["particles3D"]])
            netcdf_output = os.path.join(saveDir, "netcdf")
            os.makedirs(netcdf_output, exist_ok=True)

            L = []
            i = 0
            for xray in LagrangianReader.toNetcdf(basefiles=results_full_path,datetimeFormat=self.datetimeFormat):
                L.append(xray)

                if len(L) == 100:  # args.chunk:
                    finalxarray = xarray.concat(L, dim="datetime")
                    for var, shift in zip(["x", "y"], [xshift, yshift]):
                        finalxarray[var] -= shift
                    finalxarray.to_netcdf(os.path.join(netcdf_output, "data%s.nc" % i))
                    L = []
                    i += 1

            # save the rest.
            finalxarray = xarray.concat(L, dim="datetime")

            new_coords = dict(x=finalxarray.x-xshift,y=finalxarray.y-yshift)
            finalxarray= finalxarray.assign_coords(coords=new_coords)

            print("saved xarray in ",netcdf_output)
            if not self.forceKeep:
                machsanPath = os.path.dirname(results_full_path)
                allfiles = os.path.join(machsanPath ,"*")
                os.system(f"rm {allfiles}")


            finalxarray.to_netcdf(os.path.join(netcdf_output, "data%s.nc" % i))
            return finalxarray

    def getLSMRuns(self,**query):
        """
        get a list of SingleSimulation objects that fulfill the query
        :param query:
        :return:
        """
        docList = self.getSimulationsDocuments(type=self.LSM_RUN_TYPE,**query)
        return [SingleSimulation(doc) for doc in docList]