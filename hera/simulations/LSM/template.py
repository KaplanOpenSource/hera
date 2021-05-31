import os
from .singleSimulation import SingleSimulation
from ...datalayer import nonDBMetadataFrame,datatypes
from itertools import product
import glob
from ..utils.inputForModelsCreation import InputForModelsCreator
import xarray
import pandas
import numpy
from unum.units import *
from ... import toolkit
from ... utils.ConvertJSONtoConf import ConvertJSONtoConf

meterKeys = ["TopoXmin","TopoXmax","TopoYmin","TopoYmax","releaseHeight","inversionHeight","savedx","savedy","savedz"]
secondKeys = ["releaseDuration","savedt"]
minuteKeys = ["duration"]
velocityKeys = ["windSpeed"]

class LSMTemplate:
    _document = None
    _toolkit = None

    _config =None

    STABILITY_NEUTRAL = "neutral"
    STABILITY_STABLE = "stable"
    STABILITY_UNSTABLE = "unstable"

    @property
    def toolkit(self):
        return self._toolkit

    @property
    def doctype_simulation(self):
        return "LSM_run"

    def __init__(self, document,toolkit):
        """
        Initializes the template object.

    Parameters
    ----------

        document: DB document (or JSON)
            contains the defaults parameters and the data under document['desc']


        """

        self._document = document
        self._toolkit  = toolkit

        self.to_xarray      = toolkit.to_xarray
        self.to_database    = toolkit.to_database
        self.forceKeep      = toolkit.forceKeep
        self.datetimeFormat = "timestamp"
        self.xColumn        = "x"
        self.yColumn        = "y"
        self.uColumn        = "u"
        self.hColumn        = "h"
        self.directionColumn= "direction"
        self.timeColumn     = "datetime"
        self.stationColumn = "station"


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
    def templateName(self):
        return self._document['desc']['datasourceName']

    @property
    def modelFolder(self):
        return self._document['desc']['modelFolder']

    def run(self,topography=None, stations=None,canopy=None,params=dict(),depositionRates=None, saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE,**descriptor):
        """
        Execute the LSM simulation

        Parameters
        ----------
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
        saveDir = os.path.abspath(self.toolkit.FilesDirectory)

        # create the input file.

        updated_params = dict(self._document['desc']['params'])
        updated_params.update(params)
        updated_params.update(descriptor)
        updated_params = ConvertJSONtoConf(updated_params)
        for keys, unit in zip([minuteKeys,meterKeys,secondKeys,velocityKeys],[min,m,s,m/s]):
            for key in keys:
                updated_params[key] = updated_params[key].asNumber(unit)

        print(updated_params)

        if topography is None:
            updated_params.update(homogeneousWind=".TRUE.")
            if stations is None:
                print("setting homogeneous wind")
        else:
            updated_params.update(TopoFile="'TOPO'")

        if depositionRates is not None:
            if not isinstance(depositionRates,list):
                depositionRates = [depositionRates,depositionRates]
            updated_params.update(n_vdep=len(depositionRates))

        if stations is not None:
            updated_params.update(homogeneousWind=".FALSE.",StationsFile="'STATIONS'")
        if canopy is None:
            updated_params.update(canopy=".FALSE.")
        else:
            updated_params.update(canopy=".TRUE.")

        xshift = (updated_params["TopoXmax"] - updated_params["TopoXmin"]) * updated_params["sourceRatioX"]

        yshift = (updated_params["TopoYmax"] - updated_params["TopoYmin"]) * updated_params["sourceRatioY"]

        ifmc = InputForModelsCreator(self.dirPath) # was os.path.dupdated_paramsirname(__file__)
        ifmc.setParamsMap(updated_params)
        ifmc.setTemplate('LSM_%s' % (self.version))
        docList = self.toolkit.getSimulationsDocuments(type=self.doctype_simulation,
                                                       templateName=self.templateName,
                                                       version=self.version,**updated_params)
        print(f"Found {docList}")
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:
            if len(docList) == 0:
                if saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB:
                    raise ValueError(f"A run with requested parameters already exists in the databse; you may choose "
                                     f"{toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE} in order to replace it.")
                if saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE:
                    docList[0].delete()
            doc = self.toolkit.addSimulationsDocument(
                type=self.doctype_simulation,
                resource='None',
                dataFormat='None',
                desc=dict(version=self.version,
                          datetimeFormat=self.datetimeFormat,
                          templateName=self.templateName,
                          **updated_params)
            )

            saveDir = os.path.join(saveDir, str(doc.id))
            if self.to_xarray:
                doc['resource'] = os.path.join(saveDir, 'netcdf', '*')
                doc['dataFormat'] = datatypes.NETCDF_XARRAY
            else:
                doc['resource'] = os.path.join(saveDir)
                doc['dataFormat'] = datatypes.STRING
            doc.save()

        print(f"The saveDir is {saveDir}")

        if os.path.exists(os.path.join(saveDir,"netcdf")):
            if saveMode == toolkit.TOOLKIT_SAVEMODE_ONLYFILE:
                raise ValueError(f"The outputfile {os.path.join(saveDir,'netcdf')} exists. Either remove it or run a saveMode that ends with _REPLACE")
            elif saveMode == toolkit.TOOLKIT_SAVEMODE_ONLYFILE_REPLACE:
                os.system(f"rm -r {os.path.join(saveDir,'netcdf')}")

        ## If overwrite, or document does not exist in DB, or running without DB.
        os.makedirs(saveDir, exist_ok=True)

        os.system('cp -rf %s %s' % (os.path.join(self.modelFolder, '*'), saveDir))
        # write to file.
        ifmc.render(os.path.join(saveDir, 'INPUT'))

        os.chdir(saveDir)
        if topography is not None:
            with open("TOPO","w") as topofile:
                topofile.write(topography)
        if depositionRates is not None:
            depositionsFile = ""
            for rate in depositionRates:
                depositionsFile += f"{rate}\n"
            with open("INPUT_VDEP","w") as newDepositionFile:
                newDepositionFile.write(depositionsFile)
            # make stations files

        if stations is not None:
            stations = stations.rename(columns={self.timeColumn:"datetime"})
            onlyStations = stations.drop(columns=["datetime",self.uColumn,
                                                  self.directionColumn])
            if self.hColumn in onlyStations.columns:
                onlyStations = onlyStations.drop(columns=self.hColumn)

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
        if canopy is not None:
            with open("canopy_properties.txt","w") as newCanopyFile:
                newCanopyFile.write(canopy)
            if (topography is not None or stations is not None):
                if self.hColumn not in stations.columns:
                    raise KeyError("Height column is missing in the stations dataframe")
                else:
                    hStations = ""
                    for h in stations.drop_duplicates(self.stationColumn)[self.hColumn]:
                        hStations += f"{h} "
                    with open("h_stations.txt", "w") as newStationFile:
                        newStationFile.write(hStations)

        print("Running the model")
        # run the model.
        os.system('./a.out')
        if self.to_xarray:
            results_full_path = os.path.join(saveDir, "tozaot", "machsan", fileDict[updated_params["particles3D"]])
            netcdf_output = os.path.join(saveDir, "netcdf")
            os.makedirs(netcdf_output, exist_ok=True)

            L = []
            i = 0
            for xray in self._toNetcdf(basefiles=results_full_path, datetimeFormat=self.datetimeFormat):
                L.append(xray)

                if len(L) == 100:  # args.chunk:
                    finalxarray = xarray.concat(L, dim="datetime")
                    new_coords = dict(x=finalxarray.x - xshift, y=finalxarray.y - yshift)
                    finalxarray = finalxarray.assign_coords(coords=new_coords)
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

            if saveMode != toolkit.TOOLKIT_SAVEMODE_NOSAVE:
                finalxarray.to_netcdf(os.path.join(netcdf_output, "data%s.nc" % i))

            return nonDBMetadataFrame(finalxarray,**updated_params)
        else:
            return None

    def getLSMRuns(self,**query):
        """
        get a list of SingleSimulation objects that fulfill the query
        :param query:
        :return:
        """
        docList = self.getSimulationsDocuments(type=self.doctype_simulation,
                                               templateName = self.templateName,
                                               **query)

        return [SingleSimulation(doc) for doc in docList]


    def _toNetcdf(self, basefiles, addzero=True, datetimeFormat="timestamp"):
        """
            Converts the data to netcdf.
            The dosage are converted to s/m**3 instead of min/m**3.

            Parameters
            ----------
            basefiles: str
                Path to the directory with the netcdf files

            addZero: bool
                if true, adds a 0 file at the begining of the files (with time shift 0)

        """

        # outfilename = name
        filenameList = []
        times = []
        for infilename in glob.glob(os.path.join("%s*" % basefiles)):
            filenameList.append(infilename)
            times.append(float(infilename.split("_")[-1]))

        print("Processing the files")
        print(basefiles)
        print([x for x in glob.glob(os.path.join("%s*" % basefiles))])
        # Sort according to time.
        combined = sorted([x for x in zip(filenameList, times)], key=lambda x: x[1])
        dt = None

        for (i, curData) in enumerate(combined):
            print("\t... reading %s" % curData[0])
            cur = pandas.read_csv(curData[0], delim_whitespace=True,
                                  names=["y", "x", "z", "Dosage"])  # ,dtype={'x':int,'y':int,'z':int,'Dosage':float})

            cur['time'] = curData[1]
            if dt is None:
                dt = cur.iloc[0]['time'] * s

            cur['Dosage'] *= (s / m ** 3).asNumber(min / m ** 3)
            xray = cur.sort_values(['time', 'x', 'y', 'z']).set_index(['time', 'x', 'y', 'z']).to_xarray()
            if datetimeFormat == "timestamp":
                datetime = pandas.to_datetime("1-1-2016 12:00") + pandas.to_timedelta(xray.time.values, 's')
            elif datetimeFormat == "seconds":
                datetime = xray.time.values

            # finalxarray.to_netcdf(os.path.join(topath,name,"%s_%s.nc" % (outfilename, str(cur['time'].iloc[0]).replace(".", "_"))) )

            if (i == 0) and addzero:
                if datetimeFormat == "timestamp":
                    zdatetime = [pandas.to_datetime("1-1-2016 12:00")]
                elif datetimeFormat == "seconds":
                    zdatetime = [0]

                finalxarray = xarray.DataArray(numpy.zeros(xray['Dosage'].values.shape), \
                                               coords={'x': xray.x, 'y': xray.y, 'z': xray.z, 'datetime': zdatetime},
                                               dims=['datetime', 'y', 'x', 'z']).to_dataset(name='Dosage')

                yield finalxarray
            # finalxarray.to_netcdf(os.path.join(topath,name,"%s_0_0.nc" % outfilename) )

            finalxarray = xarray.DataArray(xray['Dosage'].values, \
                                           coords={'x': xray.x, 'y': xray.y, 'z': xray.z, 'datetime': datetime},
                                           dims=['datetime', 'y', 'x', 'z']).to_dataset(name='Dosage')

            yield finalxarray

    def getSimulations(self, **query):
        """
        get a list of SingleSimulation objects that fulfill the query

        Parameters
        ---------
        templateName : str
            The name of the template of the simulation

        query: parameters
            Parameters for the query

        Returns
        -------
            Simulation object
        """

        docList = self.getDocuments(type=self.doctype_simulation,
                                    templateName=self.templateName,
                                    **query)
        return [SingleSimulation(doc) for doc in docList]

    def getSimulationByID(self,id):
        """
        get a simulation by document id

        :param id:
        :return:
        """
        return SingleSimulation(self.getDocumentByID(id))

    def listSimulations(self,templateName, wideFormat=False, **query):
        """
            List the Simulation parameters that fulfil the query
        :param query:
        :return:
        """

        docList = self.getDocuments(type=self.doctype_simulation,
                                    templateName=self.templateName
                                    **query)
        descList = [doc.desc.copy() for doc in docList]

        for (i, desc) in enumerate(descList):
            desc.update({'id':docList[i].id})

        params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
        params_df_list = [df.rename(columns=dict([(x,"params__%s"%x) for x in df.columns])) for df in params_df_list]

        desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
        df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]

        new_df_list = []

        for df in df_list:
            id = df['id'][0]
            new_df = df.copy().drop(columns=['id']).melt()
            new_df.index = [id]*len(new_df)
            new_df_list.append(new_df)

        try:
            df = pandas.concat(new_df_list)
            if wideFormat:
                return df.pivot(columns='variable', values='value')
            else:
                return df
        except ValueError:
            raise FileNotFoundError('No simulations found')
