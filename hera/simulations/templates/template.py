import os
# from ...datalayer.document.metadataDocument import Simulations as SimulationsDoc
from ...datalayer.document import getDBObject
from ..LSM import LagrangianReader
from .inputForModelsCreation import InputForModelsCreator
import xarray


class LSMTemplate(object):
    _document = None

    def __init__(self, document):
        self._document = document
        pass

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

    def run(self, saveDir, projectName='LSM', to_xarray=True, to_database=True, **kwargs):
        """
        Execute the LSM simulation

        Parameters
        ----------
        projectName: str
            The project name

        saveDir: str
            Path of the directory to put in the model run

        to_xarray: bool
            Save the simulation results into xarray or not

        to_database: bool
            Save the simulation run in the database or not
        """

        # create the input file.
        # paramsMap['wind_dir'] = self.paramsMap['wind_dir_meteorological']
        self._document['desc']['params'].update(kwargs)
        print(self.dirPath)
        ifmc = InputForModelsCreator(self.dirPath) # was os.path.dirname(__file__)
        ifmc.setParamsMap(self.params)
        ifmc.setTemplate('%s_%s' % (self.modelName, self.version))

        if to_database:
            SimulationsDoc = getDBObject('Simulations')
            doc = dict(projectName=projectName,
                       type='%s_run' % self.modelName,
                       resource='None',
                       dataFormat='None',
                       desc=dict(params=self.params,
                                 version=self.version
                                 )
                       )
            doc = SimulationsDoc(**doc).save()
            saveDir = os.path.join(saveDir, str(doc.id))
            if to_xarray:
                doc['resource'] = os.path.join(saveDir, 'netcdf', '*')
                doc['dataFormat'] = 'netcdf_xarray'
            else:
                doc['resource'] = saveDir
                doc['dataFormat'] = 'string'

            doc.save()

        os.makedirs(saveDir, exist_ok=True)

        os.system('cp -rf %s %s' % (os.path.join(self.modelFolder, '*'), saveDir))
        # write to file.
        ifmc.render(os.path.join(saveDir, 'INPUT'))

        # run the model.
        os.chdir(saveDir)
        os.system('./a.out')

        if to_xarray:
            results_full_path = os.path.join(saveDir, "tozaot", "machsan", "OUTD2d03_3_")
            netcdf_output = os.path.join(saveDir, "netcdf")
            os.makedirs(netcdf_output, exist_ok=True)

            L = []
            i = 0
            for xray in LagrangianReader.toNetcdf(basefiles=results_full_path):
                L.append(xray)

                if len(L) == 100:  # args.chunk:
                    finalxarray = xarray.concat(L, dim="datetime")
                    finalxarray.to_netcdf(os.path.join(netcdf_output, "data%s.nc" % i))
                    L = []
                    i += 1


            # save the rest.
            finalxarray = xarray.concat(L, dim="datetime")
            finalxarray.to_netcdf(os.path.join(netcdf_output, "data%s.nc" % i))

