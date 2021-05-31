import pandas
import json
import os

from .template import LSMTemplate,meterKeys,minuteKeys,secondKeys,velocityKeys
from itertools import product
from ... import datalayer
from ... import toolkit
from .singleSimulation import SingleSimulation
from unum.units import *

class LSMToolkit(toolkit.abstractToolkit):
    """
        The LSM toolkit

        The datasources are the templates.

        The LSM template JSON structure is :


        << TO DO >>



    """

    _to_xarray = None
    _to_database = None
    _forceKeep = None

    @property
    def to_xarray(self):
        return self._to_xarray

    @to_xarray.setter
    def to_xarray(self, value):
        if not isinstance(value,bool):
            raise ValueError("to_xarray must be boolean")
        self._to_xarray = value


    @property
    def to_database(self):
        return self._to_database

    @to_database.setter
    def to_database(self, value):
        if not isinstance(value,bool):
            raise ValueError("to_xarray must be boolean")
        self._to_database = value

    @property
    def forceKeep(self):
        return self._forceKeep

    @forceKeep.setter
    def forceKeep(self, value):
        if not isinstance(value, bool):
            raise ValueError("to_xarray must be boolean")
        self._forceKeep = value

    @property
    def singleSimulation(self):
        return SingleSimulation

    def __init__(self,projectName,FilesDirectory=None,to_xarray=True,to_database=False,forceKeep=False):
        """
            Initializes the LSM toolkit

        Parameters
        ----------

        projectName: str
            The name of the project that contains the

        FilesDirectory: str
            The directory to save the simulations.


        to_xarray: bool
            Save the simulation results into xarray or not

        to_database: bool
            Save the simulation run in the database or not

        forceKeep: bool
            If to_xarray is true, determine wehter to keep the original files.
            if False, removes the Lagrnagian files.

        """
        super().__init__(projectName=projectName,toolkitName="LSM",FilesDirectory=FilesDirectory)

        self.to_xarray = to_xarray
        self.to_database = to_database
        self.forceKeep = forceKeep


    def getTemplates(self, **query):
        """
            Get a list of Template objects that fulfill the query

        Returns
        -------
            List of LSMTemplates

        """
        docList = self.getDatasourceDocumentsList(**query)
        return [LSMTemplate(doc,self) for doc in docList]

    def getTemplateByName(self,templateName,templateVersion=None,to_xarray=True,to_database=False,forceKeep=False):
        """
            Retrieve the template by its name.

        Parameters
        ----------

        templateName: str
                The name of the template.

        to_xarray: bool
                Convert the simulations to xarray
                default: True

        to_database: bool
                Save simulation results to the database
                default: False

        forceKeep: bool
                If to_xarray is true,
                Determine wehter to keep the original files.

                If False, removes the Lagrangian files.
                defaultL False

        Returns
        -------
            The template by the name
        """
        doc = self.getDatasourceDocument(datasourceName=templateName,version=templateVersion)
        return LSMTemplate(doc,self)

    def listTemplates(self, wideFormat=False, **query):
        """
            :ist the template parameters that fulfil the query
        :param query:
        :return:
        """
        docList = self.getDatasourceDocumentsList(**query)
        if len(docList) > 0:
            descList = [doc.desc.copy() for doc in docList]
            for (i, desc) in enumerate(descList):
                desc.update({'id':docList[i].id})
                desc.update({'projectName': docList[i].projectName})

            params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
            params_df_list = [df.rename(columns=dict([(x,"params__%s"%x) for x in df.columns])) for df in params_df_list]
            desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
            df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]
            ret = pandas.concat(df_list,ignore_index=True,sort=False)
            if not wideFormat:
                ret = ret.melt()
        else:
            ret = []
        return ret


    def loadData(self,fileNameOrData,saveMode=toolkit.TOOLKIT_SAVEMODE_FILEANDDB,**kwargs):
        """
            Load a template object. Possibly saves to the DB.

        Parameters
        ----------
        fileNameOrData: str or JSON str or a JSON object.

                If str , a file or an JSON str that describes a template.

        saveMode: str
                Can be either:

                    - TOOLKIT_SAVEMODE_NOSAVE   : Just create template object


                    - TOOLKIT_SAVEMODE_FILEANDDB : Creates the template and store to the DB as a source.
                                                    Raise exception if the entry exists.

                    - TOOLKIT_SAVEMODE_FILEANDDB_REPLACE: Creates the template and store to the DB as a source.
                                                    Replace the entry in the DB if it exists.

        Returns
        -------
            template.LSMTemplate object

            Return the template object.
        """

        if isinstance(fileNameOrData,str):
            if os.path.isfile(fileNameOrData):
                with open(fileNameOrData) as templateFile:
                    desc = json.load(templateFile)
            else:
                raise ValueError("fileNameOrData must be a JSON template file")
        else:
            raise ValueError("fileNameOrData must be a JSON template file")

        templateName = desc['name']
        version = kwargs.get("version",None)

        doc = None
        if saveMode in [toolkit.TOOLKIT_SAVEMODE_FILEANDDB,toolkit.TOOLKIT_SAVEMODE_FILEANDDB_REPLACE]:

            doc = self.getDatasourceDocument(templateName)

            if doc is None:

                self.addDataSource(dataSourceName=templateName,
                                   resource=fileNameOrData,
                                   dataFormat=datalayer.datatypes.STRING,
                                   version=version,
                                   **kwargs)

            else:
                if  (saveMode == toolkit.TOOLKIT_SAVEMODE_FILEANDDB):
                    raise ValueError(f"Template {templateName} already exists in the DB")

                doc.resource = fileNameOrData
                doc.desc = kwargs
                doc.desc['version'] = version

        if doc is None:
            doc = datalayer.nonDBMetadataFrame(data=fileNameOrData,**kwargs)

        return LSMTemplate(doc,self)

    def getSimulations(self, **query):
        """
        get a list of SingleSimulation objects that fulfill the query
        :param query:
        :return:
        """
        for keys, unit in zip([minuteKeys, meterKeys, secondKeys, velocityKeys], [min, m, s, m / s]):
            for key in keys:
                if key in query.keys():
                    query[key] = query[key].asNumber(unit)
        docList = self.getSimulationsDocuments(type=LSMTemplate("",self).doctype_simulation, **query)
        retList = []
        for doc in docList:
            try:
                retList.append(SingleSimulation(doc))
            except:
                print(f"Warning: could not find data with the following document: {doc.asDict()}")
        return

    def listSimulations(self, wideFormat=False, **query):
        """
            List the Simulation parameters that fulfil the query
        :param query:
        :return:
        """
        docList = self.getSimulationsDocuments(type=LSMTemplate("",self).doctype_simulation, **query)
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
