import pandas
import json
import os

from .template import LSMTemplate
from itertools import product
from ... import datalayer
from ... import toolkit
from .singleSimulation import SingleSimulation
from unum.units import *
from ..utils.coordinateHandler import coordinateHandler


class LSMToolkit(toolkit.abstractToolkit):
    """
        The LSM.old toolkit

        The datasources are the templates.

        The LSM.old template JSON structure is :


        << TO DO >>



    """
    TRUE = ".TRUE."
    FALSE = ".FALSE."
    _to_xarray = None
    _to_database = None
    _forceKeep = None
    _analysis = None

    @property
    def analysis(self):
        return self._analysis
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

    def __init__(self, projectName, filesDirectory=None, to_xarray=True, to_database=False, forceKeep=False):
        """
        Initialize the LSM toolkit.

        Parameters
        ----------
        projectName : str
            The name of the project that contains the simulation templates and results.
            The project must exist in the Hera database.

        filesDirectory : str, optional
            The directory to save simulation output files. If None, uses the project's
            default files directory. Default is None.

        to_xarray : bool, optional
            If True, automatically convert simulation results to xarray format.
            This enables easier data analysis and visualization. Default is True.

        to_database : bool, optional
            If True, save simulation run metadata to the database. This allows
            querying and tracking of simulation runs. Default is False.

        forceKeep : bool, optional
            If True and to_xarray is True, keep the original Lagrangian particle
            files after conversion. If False, removes them to save disk space.
            Default is False.

        Examples
        --------
        >>> lsm_tk = LSMToolkit(
        ...     projectName="urban_dispersion",
        ...     filesDirectory="/path/to/simulations",
        ...     to_xarray=True,
        ...     to_database=True
        ... )
        """
        super().__init__(projectName=projectName, toolkitName="LSM.old", filesDirectory=filesDirectory)

        self.to_xarray = to_xarray
        self.to_database = to_database
        self.forceKeep = forceKeep
        self._analysis = analysis(self)

    def getTemplates(self, **query):
        """
        Get a list of template objects that match the query criteria.

        This method queries the datasource documents and returns LSMTemplate
        objects for each matching template.

        Parameters
        ----------
        **query : dict
            Query filters to apply when searching for templates. These are passed
            directly to `getDataSourceDocumentsList()`. Common filters include:
            - datasourceName: Filter by template name
            - version: Filter by template version
            - Any other fields stored in the template document description

        Returns
        -------
        list of LSMTemplate
            List of LSMTemplate objects matching the query criteria.

        Examples
        --------
        >>> # Get all templates
        >>> templates = lsm_tk.getTemplates()
        >>> 
        >>> # Get templates with specific name
        >>> templates = lsm_tk.getTemplates(datasourceName="urban_release")
        >>> 
        >>> # Get templates with specific version
        >>> templates = lsm_tk.getTemplates(version=(1, 0, 0))
        """
        docList = self.getDataSourceDocumentsList(**query)
        return [LSMTemplate(doc,self) for doc in docList]

    def getTemplateByName(self, templateName, templateVersion=None):
        """
        Retrieve a template by its name and optional version.

        This method looks up a template document from the datasources and returns
        an LSMTemplate object. If version is not specified, the default version
        (or latest version) is used.

        Parameters
        ----------
        templateName : str
            The name of the template to retrieve. This must match a datasource
            name registered in the project.

        templateVersion : tuple, optional
            Specific version of the template to retrieve (e.g., (1, 0, 0)).
            If None, uses the default version configured in the project, or the
            latest version if no default is set. Default is None.

        Returns
        -------
        LSMTemplate
            The template object corresponding to the requested name and version.

        Raises
        ------
        ValueError
            If the template is not found in the project.

        Examples
        --------
        >>> # Get default version of a template
        >>> template = lsm_tk.getTemplateByName("v4-general")
        >>> 
        >>> # Get specific version
        >>> template = lsm_tk.getTemplateByName("v4-general", templateVersion=(1, 2, 0))
        """
        doc = self.getDataSourceDocument(datasourceName=templateName, version=templateVersion)
        return LSMTemplate(doc,self)

    def getTemplatesTable(self, **query):
        """
        Get a pandas DataFrame listing template parameters that match the query.

        This method returns a normalized table combining template metadata and
        parameter values. Useful for comparing templates or finding templates
        with specific parameter values.

        Parameters
        ----------
        **query : dict
            Query filters to apply when searching for templates. Same as
            `getTemplates()`.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with columns for template metadata (id, projectName, etc.)
            and parameter values. Each row represents one template with its
            parameters flattened.

        Examples
        --------
        >>> # Get table of all templates
        >>> df = lsm_tk.getTemplatesTable()
        >>> 
        >>> # Filter templates by name
        >>> df = lsm_tk.getTemplatesTable(datasourceName="urban_release")
        """
        docList = self.getDataSourceDocumentsList(**query)
        if len(docList) > 0:
            descList = [doc.desc.copy() for doc in docList]
            for (i, desc) in enumerate(descList):
                desc.update({'id':docList[i].id})
                desc.update({'projectName': docList[i].projectName})

            params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
            desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
            df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]
            ret = pandas.concat(df_list,ignore_index=True,sort=False)
        else:
            ret = pandas.DataFrame()
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

            doc = self.getDataSourceDocument(templateName)

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

    def getSimulations(self, simulationName=None, unitsTemplateVersion="v4-general", **query):
        """
        Get a list of SingleSimulation objects matching the query criteria.

        This method retrieves simulation runs from the database. Parameters in the
        query are automatically converted to the correct units based on the template's
        unit definitions.

        Parameters
        ----------
        simulationName : str, optional
            Specific simulation name to retrieve. If provided, only simulations
            with this exact name are returned. Default is None.

        unitsTemplateVersion : str, optional
            Template version to use for unit conversion. The template's unit
            definitions are used to convert query parameters. Default is "v4-general".

        **query : dict
            Query parameters to filter simulations. These can include:
            - Any simulation parameter (e.g., releaseHeight, releaseLocation)
            - Parameters are automatically converted to internal units
            - Use unum units objects for unit-aware queries

        Returns
        -------
        list of SingleSimulation
            List of SingleSimulation objects matching the query criteria.

        Examples
        --------
        >>> from unum.units import *
        >>> 
        >>> # Get simulations with specific release height
        >>> sims = lsm_tk.getSimulations(releaseHeight=10*m)
        >>> 
        >>> # Get simulations by name
        >>> sims = lsm_tk.getSimulations(simulationName="run_001")
        >>> 
        >>> # Multiple query parameters
        >>> sims = lsm_tk.getSimulations(
        ...     releaseHeight=10*m,
        ...     releaseLocation=(1000, 2000),
        ...     unitsTemplateVersion="v4-general"
        ... )
        """
        template = self.getTemplateByName(unitsTemplateVersion)

        for key in template._document['desc']["units"].keys():
            if key in query.keys():
                query[key] = query[key].asNumber(eval(template._document['desc']["units"][key]))
        queryWithParams = {}
        for key in query.keys():
            queryWithParams[f"params__{key}"] = query[key]

        if simulationName is not None:
            queryWithParams["simulationName"] = simulationName

        docList = self.getSimulationsDocuments(type=LSMTemplate("",self).doctype_simulation, **queryWithParams)
        retList = []
        for doc in docList:
            try:
                retList.append(SingleSimulation(doc))
            except:
                print(f"Warning: could not find data with the following document: {doc.asDict()}")

        return retList

    def getSimulationsList(self, wideFormat=False, **query):
        """
        List simulation parameters in a tabular format.

        Returns a DataFrame containing simulation metadata and parameters. Can be
        returned in wide format (one column per parameter) or long format
        (parameter-value pairs).

        Parameters
        ----------
        wideFormat : bool, optional
            If True, returns a wide-format DataFrame with one column per parameter.
            If False, returns a long-format DataFrame with variable-value pairs.
            Default is False.

        **query : dict
            Query filters to apply when searching for simulations. Same as
            `getSimulations()`.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing simulation parameters. Format depends on
            wideFormat parameter.

        Raises
        ------
        FileNotFoundError
            If no simulations match the query criteria.

        Examples
        --------
        >>> # Get long format list
        >>> df = lsm_tk.getSimulationsList()
        >>> 
        >>> # Get wide format list
        >>> df = lsm_tk.getSimulationsList(wideFormat=True)
        >>> 
        >>> # Filter by query
        >>> df = lsm_tk.getSimulationsList(releaseHeight=10*m)
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
            raise FileNotFoundError('No simulations.old found')

class analysis:
    """
    Analysis layer for the LSM toolkit.

    Provides analysis and post-processing capabilities for LSM simulation results.
    This includes coordinate handling and spatial analysis utilities.

    Parameters
    ----------
    dataLayer : LSMToolkit
        The LSM toolkit instance to provide analysis for.

    Attributes
    ----------
    datalayer : LSMToolkit
        Reference to the toolkit instance.
    coordinateHandler : coordinateHandler
        Handler for coordinate system conversions and transformations.
    """

    _datalayer = None
    _coordinateHandler = None

    @property
    def coordinateHandler(self):
        return self._coordinateHandler

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer
        self._coordinateHandler = coordinateHandler()
