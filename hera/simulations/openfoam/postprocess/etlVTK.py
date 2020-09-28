import pydoc
import glob
import os
import sys

version = sys.version_info[0]
if version == 3:
    from hera import datalayer
    import dask.dataframe as dd


def recurseNode(Tree,nodeName,nodeData,metadata,pipelines,path,name,projectName,JSONName):

    """
    This function is recursively passed on the pipeline tree from json file and execute the needed functions by file write format
    Parameter
    ----------

    Tree : a list of filters applied in the current tree before the current node
    nodeName : the current filter node name from the pipeline
    nodeData : the current filter node properties from the pipeline
    metadata : the metadata from the json file
    pipelines : the pipelines from the json file
    path : path to the main directory
    name : output folder/ hdf files string
    projectName : projectName string


    """
    funcName = "load_%s" % nodeData.get('write', 'None')
    mod=pydoc.locate(__name__)
    loader = getattr(mod,funcName)
    loader(Tree,nodeName,nodeData,metadata,pipelines,path,name,projectName,JSONName)
    ds = nodeData.get("downstream", {})
    for nodeName, nodeData in ds.items():
        Tree.append(nodeName)
        recurseNode(Tree, nodeName, nodeData,metadata,pipelines,path,name,projectName)

def load_None(Tree, filterName, filterPipe, metadata, pipelines, path, name, projectName):
    """
    This function is for unsaved filters. It does nothing but allow proper operation of the recursion

    """

    print("%s write label is not exist in pipeline. data will not be cached" % filterName)

def load_netcdf(Tree, filterName, filterPipe, metadata, pipelines, path, name, projectName,JSONName):

    """
    This function handles the saving of netcdf files to the database.
    all parameters are as described in recurseNode function

    parameters:
    -----------

    Tree:
    filterName:
    filterPipe:
    metadata:
    pipelines:
    path:
    name:
    projectName:

    """


    piplelineTree=str(Tree).replace("[", "").replace("]", "").replace("'", "").replace(", ", ".")

    filterProps=dict(filterPipe)
    if 'downstream' in filterProps.keys():
        del (filterProps['downstream'])


    docList=datalayer.Simulations.getDocuments(type='HermesOpenFoam',
                                               projectName=projectName,
                                               resource=path)

    OF_WorkFlow = docList[0].desc['OF_Workflow']
    
    datalayer.Simulations.addDocument(projectName=projectName,
                                      desc=dict(filter=filterName,
                                                pipeline=pipelines,
                                                pipline_metadata=metadata,
                                                filterpipeline=piplelineTree,
                                                OF_WorkFlow=OF_WorkFlow
                                                ),
                                      type="OFsimulation",
                                      resource=glob.glob(os.path.join(path, name,"netcdf","%s_%s*.nc" %(JSONName, filterName))),
                                      dataFormat="netcdf_xarray")

    print("%s_%s data cached as netcdf" %(JSONName, filterName))

def load_hdf(Tree, filterName, filterPipe, metadata, pipelines, path, name, projectName):

    """
    This function handles the saving of hdf files to the database.
    it converts them inro .parquet format
    all parameters are as described in recurseNode function

    parameters:
    -----------

    Tree:
    filterName:
    filterPipe:
    metadata:
    pipelines:
    path:
    name:
    projectName:

    """

    piplelineTree = str(Tree).replace("[", "").replace("]", "").replace("'", "").replace(", ", ".")

    hdfFiles=glob.glob("%s/%s/hdf/%s_*.hdf" % (path, name, name))
    data = dd.read_hdf(hdfFiles, key=filterName)
    data = data.reset_index()\
        .set_index("time")\
        .repartition(partition_size="100Mb")

    data.to_parquet("%s/%s/parquet/%s.parquet" % (path, name, filterName))

    # Adding a link to the parquet to the database

    filterProps=dict(filterPipe)
    if 'downstream' in filterProps.keys():
        del (filterProps['downstream'])


    datalayer.Measurements.addDocument(projectName=projectName,
                                       desc=dict(filter=filterName,
                                                 filterProperties=filterProps,
                                                 filterpipeline=piplelineTree,
                                                 pipeline=pipelines,
                                                 **metadata
                                                 ),
                                       type="OFsimulation",
                                       resource="%s/%s/parquet/%s.parquet" % (path, name, filterName),
                                       dataFormat="parquet")

    print("%s data cached as parquet" % filterName)

