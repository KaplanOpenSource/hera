from ...datalayer import Simulations
from .DataLayer import SingleSimulation
from ..templates import LSMTemplate
import pandas
from itertools import product

DEFAULTPROJET = "LSM"

DOCTYPE_TEMPLATE = "LSM_Template"
DOCTYPE_RUN      = "LSM_run"

def getTemplates(**query):
    """
    get a list of Template objects that fulfill the query
    :param query:
    :return:
    """

    docList = Simulations.getDocuments(projectName=DEFAULTPROJET , type=DOCTYPE_TEMPLATE , **query)
    return [LSMTemplate(doc) for doc in docList]

def getTemplateByID(id):
    """
    get a teamplate by document id

    :param id:
    :return:
    """
    return LSMTemplate(Simulations.getDocumentByID(id))

def listTemplates(projectName=DEFAULTPROJET , wideFormat=False, **query):
    """
    list the template parameters that fulfil the query
    :param query:
    :return:
    """
    docList = Simulations.getDocuments(projectName=projectName, type=DOCTYPE_TEMPLATE , **query)
    descList = [doc.desc.copy() for doc in docList]
    for (i, desc) in enumerate(descList):
        desc.update({'id':docList[i].id})
    params_df_list = [pandas.DataFrame(desc.pop('params'), index=[0]) for desc in descList]
    params_df_list = [df.rename(columns=dict([(x,"params__%s"%x) for x in df.columns])) for df in params_df_list]
    desc_df_list = [pandas.DataFrame(desc, index=[0]) for desc in descList]
    df_list = [desc.join(params) for (desc,params) in product(desc_df_list, params_df_list)]

    ret = pandas.concat(df_list,ignore_index=True,sort=False)
    if wideFormat:
        ret = ret.melt()

    return ret

def getSimulations(projectName=DEFAULTPROJET , **query):
    """
    get a list of SingleSimulation objects that fulfill the query
    :param query:
    :return:
    """

    docList = Simulations.getDocuments(projectName=projectName, type=DOCTYPE_RUN, **query)
    return [SingleSimulation(doc) for doc in docList]

def getSimulationByID(id):
    """
    get a simulation by document id

    :param id:
    :return:
    """
    return SingleSimulation(Simulations.getDocumentByID(id))

def listSimulations(projectName=DEFAULTPROJET , wideFormat=False, **query):
    """
        List the Simulation parameters that fulfil the query
    :param query:
    :return:
    """
    docList = Simulations.getDocuments(projectName=projectName, type=DOCTYPE_RUN, **query)
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