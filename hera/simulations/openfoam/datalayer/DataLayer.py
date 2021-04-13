from .... import datalayer
import xarray
import os
import dask

class openfoam_Datalayer(object):


    def getDocFromFile(self,path,dataFormat='xarray',**kwargs):
        if dataFormat=='xarray':
            dataFormat='netcdf_xarray'
            data= xarray.open_mfdataset(os.path.join(path,'*.nc'), combine='by_coords')

        elif dataFormat =='pandas':
            dataFormat = 'parquet'
            data = dask.dataframe.read_parquet(path)
            pass
        else:
            print("data format must be one of the follows: [xarray ; pandas]")

        return [datalayer.document.metadataDocument.nonDBMetadataFrame(data=data, resource=path, dataFormat=dataFormat, type='OFsimulation', **kwargs)]



    def getDocFromDB(self):

        pass

    def getNodeList(self,doc):
        try:
            nodeList=doc.desc['OF_workflow']['workflow']['nodeList']
            return list(nodeList)
        except:
            print(f"doc is from type {doc.type} and does not have 'OF_workflow' key")