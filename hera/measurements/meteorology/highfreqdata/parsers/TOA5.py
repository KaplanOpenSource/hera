import csv
import pandas as pd


class Parser:
    def __init__(self):
        pass

    def parse(self, path :str, fromTime=None, toTime=None):                             ## Will return N parquet files as number of devices inside a binary file
        if path[len(path)-4:] == ".dat":                                                ## File path includes .dat at the end
            dfs = self.getPandasFromFile(self, path, fromTime, toTime)
        else:
            dfs = self.getPandasFromDir(self, path, fromTime, toTime)                   ## Folder Path if it does not include .dat at the  end
        return dfs

    def getPandasFromFile(self, path, fromTime, toTime):
        cols,number_of_devices = self.get_columns(path)                                 ## Read metadata for detecting Raw Sonic or TCT and Number of Devices
        dfs = []
        allDevices = pd.read_csv(path)                                                  ## Read all csv
        allDevices = allDevices.rename(columns=df.iloc[0])[3:].reset_index(drop=True)   ## Rename false metadata columns to real ones and remove metadata (first 3 rows)
        for devince_ID in range(1,number_of_devices+1):                                 ## iterate in all cols and devices
            columnsName = ['TIMESTAMP','RECORD'] + [f"{x}_{devince_ID}" for x in cols]
            columnNameMapping = dict([(f"{x}_{i}", f"{x}") for x in cols])
            columnNameMapping['TIMESTAMP'] = 'TIMESTAMP'
            columnNameMapping['RECORD'] = 'RECORD'

            dfs.append(allDevices[columnsName].copy().rename(columns=columnNameMapping))

        return dfs

    def get_columns(self,path):
        with open(path) as meta_data:
            meta_data_reader = csv.reader(meta_data)
            device_type = next(meta_data_reader)[-1]
            if device_type=="Raw_Sonic":
                cols = ['U','V','W','T']
            else:
                cols = ['TC_T', 'TRH', 'RH']
            number_of_devices = len([word for word in next(meta_data_reader) if word.startswith(cols[0])])

            return cols,number_of_devices

    def getPandasFromDir(self,path, fromTime, toTime):
        pass