import csv
import pandas as pd
from tqdm import tqdm
import os


class ASCIIParser:
    def __init__(self):
        pass

    def parse(self, path :str, fromTime=None, toTime=None):                                                             ## Will return N parquet files as number of devices inside a binary file
        if path[len(path)-4:] == ".dat":                                                                                ## File path includes .dat at the end
            dfs = self.getPandasFromFile(path, fromTime, toTime)
        else:
            dfs = self.getPandasFromDir(path, fromTime, toTime)                                                         ## Folder Path if it does not include .dat at the  end
        return dfs

    def getPandasFromFile(self, path, fromTime, toTime):
        """
         Reads the ascii file from fromTime to toTime and returns dictionary where
         Key is Device ID and value is dataframe for the corresponding device.
         Parameters
         ---------
           path: str
           fromTime:  None, str - should include date and hour seperated with space (order does not matter)
           toTime: None, str - should include date and hour seperated with space (order does not matter)

           Returns
           -------
           Dictionary where Key is Device ID and value is dataframe for the corresponding device.
        """

        cols,number_of_devices,include_metadata = self.get_columns(path)                                                                 ## Read metadata for detecting Raw Sonic or TCT and Number of Devices
        dfs = {}
        allDevices = pd.read_csv(path)

        if include_metadata:
            allDevices = allDevices.rename(columns=allDevices.iloc[0])[3:].reset_index(drop=True)                           ## Rename false metadata columns to real ones and remove metadata (first 3 rows)
        else:
            allDevices.columns = ['TIMESTAMP','RECORD'] + [f"{x}_{i+1}" for i in range(number_of_devices) for x in cols]

        if fromTime or toTime:                                                                                          ##Handle from time and to time
            allDevices = self.fromTime_toTime_handler(allDevices,fromTime,toTime)

        for device_ID in range(1,number_of_devices+1):                                                                  ## iterate in all cols and devices

            if cols[0] =="U":
                columnsName = ['TIMESTAMP','RECORD'] + [f"{x}_{device_ID}" for x in cols]                               ##For raw sonic there is _
                columnNameMapping = dict([(f"{x}_{device_ID}", f"{x}") for x in cols])
                device = "Raw_Sonic"
            else:
                columnsName = ['TIMESTAMP', 'RECORD'] + [f"{x}{device_ID}" for x in cols]                               ## For TCT there is no _
                columnNameMapping = dict([(f"{x}{device_ID}", f"{x}") for x in cols])
                device = "TCT_TRH"

            columnNameMapping['TIMESTAMP'] = 'TIMESTAMP'
            columnNameMapping['RECORD'] = 'RECORD'

            dfs[f"{device}_{device_ID}"] = allDevices[columnsName].copy().rename(columns=columnNameMapping)

        return dfs

    def get_columns(self,path):
        """
         Reads the ascii file columns and returns columns names and number of devices
         ---------
           path: str

           Returns
           -------
           columns names and number of devices.
           - ['U','V','W','T'] for Raw Sonic
           - ['TC_T', 'TRH', 'RH'] for TCT

        """
        include_metadata = True
        with open(path) as meta_data:
            meta_data_reader = csv.reader(meta_data)
            device_type = next(meta_data_reader)
            if device_type[-1]=="Raw_Sonic":
                cols = ['U','V','W','T']
                number_of_devices = len([word for word in next(meta_data_reader) if word.startswith(cols[0])])
            elif device_type[-1]=="Tct_Trh":
                cols = ['TC_T', 'TRH', 'RH']
                number_of_devices = len([word for word in next(meta_data_reader) if word.startswith(cols[0])])
            else:                                       ##No metadata in file
                include_metadata = False
                if 'Raw_Sonic' in path:
                    cols = ['U', 'V', 'W', 'T']
                elif 'Tct_Trh' in path:
                    cols = ['TC_T', 'TRH', 'RH']
                else:
                    raise ValueError("No Device could be detected by given file. File should include metadata or device name in file name.")

                number_of_devices = int(len(device_type[2:])/4)


            return cols,number_of_devices,include_metadata

    def getPandasFromDir(self,path, fromTime, toTime):
        """
         Recives Directory of ASCII files and returns dictionary where Key is Device ID and value is dataframe for the corresponding device.
         ---------
           path: str
           fromTime:  None, str - should include date and hour seperated with space (order does not matter)
           toTime: None, str - should include date and hour seperated with space (order does not matter)

           Returns
           -------
           Dictionary where Key is Device ID and value is dataframe for the corresponding device.
        """

        allDevices_in_directory = {}
        for file_path in tqdm(os.listdir(path)):
            file_path = os.path.join(path,file_path)
            file_dict = self.getPandasFromFile(file_path,fromTime,toTime)
            for device in file_dict.keys():
                if device in allDevices_in_directory:
                    allDevices_in_directory[device] = pd.concat([allDevices_in_directory[device],file_dict[device]],axis=0).reset_index(drop=True)
                else:
                    allDevices_in_directory[device] = file_dict[device]

        for device in allDevices_in_directory.keys():
            allDevices_in_directory[device] = allDevices_in_directory[device].set_index("TIMESTAMP").sort_index().reset_index()

        return allDevices_in_directory

    def fromTime_toTime_handler(self,df,fromTime,toTime):
        """
         Recives dataframe and time intervals, returns dataframe within time interval .

         Parameters
         ---------
           df: pd.DataFrame
           fromTime: None, str - should include date and hour seperated with space (order does not matter)
           toTime: None, str - should include date and hour seperated with space (order does not matter)

           Returns
           -------
           return dataframe between fromTime until toTime
        """

        df = df.set_index('TIMESTAMP')
        start = None
        end = None

        #assert len(fromTime.split(" "))==2 and len(toTime.split(" "))==2, "Please enter date and time with space seperation"
        if fromTime and toTime:
            assert pd.to_datetime(fromTime) <= pd.to_datetime(toTime), f"fromTime {fromTime} is larger than toTime {toTime}"
            start = str(pd.to_datetime(fromTime))
            end = str(pd.to_datetime(toTime))
        elif fromTime and not toTime:
            start = str(pd.to_datetime(fromTime))
            end = str(pd.to_datetime(df.iloc[len(df)-1].name))
        else:
            start = str(pd.to_datetime(df.iloc[0].name))
            end = str(pd.to_datetime(toTime))


        return df[start:end].reset_index()


