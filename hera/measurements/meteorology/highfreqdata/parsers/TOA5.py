import csv

class Parser:
    def __init__(self):
        pass

    def parse(self, path :str, fromTime=None, toTime=None):                         ## Will return N parquet files as number of devices inside a binary file
        if path[len(path)-4:] == ".dat":                                            ## File path includes .dat at the end
            dfs = self.getPandasFromFile(self, path, fromTime, toTime)
        else:
            dfs = self.getPandasFromDir(self, path, fromTime, toTime)               ## Folder Path if it does not include .dat at the  end
        return dfs


    def getPandasFromFile(self, path, fromTime, toTime):
        # abi = ASCIIBinaryInterface(file=path)
        ts, cols, data = self.getData(path, fromTime=fromTime, toTime=toTime)
        dfList = []
        for i, key in enumerate(data.keys()):
            columns = cols[i]
            tmp_df = pandas.DataFrame(data[key], index=ts, columns=columns)
            # tmp_df['height'] = int(key)
            # tmp_df['station'] = cbi.headers[0].split(',')[1]
            # tmp_df['instrument'] = cbi.headers[0].split(',')[-1]
            dfList.append(tmp_df)
        return dfList

    def getPandasFromDir(self, path, fromTime, toTime):
        pass

    def getData(self, file, fromTime, toTime):
        abi = ASCIIBinaryInterface(file=file)
        retVal = {}

        # for i in cbi.heights:
        #     retVal[i] = []


        # recordIndex = 0 if fromTime is None else cbi.getRecordIndexByTime(fromTime)
        # endIndex = cbi.recordsNum if toTime is None else cbi.getRecordIndexByTime(toTime) + 1

        recordIndex = 0
        endIndex = abi.recordsNum

        ts = []

        while recordIndex < endIndex:
            time, line = abi.getRecordByIndex(recordIndex)
            ts.append(time)

            # for i, key in enumerate(retVal):





class ASCIIBinaryInterface(object):
    _file = None
    _data = None
    _headers = None
    _columnsNames = None
    _recordsNum = None


    def __init__(self, file):
        self._file = file
        with open(file) as fp:
            self._data = list(csv.reader(fp))

    @property
    def headers(self):
        if self._headers is None:
            self._headers = self._getHeaders()
        return self._headers

    @property
    def columnsNames(self):
        if self._columnsNames is None:
            self._columnsNames = self._getColumnNames()
        return self._columnsNames

    @property
    def recordsNum(self):
        return len(self._data[4:])


    def _getHeaders(self):
        headers = []
        for i in range(4):
            headers.append(self._data[i])
        return headers

    def _getColumnNames(self):
        colheader = str(self.headers[1]).upper()
        cols = []

        if colheader.find("U_") != -1:                                     #U,V,W,T
            # Raw Sonic Binary data file
            for i in range(3):
                if colheader.find("U_{}".format(i + 1)) != -1:
                    cols.append(['u', 'v', 'w', 'T'])

        elif colheader.find("TC_T") != -1:                                 #TCT
            if colheader.find("TC_T1") != -1:
                cols.append(['TcT'])
            else:
                for i in range(3):
                    if colheader.find("TC_T({})".format(i + 1)) != -1:
                        cols.append(['TcT'])

            cols[len(cols) - 1].append('TRH')
            cols[len(cols) - 1].append('RH')
        return cols


    def getRecordByIndex(self, i):
        i = i+4
        return self._data[i][0] , self._data[i][1:]

