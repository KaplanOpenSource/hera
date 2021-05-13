import os
import glob
import struct
import pandas
import dask.dataframe as dd


class Parser_CampbellBinary(object):
    _lut = None
    _dataContent = None # The array of the data.
    _basenum = None

    byteSize = None # the size of one record in bytes.

    _chunkSize = None # The number of records to read in each batch.

    def __init__(self, chunkSize = 10000):
        self._lut = {}
        self._dataContent = 0
        self._basenum = 0
        self._chunkSize = chunkSize

    def parse(self, path, fromTime=None, toTime=None, **metadata):
        if os.path.isfile(path):
            df = self.getPandasFromFile(path, fromTime=fromTime, toTime=toTime)
        else:
            df = self.getPandasFromDir(path, fromTime=fromTime, toTime=toTime)

        metadata_dict = dict()
        stations = df['station'].unique()
        for station in stations:
            station_metadata = metadata_dict.setdefault(station, dict())
            station_df = df.query("station==@station")
            instruments = station_df['instrument'].unique()
            for instrument in instruments:
                instrument_df = station_df.query("instrument==@instrument")
                heights = list(instrument_df['height'].unique())
                instrument_metadata = station_metadata.setdefault(instrument, dict())
                for height in heights:
                    metadata.update(dict(station=station, instrument=instrument, height=int(height)))
                    instrument_metadata.setdefault(int(height), metadata.copy())

        loaded_dask = dd.from_pandas(df, npartitions=1)
        return loaded_dask, metadata_dict

    def getPandasFromFile(self, path, fromTime, toTime):
        cbi = CampbellBinaryInterface(file=path)
        ts, cols, data = self.getData(path, fromTime=fromTime, toTime=toTime)
        dfList = []
        for i, key in enumerate(data.keys()):
            columns = cols[i]
            tmp_df = pandas.DataFrame(data[key], index=ts, columns=columns)
            tmp_df['height'] = int(key)
            tmp_df['station'] = cbi.headers[0].split(',')[1]
            tmp_df['instrument'] = cbi.headers[0].split(',')[-1]
            dfList.append(tmp_df)
        return pandas.concat(dfList, sort=True)

    def getPandasFromDir(self, path, fromTime, toTime):
        dfList = []
        for file in glob.glob(os.path.join(path, '*.dat')):
            tmp_df = self.getPandasFromFile(file, fromTime=fromTime, toTime=toTime)
            dfList.append(tmp_df)
        return pandas.concat(dfList, sort=True)

    def getData(self, file, fromTime, toTime):
        cbi = CampbellBinaryInterface(file=file)
        retVal = {}

        for i in cbi.heights:
            retVal[i] = []

        if type(fromTime) == str:
            fromTime = pandas.Timestamp(fromTime)

        if type(toTime) == str:
            toTime = pandas.Timestamp(toTime)

        recordIndex = 0 if fromTime is None else cbi.getRecordIndexByTime(fromTime)
        endIndex = cbi.recordsNum if toTime is None else cbi.getRecordIndexByTime(toTime)+1

        ts = []

        while recordIndex < endIndex:
            time, line = cbi.getRecordByIndex(recordIndex)
            ts.append(time)

            for i, key in enumerate(retVal):
                retVal[key].append(line[cbi.columnsIndexes[i][0]: cbi.columnsIndexes[i][1]])
            recordIndex += 1

        return ts, cbi.columnsNames, retVal

class Parser_TOA5(object):

    def __init__(self):
        pass

    def parse(self, file):
        pass


#############################################################



class CampbellBinaryInterface(object):
    _file = None
    _binData = None
    _headersSize = None
    _headers = None
    _recordSize = None
    _format = None
    _rawFormat = None
    _lut = None
    _firstTime = None
    _lastTime = None
    _columnsNames = None
    _columnsIndexes = None

    @property
    def headersSize(self):
        if self._headersSize is None:
            self._headersSize = self._getHeadersSize()
        return self._headersSize

    @property
    def headers(self):
        if self._headers is None:
            self._headers = self._getHeaders()
        return self._headers

    @property
    def station(self):
        return self.headers[0].split(',')[1]

    @property
    def instrument(self):
        return self.headers[0].split(',')[-1]

    @property
    def heights(self):
        if len(self.columnsNames) == 1:
            return [10]
        else:
            return [6, 11, 16]

    @property
    def recordSize(self):
        if self._recordSize is None:
            if self.headers[4].find(",") == -1:
                raise Exception("Missing Format Descriptor in line 4....")
            self._recordSize = struct.calcsize(self.format)

        return self._recordSize

    @property
    def recordsNum(self):
        return (len(self._binData)-self.headersSize)//self.recordSize

    @property
    def rawFormat(self):
        if self._rawFormat is None:
            self._getFormat()
        return self._rawFormat

    @property
    def format(self):
        if self._format is None:
            self._format = self._getFormat()
        return self._format

    @property
    def firstTime(self):
        if self._firstTime is None:
            self._firstTime = self._getFirstTime()
        return self._firstTime

    @property
    def lastTime(self):
        if self._lastTime is None:
            self._lastTime = self._getLastTime()
        return self._lastTime

    @property
    def columnsNames(self):
        if self._columnsNames is None:
            self._columnsNames = self._getColumnNames()
        return self._columnsNames

    @property
    def columnsIndexes(self):
        if self._columnsIndexes is None:
            self._columnsIndexes = self._getColumnIndexes()
        return self._columnsIndexes

    @property
    def binData(self):
        return self._binData

    def __init__(self, file):
        self._file = file

        with open(file, 'rb') as binFile:
            self._binData = binFile.read()

        self._lut = {}

    def _getFormat(self):
        rawFormata = self.headers[4].split(",")
        self._rawFormat = len(rawFormata) * ['']

        format = "<"
        for i in range(len(rawFormata)):
            self._rawFormat[i] = rawFormata[i].strip('"')
            if rawFormata[i] == 'ULONG':
                format += "I"
            elif rawFormata[i] == 'FP2':
                format += "H"
            elif rawFormata[i] == 'IEEE4':
                format += "f"
            elif rawFormata[i] == 'IEEE8':
                format += "d"
            elif rawFormata[i] == 'USHORT':
                format += "H"
            elif rawFormata[i] == 'LONG':
                format += "l"
            elif rawFormata[i] == 'BOOL':
                format += "?"
            elif rawFormata[i].find("ASCII(") != -1:
                format += rawFormata[i][6: -1] + 's'
            else:
                raise Exception("Unknown {} Format....".format(rawFormata[i]))
        return format

    def _getHeaders(self):
        headers = []
        numlf = 5
        basenum = 0
        tempstr = ''

        while numlf > 0:
            tempstr += chr(self._binData[basenum])
            if self._binData[basenum] == 10:
                headers.append(tempstr)
                tempstr = ''
                numlf -= 1
            basenum += 1

        for i in range(len(headers)):
            headers[i] = headers[i].replace('"', '').replace('\r\n', '')

        # headers[0] = headers[0].replace('TOB1', 'TOA5')
        # headers[1] = headers[1].replace('SECONDS,NANOSECONDS', 'TIMESTAMP')
        # headers[2] = headers[2].replace('SECONDS,NANOSECONDS', 'TS')
        # headers[3] = headers[3][1:]
        # headers = headers[:-1]

        return headers

    def _getColumnNames(self):
        colheader = self.headers[1].upper()
        cols = []

        if colheader.find("U_") != -1:
            # Raw Sonic Binary data file
            for i in range(3):
                if colheader.find("U_{}".format(i + 1)) != -1:
                    cols.append(['u', 'v', 'w', 'T'])

        elif colheader.find("TC_T") != -1:
            if colheader.find("TC_T1") != -1:
                cols.append(['TcT'])
            else:
                for i in range(3):
                    if colheader.find("TC_T({})".format(i + 1)) != -1:
                        cols.append(['TcT'])

            cols[len(cols) - 1].append('TRH')
            cols[len(cols) - 1].append('RH')
        return cols

    def _getColumnIndexes(self):
        colheader = self.headers[1].upper()
        Indexes = []

        if colheader.find("U_") != -1:
            # Raw Sonic Binary data file
            for i in range(3):
                if colheader.find("U_{}".format(i + 1)) != -1:
                    Indexes.append([1 + 4 * i, 5 + 4 * i])


        elif colheader.find("TC_T") != -1:
            if colheader.find("TC_T1") != -1:
                Indexes.append([1, 2])
            else:
                for i in range(3):
                    if colheader.find("TC_T({})".format(i + 1)) != -1:
                        Indexes.append([i + 1, i + 2])

            Indexes[len(Indexes) - 1][1] += 2
        return Indexes

    def _getHeadersSize(self):
        numlf = 5
        header_size = 0

        while numlf > 0:
            if self._binData[header_size] == 10:
                numlf -= 1
            header_size += 1

        return header_size

    def _getFirstTime(self):
        time, _ = self.getRecordByIndex(0)
        return time

    def _getLastTime(self):
        time, _ = self.getRecordByIndex(self.recordsNum-1)
        return time

    def getRecordByIndex(self, i):
        index = self.headersSize+i*self.recordSize
        lastSec, lastmili, line = self._getDataFromStream(self._binData[index: index+self.recordSize])
        time = pandas.Timestamp(1990, 1, 1) + pandas.Timedelta(days=lastSec / 86400.0, milliseconds=lastmili)
        return time, line

    def getRecordByTime(self, time):
        i = self.getRecordIndexByTime(time)
        return self.getRecordByIndex(i)

    def _getDataFromStream(self, partStream):
        retval = list(struct.unpack(self.format, partStream))
        for i in range(3, len(retval)):
            if self.rawFormat[i] == 'FP2':
                retval[i] = self._newfloatConvert(retval[i])
            elif self.rawFormat[i].find("ASCII(") != -1:
                retval[i] = self._byteToStr(retval[i])
        return retval[0], retval[1] / 1000000, retval[2:]

    def _byteToStr(self,inpbyte):
        retval = ''
        for i in range(len(inpbyte)):
            retval += chr(inpbyte[i])
        return retval.strip('\0')

    def _floatConvert(self, hbyte, lowbyte):
        if (hbyte & 0x80) > 0:
            sign = -1.0
        else:
            sign = 1.0

        shorti = hbyte & 0x60
        if shorti == 0x60:
            factor = 1000.0
        elif shorti == 0x40:
            factor = 100.0
        elif shorti == 0x20:
            factor = 10.0
        else:
            factor = 1.0

        val = sign * ((hbyte & 0x1f) * 256.0 + lowbyte) / factor
        return val

    def _newfloatConvert(self, key):
        try:
            return self._lut[key]
        except:
            if key == 65183:
                self._lut[key] = float('nan')
                return
            val = self._floatConvert(int(key % 256), key / 256)
            self._lut[key] = val
            return val

    def getTimeByRecordIndex(self, i):
        time, _ = self.getRecordByIndex(i)
        return time

    def getRecordIndexByTime(self, time):
        upperIndex = self.recordsNum
        lowerIndex = 0
        recordTime = self.getTimeByRecordIndex((lowerIndex+upperIndex)//2)
        tmpUpperIndex = -1
        tmpLowerIndex = -1

        while recordTime != time and (tmpLowerIndex != lowerIndex or tmpUpperIndex != upperIndex):
            tmpUpperIndex = upperIndex
            tmpLowerIndex = lowerIndex
            if time > recordTime:
                lowerIndex = (lowerIndex+upperIndex)//2
            else:
                upperIndex = (lowerIndex+upperIndex)//2
            recordTime = self.getTimeByRecordIndex((lowerIndex+upperIndex)//2)

        if recordTime!=time:
            raise IndexError("There is no record at %s" % time)
        else:
            return (lowerIndex+upperIndex)//2


