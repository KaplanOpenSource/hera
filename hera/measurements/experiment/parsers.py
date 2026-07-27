import os
import glob
import struct
import pandas
import json
import dask.dataframe as dd


class Parser_OldStyleMetaDataParquet():
    """Parser for old-style metadata stored alongside Parquet data files."""

    def __init__(self):
        """Initialize the parser."""

        pass


    def parse(self, pathToData):

        """

        This format has 2 meta data files:

            * metadata.json - holds the information on all the stations and devices in the campaign.
            * campaignDescription.json - holds the path to the parquet files of the campaign.


        This function changes the format of the JSON to the new format
        that is easier to load to the Database.

        Parameters
        -----------
        pathToData: str
                The path to the directory that holds the meta data files.

        Returns
        -------

            A map of:

                Stations: { the metadata for the stations } (JSON, or pandas)
                Devices: [ A list of the metadata of each device, remember to add the folder name of the parquet  (the resource). ].

        """

        MDpath = os.path.join(pathToData,'metadata.json')
        with open(MDpath) as f:
            metadata = json.load(f)

        Descriptionpath = os.path.join(pathToData, 'campaignDescription.json')
        with open(Descriptionpath) as f:
            descriptionData = json.load(f)


        metadataNew= self.getExperimentDict(metadata,descriptionData)

        Experiment=dict()

        Experiment[descriptionData['campaignName']]=dict(Stations=descriptionData['stationLocations'],
                                                     devices=metadataNew['devices'],
                                                     trials=metadataNew['trials'],
                                                     assets=metadataNew['assets'],
                                                    toolkitHandler=descriptionData['toolkitHandler'])


        return Experiment

    def getExperimentDict(self,metadata,descriptionData):
        """
            creates a dictionary from trials and devices data.

        parameters
        ----------

        Returns
        -------
        expDict
            the resulted dictionary
        """

        devicesList,trialsList, assetsDict = self._getLists(metadata,descriptionData)

        # expDict=dict()
        expDict=dict(trials  = trialsList,
                     devices = devicesList,
                     assets = assetsDict)

        return expDict


    def _getLists(self,metadata,descriptionData):

        """
            creates a devices and trials list from metadata.

        parameters
        ----------

        Returns
        -------
        devicesList
            the devices list data
        trialsList
            the trials list data
        """

        trialsList = []
        devicesList = []
        assetsDict=dict()

        stationsData = pandas.DataFrame.from_dict(descriptionData['stationLocations'])\
            .reset_index()\
            .rename(columns={"index": "station_name"})

        trialName = 'measurement1'
        trialSet  = 'measurement'

        devicePropertiesDict = dict()
        trialPropertiesDict = dict()

        inst_counters = {}

        for stn in metadata['stations'].keys():
            stnmd = metadata['stations'][stn]
            stndata = stationsData.query("station_name==@stn")

            entityList = []
            typeList = []

            for i,inst in enumerate(stnmd['instruments'].keys()):

                for devHgt in stnmd['instruments'][inst]:

                    inst_counters[inst] = inst_counters.get(inst, 0) + 1
                    counter = inst_counters[inst]
                    deviceName= f'{inst}_{counter}'
                    deviceDataPath=os.path.join(descriptionData['pathToData'],stn,inst,devHgt)

                    entityList.append(deviceName)
                    typeList.append(inst)


                    deviceDesignDict = dict(name=deviceName,
                                            deviceLocation=stn,
                                            deviceCoords=[stndata['lat'].item(), stndata['lon'].item()],
                                            deviceHeight=devHgt)


                    devicesList.append(dict(deviceName=deviceName,
                                            deviceType=inst,
                                            properties=devicePropertiesDict,
                                            trials={trialName : dict(design=deviceDesignDict,
                                                                     deploy=deviceDesignDict)
                                                   },
                                            deviceDataPath=deviceDataPath
                                            )
                                       )

            assetTrialDict = dict()
            assetTrialDict[trialName] = dict(entities=entityList,
                                                 types=typeList)

            assetsDict[stn]=dict(name=stn,
                                 codeName=stndata['station_code'].item(),
                                 properties=dict(coords=[stndata['lat'].item(), stndata['lon'].item()],
                                                 avgAreaBuildingsHgt=stnmd['averagedHeight'],
                                                 buildingHgt=stnmd['buildingHeight']),
                                 trials=assetTrialDict)


        trialsList.append(dict(trialSet=trialSet,
                               trialName=trialName,
                               properties = dict(
                                   trialStnProps=trialPropertiesDict)
                               )
                          )

        return devicesList, trialsList, assetsDict

class Parser_CampbellBinary(object):
    """Parser for Campbell Scientific binary data files."""

    _lut = None
    _dataContent = None # The array of the data.
    _basenum = None

    byteSize = None # the size of one record in bytes.

    _chunkSize = None # The number of records to read in each batch.

    def __init__(self, chunkSize = 10000):
        """
        Initialize the Campbell binary parser.

        Parameters
        ----------
        chunkSize : int, optional
            Number of records to read per batch. Default is 10000.
        """
        self._lut = {}
        self._dataContent = 0
        self._basenum = 0
        self._chunkSize = chunkSize

    def parse(self, path, fromTime=None, toTime=None, **metadata):
        """
        Parse binary data from a file or directory.

        Parameters
        ----------
        path : str
            Path to a binary file or directory of binary files.
        fromTime : str or pandas.Timestamp, optional
            Start time for filtering records.
        toTime : str or pandas.Timestamp, optional
            End time for filtering records.
        **metadata
            Additional metadata to attach to each record.

        Returns
        -------
        loaded_dask : dask.dataframe.DataFrame
            Parsed data as a Dask DataFrame.
        metadata_dict : dict
            Nested dict of station/instrument/height metadata.
        """
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
        """
        Parse a single binary file into a pandas DataFrame.

        Parameters
        ----------
        path : str
            Path to the binary file.
        fromTime : str or pandas.Timestamp or None
            Start time filter.
        toTime : str or pandas.Timestamp or None
            End time filter.

        Returns
        -------
        pandas.DataFrame
            Concatenated DataFrame with height, station, and instrument columns.
        """
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
        """
        Parse all .dat files in a directory into a single DataFrame.

        Parameters
        ----------
        path : str
            Directory containing .dat binary files.
        fromTime : str or pandas.Timestamp or None
            Start time filter.
        toTime : str or pandas.Timestamp or None
            End time filter.

        Returns
        -------
        pandas.DataFrame
            Concatenated DataFrame from all files in the directory.
        """
        dfList = []
        for file in glob.glob(os.path.join(path, '*.dat')):
            tmp_df = self.getPandasFromFile(file, fromTime=fromTime, toTime=toTime)
            dfList.append(tmp_df)
        return pandas.concat(dfList, sort=True)

    def getData(self, file, fromTime, toTime):
        """
        Extract timestamps, column names, and data from a binary file.

        Parameters
        ----------
        file : str
            Path to the binary file.
        fromTime : str or pandas.Timestamp or None
            Start time filter.
        toTime : str or pandas.Timestamp or None
            End time filter.

        Returns
        -------
        ts : list
            List of pandas.Timestamp for each record.
        columnsNames : list
            Column name lists per height group.
        retVal : dict
            Data values keyed by height.
        """
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

############################## Private ###############################
class CampbellBinaryInterface(object):
    """Low-level interface for reading Campbell Scientific TOB1 binary files."""

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
        """int : Size of the header section in bytes."""
        if self._headersSize is None:
            self._headersSize = self._getHeadersSize()
        return self._headersSize

    @property
    def headers(self):
        """list of str : Parsed header lines from the binary file."""
        if self._headers is None:
            self._headers = self._getHeaders()
        return self._headers

    @property
    def station(self):
        """str : Station name extracted from the first header line."""
        return self.headers[0].split(',')[1]

    @property
    def instrument(self):
        """str : Instrument name extracted from the first header line."""
        return self.headers[0].split(',')[-1]

    @property
    def heights(self):
        """list of int : Measurement heights derived from column layout."""
        if len(self.columnsNames) == 1:
            return [10]
        else:
            return [6, 11, 16]

    @property
    def recordSize(self):
        """int : Size of a single data record in bytes."""
        if self._recordSize is None:
            if self.headers[4].find(",") == -1:
                raise Exception("Missing Format Descriptor in line 4....")
            self._recordSize = struct.calcsize(self.format)

        return self._recordSize

    @property
    def recordsNum(self):
        """int : Total number of data records in the file."""
        return (len(self._binData)-self.headersSize)//self.recordSize

    @property
    def rawFormat(self):
        """list of str : Raw Campbell format type strings (e.g., 'FP2', 'IEEE4')."""
        if self._rawFormat is None:
            self._getFormat()
        return self._rawFormat

    @property
    def format(self):
        """str : Struct format string for unpacking binary records."""
        if self._format is None:
            self._format = self._getFormat()
        return self._format

    @property
    def firstTime(self):
        """pandas.Timestamp : Timestamp of the first record."""
        if self._firstTime is None:
            self._firstTime = self._getFirstTime()
        return self._firstTime

    @property
    def lastTime(self):
        """pandas.Timestamp : Timestamp of the last record."""
        if self._lastTime is None:
            self._lastTime = self._getLastTime()
        return self._lastTime

    @property
    def columnsNames(self):
        """list of list of str : Column names grouped by height."""
        if self._columnsNames is None:
            self._columnsNames = self._getColumnNames()
        return self._columnsNames

    @property
    def columnsIndexes(self):
        """list of list of int : Start/end index pairs for each height group."""
        if self._columnsIndexes is None:
            self._columnsIndexes = self._getColumnIndexes()
        return self._columnsIndexes

    @property
    def binData(self):
        """bytes : Raw binary content of the file."""
        return self._binData

    def __init__(self, file):
        """
        Initialize the interface and read the binary file.

        Parameters
        ----------
        file : str
            Path to the Campbell binary file.
        """
        self._file = file

        with open(file, 'rb') as binFile:
            self._binData = binFile.read()

        self._lut = {}

    def _getFormat(self):
        """
        Build the struct format string from header format descriptors.

        Returns
        -------
        str
            Struct-compatible format string.
        """
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
        """
        Parse the five header lines from the binary data.

        Returns
        -------
        list of str
            Cleaned header lines.
        """
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
        """
        Derive column name lists from the second header line.

        Returns
        -------
        list of list of str
            Column names grouped by measurement height.
        """
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
        """
        Derive column index ranges from the second header line.

        Returns
        -------
        list of list of int
            Start/end index pairs for slicing record data per height.
        """
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
        """
        Calculate the byte offset where data records begin.

        Returns
        -------
        int
            Number of bytes occupied by the header section.
        """
        numlf = 5
        header_size = 0

        while numlf > 0:
            if self._binData[header_size] == 10:
                numlf -= 1
            header_size += 1

        return header_size

    def _getFirstTime(self):
        """
        Return the timestamp of the first record.

        Returns
        -------
        pandas.Timestamp
        """
        time, _ = self.getRecordByIndex(0)
        return time

    def _getLastTime(self):
        """
        Return the timestamp of the last record.

        Returns
        -------
        pandas.Timestamp
        """
        time, _ = self.getRecordByIndex(self.recordsNum-1)
        return time

    def getRecordByIndex(self, i):
        """
        Retrieve a single record by its zero-based index.

        Parameters
        ----------
        i : int
            Record index.

        Returns
        -------
        time : pandas.Timestamp
            Record timestamp.
        line : list
            Parsed data values for the record.
        """
        index = self.headersSize+i*self.recordSize
        lastSec, lastmili, line = self._getDataFromStream(self._binData[index: index+self.recordSize])
        time = pandas.Timestamp(1990, 1, 1) + pandas.Timedelta(days=lastSec / 86400.0, milliseconds=lastmili)
        return time, line

    def getRecordByTime(self, time):
        """
        Retrieve a record matching the given timestamp.

        Parameters
        ----------
        time : pandas.Timestamp
            Target timestamp.

        Returns
        -------
        time : pandas.Timestamp
            Record timestamp.
        line : list
            Parsed data values.
        """
        i = self.getRecordIndexByTime(time)
        return self.getRecordByIndex(i)

    def _getDataFromStream(self, partStream):
        """
        Unpack a binary record and convert special format types.

        Parameters
        ----------
        partStream : bytes
            Raw bytes for one record.

        Returns
        -------
        seconds : int
            Seconds since 1990-01-01.
        nanoseconds : float
            Sub-second nanosecond fraction in milliseconds.
        data : list
            Converted data values.
        """
        retval = list(struct.unpack(self.format, partStream))
        for i in range(3, len(retval)):
            if self.rawFormat[i] == 'FP2':
                retval[i] = self._newfloatConvert(retval[i])
            elif self.rawFormat[i].find("ASCII(") != -1:
                retval[i] = self._byteToStr(retval[i])
        return retval[0], retval[1] / 1000000, retval[2:]

    def _byteToStr(self,inpbyte):
        """
        Convert a byte sequence to a null-stripped string.

        Parameters
        ----------
        inpbyte : bytes
            Raw byte sequence.

        Returns
        -------
        str
            Decoded string with null characters removed.
        """
        retval = ''
        for i in range(len(inpbyte)):
            retval += chr(inpbyte[i])
        return retval.strip('\0')

    def _floatConvert(self, hbyte, lowbyte):
        """
        Convert Campbell FP2 two-byte value to a float.

        Parameters
        ----------
        hbyte : int
            High byte.
        lowbyte : int
            Low byte.

        Returns
        -------
        float
            Decoded floating-point value.
        """
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
        """
        Convert an FP2 value using a lookup-table cache.

        Parameters
        ----------
        key : int
            Raw FP2 integer value.

        Returns
        -------
        float or None
            Converted float, or NaN for the sentinel value 65183.
        """
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
        """
        Return the timestamp for a given record index.

        Parameters
        ----------
        i : int
            Record index.

        Returns
        -------
        pandas.Timestamp
        """
        time, _ = self.getRecordByIndex(i)
        return time

    def getRecordIndexByTime(self, time):
        """
        Find the record index for a given timestamp using binary search.

        Parameters
        ----------
        time : pandas.Timestamp
            Target timestamp.

        Returns
        -------
        int
            Index of the matching record.

        Raises
        ------
        IndexError
            If no record matches the given time.
        """
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


