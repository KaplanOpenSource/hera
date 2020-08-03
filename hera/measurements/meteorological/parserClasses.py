import os
import glob
import struct
import pandas
import dask.dataframe as dd


class Parser_IMS(object):
    _removelist = None

    def __init__(self):
        self._removelist = ['BET DAGAN RAD', 'SEDE BOQER UNI', 'BEER SHEVA UNI']

    def _process_HebName(self, Station):
        HebName = Station.Stn_name_Heb.item()
        return HebName

    def _process_ITM_E(self, Station):
        ITM_E = Station.ITM_E.item()
        return ITM_E

    def _process_ITM_N(self, Station):
        ITM_N = Station.ITM_N.item()
        return ITM_N

    def _process_LAT_deg(self, Station):
        LAT_deg = float(Station.Lat_deg.item()[:-1])
        return LAT_deg

    def _process_LON_deg(self, Station):
        LON_deg = float(Station.Lon_deg.item()[:-1])
        return LON_deg

    def _process_MASL(self, Station):
        MASL = float(Station.MASL.item().replace("~", "")) if not Station.MASL.size == 0 else None
        return MASL

    def _process_Station_Open_date(self, Station):
        Station_Open_date = pandas.to_datetime(Station.Open_Date.item())
        return Station_Open_date

    def _process_Rain_instrument(self, Station):
        Rain_instrument = True if "גשם" in Station.vars.item() else False
        return Rain_instrument

    def _process_Temperature_instrument(self, Station):
        Temperature_instrument = True if "טמפ'" in Station.vars.item() else False
        return Temperature_instrument

    def _process_Wind_instrument(self, Station):
        Wind_instrument = True if "רוח" in Station.vars.item() else False
        return Wind_instrument

    def _process_Humidity_instrument(self, Station):
        Humidity_instrument = True if "לחות" in Station.vars.item() else False
        return Humidity_instrument

    def _process_Pressure_instrument(self, Station):
        Pressure_instrument = True if "לחץ" in Station.vars.item() else False
        return Pressure_instrument

    def _process_Radiation_instrument(self, Station):
        Radiation_instrument = True if "קרינה" in Station.vars.item() else False
        return Radiation_instrument

    def _process_Screen_Model(self,Station):
        Screen_Model=Station.Screen_Model.item()
        return Screen_Model

    def _process_InstLoc_AnemometeLoc(self,Station):
        InstLoc_AnemometeLoc=Station.Instruments_loc_and_Anemometer_loc.item()
        return InstLoc_AnemometeLoc

    def _process_Anemometer_h(self,Station):
        Anemometer_h=Station.Anemometer_height_m.item()
        return Anemometer_h

    def _process_comments(self,Station):
        comments=Station.comments.item()
        return comments

    def parse(self, path, station_column, time_column, metadatafile=None, **metadata):
        if os.path.isfile(path):
            df = pandas.read_json(path)
        else:
            all_files = glob.glob(os.path.join(path, "*.json"))

            L = []

            for filename in all_files:
                df = pandas.read_json(filename)
                L.append(df)

            df = pandas.concat(L, axis=0, ignore_index=True)

        df[time_column] = pandas.to_datetime(df[time_column])
        df = df.set_index(time_column)

        stations = [x for x in df[station_column].unique() if x not in self._removelist]

        metadata_dict = dict()
        for station in stations:
            filtered_stnname = "".join(filter(lambda x: not x.isdigit(), station)).strip()
            metadata_dict.setdefault(station, self._createMD(metadatafile, filtered_stnname, **metadata))

        # Remove problematic stations
        df = df.query('%s in @stations' % station_column)

        loaded_dask = dd.from_pandas(df, npartitions=1)
        return loaded_dask, metadata_dict

    def _createMD(self, metadatafile, StationName, **metadata):
        columns_dict = dict(BP='Barometric pressure[hPa]',
                            DiffR='Scattered radiation[W/m^2]',
                            Grad='Global radiation[W/m^2]',
                            NIP='Direct radiation[W/m^2]',
                            RH='Relative Humidity[%]',
                            Rain='Accumulated rain[mm/10minutes]',
                            STDwd='Wind direction standard deviation[degrees]',
                            TD='Average temperature in 10 minutes[c]',
                            TDmax='Maximum temperature in 10 minutes[c]',
                            TDmin='Minimum temperature in 10 minutes[c]',
                            TG='Average near ground temperature in 10 minutes[c]',
                            Time="End time of maximum 10 minutes wind running average[hhmm], see 'Ws10mm'",
                            WD='Wind direction[degrees]',
                            WDmax='Wind direction of maximal gust[degrees]',
                            WS='Wind speed[m/s]',
                            WS1mm='Maximum 1 minute average Wind speed[m/s]',
                            WSmax='Maximal gust speed[m/s]',
                            Ws10mm="Maximum 10 minutes wind running average[m/s], see 'Time''"
                            )

        columns_dict['StationName'] = StationName

        vals = dict()

        if metadatafile is not None:

            F = ['HebName', 'ITM_E', 'ITM_N', 'LAT_deg', 'LON_deg', 'MASL', 'Station_Open_date', 'Rain_instrument',
                 'Temperature_instrument', 'Wind_instrument', 'Humidity_instrument', 'Pressure_instrument',
                 'Radiation_instrument', 'Screen_Model', 'InstLoc_AnemometeLoc', 'Anemometer_h', 'comments'
                 ]

            MD = pandas.read_csv(metadatafile, delimiter="\t", names=["Serial_Num", "ENVISTA_ID", "Stn_name_Heb",
                                                                      "Stn_name_Eng", "ITM_E", "ITM_N", "Lon_deg",
                                                                      "Lat_deg", "MASL", "Open_Date", "vars",
                                                                      "Screen_Model", "Instruments_loc_and_Anemometer_loc",
                                                                      "Anemometer_height_m", "comments"
                                                                      ]
                                 )

            station = MD.query("Stn_name_Eng==@StationName")

            for x in F:
                updator = getattr(self, "_process_%s" % x)
                vals[x] = updator(station)

        vals.update(columns_dict)
        vals.update(metadata)
        return vals


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

    def parse(self, path, fromTime=None, **metadata):
        if os.path.isfile(path):
            df = self.getPandasFromFile(path, fromTime=fromTime)
        else:
            df = self.getPandasFromDir(path, fromTime=fromTime)

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

    def getPandasFromFile(self, path, fromTime):
        cbi = CampbellBinaryInterface(file=path)
        ts, cols, data = self.getData(path, fromTime=fromTime)
        dfList = []
        for i, key in enumerate(data.keys()):
            columns = cols[i]
            tmp_df = pandas.DataFrame(data[key], index=ts, columns=columns)
            tmp_df['height'] = int(key)
            tmp_df['station'] = cbi.headers[0].split(',')[1]
            tmp_df['instrument'] = cbi.headers[0].split(',')[-1]
            dfList.append(tmp_df)
        return pandas.concat(dfList, sort=True)

    def getPandasFromDir(self, path, fromTime):
        dfList = []
        for file in glob.glob(os.path.join(path, '*.dat')):
            tmp_df = self.getPandasFromFile(file, fromTime=fromTime)
            dfList.append(tmp_df)
        return pandas.concat(dfList, sort=True)

    def getData(self, file, fromTime):
        cbi = CampbellBinaryInterface(file=file)
        retVal = {}

        for i in range(len(cbi.columnsNames)):
            if len(cbi.columnsNames) == 1:
                retVal[10] = []
            else:
                retVal[6 + 5 * i] = []

        if type(fromTime)==str:
            fromTime = pandas.Timestamp(fromTime)

        recordIndex = 1 if fromTime is None else cbi.getRecordIndexByTime(fromTime)

        ts = []

        while recordIndex+1 <= cbi.recordsNum:
            time, line = cbi.getRecordByIndex(recordIndex)
            ts.append(time)

            for i, key in enumerate(retVal):
                retVal[key].append(line[cbi.columnsIndexes[i][0]: cbi.columnsIndexes[i][1]])
            recordIndex += 1

        return ts, cbi.columnsNames, retVal


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

    def _getTimeByIndex(self, i):
        time, _ = self.getRecordByIndex(i)
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

    def getRecordIndexByTime(self, time):
        upperIndex = self.recordsNum-1
        lowerIndex = 0
        recordTime = self._getTimeByIndex((lowerIndex+upperIndex)//2)
        tmpUpperIndex = -1
        tmpLowerIndex = -1

        while recordTime != time and (tmpLowerIndex != lowerIndex or tmpUpperIndex != upperIndex):
            tmpUpperIndex = upperIndex
            tmpLowerIndex = lowerIndex
            if time > recordTime:
                lowerIndex = (lowerIndex+upperIndex)//2
            else:
                upperIndex = (lowerIndex+upperIndex)//2
            recordTime = self._getTimeByIndex((lowerIndex+upperIndex)//2)

        return (lowerIndex+upperIndex)//2 if recordTime==time else None


class Parser_Radiosonde(object):

    def __init__(self):
        pass

    def parse(self, file):
        pass


class Parser_TOA5(object):

    def __init__(self):
        pass

    def parse(self, file):
        pass
