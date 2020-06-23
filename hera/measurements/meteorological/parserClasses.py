import os
import glob
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

    def __init__(self):
        pass

    def parse(self, file):
        pass


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
