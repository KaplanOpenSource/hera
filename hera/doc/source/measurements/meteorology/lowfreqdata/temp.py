stationName = 'BET DAGAN'
from hera.measurements.meteorology import lowfreqdata 
datadb = lowfreqdata.lowfreqDataLayer.getStationDataFromDB(StationName=stationName)
lowfreqdata.analysis.addDatesColumns(datadb.getData())
