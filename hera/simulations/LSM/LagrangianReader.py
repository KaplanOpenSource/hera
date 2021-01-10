import pandas
import numpy 
import xarray
import glob 
import os
from unum.units import *


def toNetcdf(basefiles, addzero=True, datetimeFormat="timestamp"):
    """
        Converts the data to netcdf.
        The dosage are converted to s/m**3 instead of min/m**3.

        Parameters
        ----------
        basefiles: str
            Path to the directory with the netcdf files

        addZero: bool
            if true, adds a 0 file at the begining of the files (with time shift 0)

    """

    #outfilename = name
    filenameList = []
    times = []
    for infilename in glob.glob(os.path.join("%s*" % basefiles)):
        filenameList.append(infilename)
        times.append(float(infilename.split("_")[-1]))

    print("Processing the files")
    print(basefiles)
    print([x for x in glob.glob(os.path.join("%s*" % basefiles))])
    # Sort according to time.
    combined = sorted([x for x in zip(filenameList, times)], key=lambda x: x[1])
    dt = None

    for (i,curData) in enumerate(combined):
        print("\t... reading %s" % curData[0])
        cur = pandas.read_csv(curData[0], delim_whitespace=True, names=["y", "x", "z", "Dosage"]) #,dtype={'x':int,'y':int,'z':int,'Dosage':float})
        cur['time'] = curData[1]

        if dt is None:
            dt = cur.iloc[0]['time'] * s

        cur['Dosage'] *= (s / m ** 3).asNumber(min / m ** 3)
        xray = cur.sort_values(['time', 'x', 'y', 'z']).set_index(['time', 'x', 'y', 'z']).to_xarray()
        if datetimeFormat== "timestamp":
            datetime = pandas.to_datetime("1-1-2016 12:00") + pandas.to_timedelta(xray.time.values, 's')
        elif datetimeFormat== "seconds":
            datetime = xray.time.values

        #finalxarray.to_netcdf(os.path.join(topath,name,"%s_%s.nc" % (outfilename, str(cur['time'].iloc[0]).replace(".", "_"))) )

        if (i==0) and addzero:
            if datetimeFormat== "timestamp":
                zdatetime = [pandas.to_datetime("1-1-2016 12:00")]
            elif datetimeFormat == "seconds":
                zdatetime = [0]
      
            finalxarray = xarray.DataArray(numpy.zeros(xray['Dosage'].values.shape), \
                                       coords={'x': xray.x, 'y': xray.y, 'z': xray.z, 'datetime': zdatetime},
                                       dims=['datetime', 'y', 'x', 'z']).to_dataset(name='Dosage')

            yield finalxarray
           #finalxarray.to_netcdf(os.path.join(topath,name,"%s_0_0.nc" % outfilename) )

        finalxarray = xarray.DataArray(xray['Dosage'].values, \
                                       coords={'x': xray.x, 'y': xray.y, 'z': xray.z, 'datetime': datetime},
                                       dims=['datetime', 'y', 'x', 'z']).to_dataset(name='Dosage')

        yield finalxarray