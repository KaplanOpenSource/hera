import os
import xarray
import numpy
from unum.units import *

from ..utils import toUnum,toNumber

# def toVTK(self, data, outputdir, name, fields):
#     from pyevtk.hl import gridToVTK
#
#     for ii, tt in enumerate(data.datetime):
#         curdata = data.sel(datetime=tt)
#         outputpath = os.path.join(outputdir, "%s_%s" % (name, ii))
#         fieldsmap = dict([(key, curdata[key].values) for key in fields])
#
#         gridToVTK(outputpath, \
#                   curdata.x.values, \
#                   curdata.y.values, \
#                   curdata.z.values, \
#                   pointData=fieldsmap)

class SingleSimulation(object):
    _finalxarray = None
    _document = None

    @property
    def params(self):
        return self._document['desc']['params']

    @property
    def version(self):
        return self._document['desc']['version']

    def __init__(self, resource):
        if type(resource) is str:
            # pattern = os.path.join(resource, '*.nc')
            #
            # filenameList = []
            # times = []
            # for infilename in glob.glob(pattern):
            #     filenameList.append(infilename)
            #     timepart = float(".".join(infilename.split(".")[0].split("_")[-2:]))
            #     times.append(timepart)
            #
            # # Sort according to time.
            # combined = sorted([x for x in zip(filenameList, times)], key=lambda x: x[1])
            #
            # self._xray = xarray.open_mfdataset([x[0] for x in combined])
            self._finalxarray = xarray.open_mfdataset(os.path.join(resource, '*.nc'), combine='by_coords')
        else:
            self._document = resource
            self._finalxarray = resource.getData()
            if type(self._finalxarray) is str:
                self._finalxarray = xarray.open_mfdataset(self._finalxarray, combine='by_coords')

    def getDosage(self, Q=1 * kg, time_units=min, q_units=mg):
        """
        Calculates the dosage

        Parameters
        ----------
        Q : unum.units
            Default value is 1*kg

        time_units: unum.units
            Default value is min

        q_units: unum.units
            Default value is mg

        Returns
        -------
        self._finalxarray: xarray
            The calculated dosage in 'Dosage' key
        """
        if 'dt' not in self._finalxarray.attrs.keys():
            dt_minutes = (self._finalxarray.datetime.diff('datetime')[0].values / numpy.timedelta64(1, 'm')) * min
            self._finalxarray.attrs['dt'] = toUnum(dt_minutes, time_units)
            self._finalxarray.attrs['Q']  = toUnum(Q, q_units)
            self._finalxarray.attrs["C"] = toUnum(1, q_units/(m**3))

            Qfactor = toNumber(self._finalxarray.attrs['Q'] * min / m ** 3,
                               q_units * time_units / m ** 3)

            self._finalxarray['Dosage']   = Qfactor*self._finalxarray['Dosage']

        return self._finalxarray.copy()

    def getConcentration(self, Q=1*kg, time_units=min, q_units=mg):
        """
        Calculates the concentration

        Parameters
        ----------
        Q : unum.units
            Default value is 1*kg

        time_units: unum.units
            Default value is min

        q_units: unum.units
            Default value is mg

        Returns
        -------
        dDosage: xarray
            The calculated concentration in 'C' key
        """
        if 'dt' not in self._finalxarray.attrs.keys():
            self.getDosage(Q=Q, time_units=time_units, q_units=q_units)

        dDosage = self._finalxarray['Dosage'].diff('datetime').to_dataset().rename({'Dosage': 'dDosage'})
        dDosage['C'] = dDosage['dDosage'] / self._finalxarray.attrs['dt'].asNumber()
        dDosage.attrs = self._finalxarray.attrs

        return dDosage.copy()

    def getConcentrationAtPoint(self, x, y, datetime, Q=1*kg, time_units=min, q_units=mg):
        """
        Calculates the concentration at requested point and time

        Parameters
        ----------
        x: int
            latitude of the point

        y: int
            longitude of the point

        datetime: datetime
            time of the calculation

        Returns
        -------
        con: float
            The concentration at the requested point and time
        """
        return self.getConcentration(Q=Q, time_units=time_units, q_units=q_units)['C'].interp(x=x, y=y, datetime=datetime).values[0]