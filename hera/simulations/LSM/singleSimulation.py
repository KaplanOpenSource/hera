import os
import xarray
import numpy
import os
from hera.utils.unitHandler import  *

from ...utils import tounit,tonumber


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

        if isinstance(resource,str):
            try:
                self._finalxarray = xarray.open_mfdataset(os.path.join(resource, '*.nc'), combine='by_coords')
            except OSError:
                self._finalxarray = xarray.open_mfdataset(os.path.join(resource,"netcdf", '*.nc'), combine='by_coords')
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
        finalxarray = self._finalxarray.copy()
        if type(finalxarray.datetime.diff('datetime')[0].values.item())==float:
            dt_minutes = finalxarray.datetime.diff('datetime')[0].values.item()*s #temporary solution!!!!!
        else:
            dt_minutes = (finalxarray.datetime.diff('datetime')[0].values / numpy.timedelta64(1, 'm')) * min
        finalxarray.attrs['dt'] = tounit(dt_minutes, time_units)
        finalxarray.attrs['Q']  = tounit(Q, q_units)
        finalxarray.attrs['C']  = tounit(1, q_units/ m ** 3)

        Qfactor = tonumber(finalxarray.attrs['Q'] * min / m ** 3,
                           q_units * time_units / m ** 3)

        finalxarray['Dosage']   = Qfactor*finalxarray['Dosage']

        return finalxarray

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

        finalxarray = self.getDosage(Q=Q, time_units=time_units, q_units=q_units)

        dDosage = finalxarray['Dosage'].diff('datetime').to_dataset().rename({'Dosage': 'dDosage'})
        dDosage['C'] = dDosage['dDosage'] / finalxarray.attrs['dt'].asNumber()
        dDosage.attrs = finalxarray.attrs

        return dDosage

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
