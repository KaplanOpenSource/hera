import os
import xarray
import numpy
import os
from hera.utils.unitHandler import ureg, unumToPint

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

    def getDosage(self, Q=1 * ureg.kg, time_units=ureg.min, q_units=ureg.mg):
        """
        Calculates the dosage from the LSM simulation.
        The default are pint units, but we still support unum inputs.

        The outputs (in the attrs) are always pint.

        Parameters
        ----------
        Q : pint|unum.units
            Default value is 1*kg

        time_units: pint|unum.units
            Default value is min

        q_units: pint|unum.units
            Default value is mg

        Returns
        -------
        self._finalxarray: xarray
            The calculated dosage in 'Dosage' key
        """
        Q = unumToPint(Q)
        time_units = unumToPint(time_units)
        q_units = unumToPint(q_units)

        final_xarray = self._finalxarray.copy()
        if type(final_xarray.datetime.diff('datetime')[0].values.item())==float:
            dt_minutes = final_xarray.datetime.diff('datetime')[0].values.item()*ureg.sec #temporary solution!!!!!
        else:
            dt_minutes = (final_xarray.datetime.diff('datetime')[0].values / numpy.timedelta64(1, 'm')) * ureg.min
        final_xarray.attrs['dt'] = dt_minutes.to(time_units)
        final_xarray.attrs['Q']  = Q.to(q_units)
        final_xarray.attrs['C']  = q_units/ ureg.m ** 3

        Qfactor = (Q.to(q_units) * ureg.min / ureg.m ** 3).m_as(q_units * time_units / ureg.m ** 3)

        final_xarray['Dosage']   = Qfactor*final_xarray['Dosage']

        return final_xarray

    def getConcentration(self, Q=1*ureg.kg, time_units=ureg.min, q_units=ureg.mg):
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
        time_units=unumToPint(time_units)
        Q = unumToPint(Q)
        q_units=unumToPint(q_units)
        finalxarray = self.getDosage(Q=Q, time_units=time_units, q_units=q_units)

        dDosage = finalxarray['Dosage'].diff('datetime').to_dataset().rename({'Dosage': 'dDosage'})
        dDosage['C'] = dDosage['dDosage'] / finalxarray.attrs['dt'].m_as(time_units)
        dDosage.attrs = finalxarray.attrs

        return dDosage

    def getConcentrationAtPoint(self, x, y, datetime, Q=1*ureg.kg, time_units=ureg.min, q_units=ureg.mg):
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
