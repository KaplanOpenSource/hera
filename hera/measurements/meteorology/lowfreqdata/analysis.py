from hera.utils.statistics import calcDist2d
import numpy
import pandas
import dask
import pandas as pd


WINTER = 'Winter'
SPRING = 'Spring'
SUMMER = 'Summer'
AUTUMN = 'Autumn'

seasonsdict = {    WINTER : dict(months=[12,1,2],strmonths='[DJF]'),
                   SPRING : dict(months=[3,4,5],strmonths='[MAM]'),
                   SUMMER : dict(months=[6,7,8],strmonths='[JJA]'),
                   AUTUMN : dict(months=[9,10,11],strmonths='[SON]')}



class analysis:


    _datalayer = None

    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self, dataLayer):
        self._datalayer = dataLayer

    def addDatesColumns(self, data, datecolumn=None, monthcolumn=None):
        """
            This class adds the year, month, date, time, and the season to the dataframe.

        parameters
        -----------
        data: dask.dataframe, pandas.DataFrame or str
            The data to add to. If str, it is assumed to be a path to a Parquet file.

        datecolumn: str
            The name of the column with the date.
            If None, use the index.

        monthcolumn:

        returns
        --------
            dask.dataframe, pandas.DataFrame

        Notes
        ------
        ✨ Added support: if `data` is a string path, we automatically read it as a Parquet file using pandas.
        """

        # NEW: Support for string input (e.g., path to parquet file)
        if isinstance(data, str):
            import pandas as pd
            data = pd.read_parquet(data)

        curdata = data

        if datecolumn is None:
            curdata = curdata.assign(curdate=curdata.index)
            datecolumn = 'curdate'

        curdata = curdata.assign(yearonly=curdata[datecolumn].dt.year)

        if monthcolumn is None:
            curdata = curdata.assign(monthonly=curdata[datecolumn].dt.month)
            monthcolumn = 'monthonly'

        curdata = curdata.assign(dayonly=curdata[datecolumn].dt.day) \
            .assign(timeonly=curdata[datecolumn].dt.time)

        tm = lambda x, field: pandas.cut(x[field], [0, 2, 5, 8, 11, 12],
                                         labels=['Winter', 'Spring', 'Summer', 'Autumn', 'Winter1']).replace('Winter1',
                                                                                                             'Winter')

        if isinstance(data, dask.dataframe.DataFrame):
            curdata = curdata.map_partitions(lambda df: df.assign(season=tm(df, monthcolumn)))
        else:
            curdata = curdata.assign(season=tm(curdata, monthcolumn))

        return curdata

    def calcHourlyDist(self, data, Field, bins=30, normalization='density'):
        """
            Calculates hours distribution of the field.

        Parameters
        ----------

        data: pandas.DataFrame or dask.DataFrame or str
                The data to calculate.
                We assume that the index is a datetime object.
                If str, it is assumed to be a path to a Parquet file.

        Field: str
                The name of the column to calculate the statistics on.

        bins: int
                The number of bins to calculate

        normalization: str
            max_normalized - normalize the data by the maximal value of the histogram to make 1 the maximum value.
            y_normalized   - normalize the data by group of x values to make the data proportional to the rest of the group values.
            density        - normalize the data by the dXdY of the data. assume the data is equidistant.

        Returns
        --------

        tuple with 3 values:
        (x_mid,y_mid,M.T)

         x_mid: The bin center along the x axis.
         y_mid: The bin center along the y axis.
         M.T:   Transpose of the 2D histogram.

        Notes
        ------
        ✅ Added support for string input: if `data` is a string (path), we automatically read it as a Parquet file.
        """

        # ✅ New: Support for string input (path to Parquet file)
        if isinstance(data, str):
            import pandas as pd
            data = pd.read_parquet(data)

        curdata = data.dropna(subset=[Field])

        curdata[Field] = curdata[Field].where(curdata[Field] > -5000)
        # curdata = curdata.query("%s > -9990" % Field)

        curdata = curdata.assign(curdate=curdata.index)
        curdata.curdate = pd.to_datetime(curdata.curdate, utc=True)
        curdata = curdata.assign(houronly=curdata.curdate.dt.hour + curdata.curdate.dt.minute / 60.)

        curdata = curdata.dropna()
        y = curdata[Field]
        x = curdata['houronly']

        # ✅ Added fixed range support to avoid histogram truncation
        # This ensures consistent x-axis (0–24 hours) and y-axis range that includes all values
        x_range = (0, 24)
        y_min = min(y.min(), 0)
        y_max = max(y.max(), 1)  # fallback to 1 if all values are near 0
        y_range = (y_min, y_max)

        return calcDist2d(x=x, y=y, bins=bins, normalization=normalization, x_range=x_range, y_range=y_range)

    def _calculateCov(self, data, data_resampled, x, y, SamplingWindow):

        data_resampled[x + y] = data[x + y].resample(SamplingWindow).mean() + \
                              (data[x + "_bar"] * data[y + "_bar"]).resample(SamplingWindow).mean() \
                              - data_resampled[x + "_bar"] * data_resampled[y + "_bar"]

        return self

    def resampleSecondMoments(self, data, SamplingWindow, fieldsFirstMoments, fieldsSecondMoments):
        """

        :param data:
        :param samplingWindow:
        :param fieldsFirstMoments:
        :param fieldsSecondMoments:
        :return:
        """

        data_resampled = data[fieldsFirstMoments].resample(SamplingWindow).mean()

        for i in range(len(fieldsSecondMoments)):
            for j in range(i, len(fieldsSecondMoments)):
                self._calculateCov(data, data_resampled, fieldsSecondMoments[i], fieldsSecondMoments[j], SamplingWindow)

        return data_resampled