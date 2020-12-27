from hera.analytics.statistics import calcDist2d
import numpy
import pandas
import dask

WINTER = 'Winter'
SPRING = 'Spring'
SUMMER = 'Summer'
AUTUMN = 'Autumn'

seasonsdict = {    WINTER : dict(months=[12,1,2],strmonths='[DJF]'),
                   SPRING : dict(months=[3,4,5],strmonths='[MAM]'),
                   SUMMER : dict(months=[6,7,8],strmonths='[JJA]'),
                   AUTUMN : dict(months=[9,10,11],strmonths='[SOM]')}


def addDatesColumns(data, datecolumn=None, monthcolumn=None):
    """
        This class adds the year,month, date, time, and the season to the dataframe.

    parameters
    -----------
    data: dask.dataframe, pandas.Datafram
        The data to add to.

    datecolumn: str
        The name of the column with the date.
        If None, use the index.
    monthcolumn:

    returns
    --------
        dask.dataframe, pandas.Datafram

    """

    curdata = data

    if datecolumn is None:
        curdata = curdata.assign(curdate=curdata.index)
        datecolumn='curdate'

    curdata = curdata.assign(yearonly=curdata[datecolumn].dt.year)

    if monthcolumn==None:
        curdata = curdata.assign(monthonly=curdata[datecolumn].dt.month)
        monthcolumn = 'monthonly'

    curdata = curdata.assign(dayonly=curdata[datecolumn].dt.day)\
                     .assign(timeonly=curdata[datecolumn].dt.time)

    tm = lambda x, field: pandas.cut(x[field], [0, 2, 5, 8, 11, 12], labels=['Winter', 'Spring', 'Summer', 'Autumn', 'Winter1']).replace('Winter1', 'Winter')

    if isinstance(data,dask.dataframe.DataFrame):
        curdata = curdata.map_partitions(lambda df: df.assign(season=tm(df,monthcolumn)))
    else:
        curdata=curdata.assign(season=tm(curdata, monthcolumn))

    return curdata

def calcHourlyDist(data, Field,bins=30, normalization='density'):
    """
            Calculates hours distribution of the field.

    Parameters
    ----------

    data: pandas.DataFrame or dask.DataFrame
            The data to calculate.
            We assume that the index is a datetime object.

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

    """

    curdata = data.dropna(subset=[Field])

    curdata = curdata.query("%s > -9990" % Field)
    curdata=curdata.assign(curdate=curdata.index)
    curdata=curdata.assign(houronly=curdata.curdate.dt.hour + curdata.curdate.dt.minute / 60.)

    y = curdata[Field]
    x = curdata['houronly']

    return calcDist2d(x=x,y=y,bins=bins,normalization=normalization)


# def percentiledask(data, field, lower=None,upper=None):
#     """
#
#     Parameters
#     ----------
#
#     data:
#
#     method:
#
#     quantilefield:
#
#     quantiles:
#
#     Returns
#     -------
#         Return
#
#     """
#
#     quantiles=[]
#
#     lowerExists = False
#     upperExists = False
#
#     if lower is not None:
#         quantiles.append(lower)
#         lowerExists = True
#
#     if upper is not None:
#         quantiles.append(upper)
#         upperExists = True
#
#     if (not lowerExists and not upperExists):
#         raise ValueError("Must supply a lower or an upper limit.")
#
#     vals=numpy.atleast1d(quantiles).compute()
#
#
#     if (lowerExists and upperExists):
#         qry = "%s<%s<%s" % (vals[lower],field,vals[upper])
#     elif(lowerExists and not upperExists):
#         qry = "%s<%s" % (vals[lower], field)
#     else:
#         qry = "%s<%s" % (field,vals[upper])
#
#     return data.query(qry)
#





