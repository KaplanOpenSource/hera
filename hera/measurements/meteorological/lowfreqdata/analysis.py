from hera.analytics.statistics import calcDist2d
import numpy
import pandas
import dask

WINTER = 'Winter'
SPRING = 'Spring'
SUMMER = 'Summer'
AUTUMN = 'Autumn'

seasonsdict = dict(WINTER=dict(monthes=[12,1,2],strmonthes='[DJF]'),
                   SPRING=dict(monthes=[3,4,5],strmonthes='[MAM]'),
                   SUMMER=dict(monthes=[6,7,8],strmonthes='[JJA]'),
                   AUTUMN=dict(monthes=[9,10,11],strmonthes='[SOM]')
                   )


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

    if isinstance(data,dask.dataframe):
        curdata = curdata.map_partitions(lambda df: df.assign(season=tm(df,monthcolumn)))
    else:
        curdata=curdata.assign(season=tm(curdata, monthcolumn))

    return curdata

def calcHourlyDist(data, Field,bins=30, normalization='density'):
    """

    :param data:
    :param Field:
    :param bins:
    :param normalization:
    :return:
    """

    curdata = data.dropna(subset=[Field])

    curdata = curdata.query("%s > -9990" % Field)
    curdata=curdata.assign(curdate=curdata.index)
    curdata=curdata.assign(houronly=curdata.curdate.dt.hour + curdata.curdate.dt.minute / 60.)

    y = curdata[Field]
    x = curdata['houronly']

    return calcDist2d(x=x,y=y,bins=bins,normalization=normalization)

def percentiledask(data, field, lower=None,upper=None):
    """

    :param data:
    :param method:
    :param quantilefield:
    :param quantiles:
    :return:
    """

    quantiles=[]

    lowerExists = False
    upperExists = False

    if lower is not None:
        quantiles.append(lower)
        lowerExists = True

    if upper is not None:
        quantiles.append(upper)
        upperExists = True

    if (not lowerExists and not upperExists):
        raise ValueError("Must supply a lower or an upper limit.")

    vals=numpy.atleast1d(quantiles).compute()


    if (lowerExists and upperExists):
        qry = "%s<%s<%s" % (vals[lower],field,vals[upper])
    elif(lowerExists and not upperExists):
        qry = "%s<%s" % (vals[lower], field)
    else:
        qry = "%s<%s" % (field,vals[upper])

    return data.query(qry)






