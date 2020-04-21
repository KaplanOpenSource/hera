from ....analytics.statistics import calcDist2d
import numpy
from datetime import datetime

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










