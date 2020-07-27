import pandas


WINTER = 'Winter'
SPRING = 'Spring'
SUMMER = 'Summer'
AUTUMN = 'Autumn'

seasonsdict = dict(Winter=dict(monthes=[12,1,2]),
                   Spring=dict(monthes=[3,4,5]),
                   Summer=dict(monthes=[6,7,8]),
                   Autumn=dict(monthes=[9,10,11])
                   )


def addDatesColumns(data, dataType='dask', datecolumn=None, monthcolumn=None):
    """
    parameters
    -----------
    data:
    dataType:
    datecolumn:
    monthcolumn:

    returns
    --------
    curdata:
    """

    curdata = data

    if datecolumn is None:
        curdata = curdata.assign(curdate=curdata.index)
        datecolumn='curdate'

    curdata = curdata.assign(yearonly=curdata[datecolumn].dt.year)

    if monthcolumn==None:
        curdata = curdata.assign(monthonly=curdata[datecolumn].dt.month)
        monthcolumn = 'monthonly'

    curdata = curdata.assign(dayonly=curdata[datecolumn].dt.day)

    curdata = curdata.assign(timeonly=curdata[datecolumn].dt.time)

    tm = lambda x, field: pandas.cut(x[field], [0, 2, 5, 8, 11, 12], labels=['Winter', 'Spring', 'Summer', 'Autumn', 'Winter1']).replace('Winter1', 'Winter')

    if dataType=='dask':
        curdata = curdata.map_partitions(lambda df: df.assign(season=tm(df,monthcolumn)))

    else:
        curdata=curdata.assign(season=tm(curdata, monthcolumn))

    return curdata