import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

class WindPlots(object):
    def plot_hist2d(self, data = None, y = None, x = "time_within_day", axis = None, ax_props = None, bins = 30, levels = 20, max_normalized = True, **kwargs):
        """

            plot_hist2d_windVec2(data=..,y="WD",month=3,hour=1)
                        query field is month,hour

                        querystr = ["month==3","hour==1"] ==> "month==3 and hour==1"

            plot_hist2d_windVec2(data=..,y="WD",month=3)
                        query field is month
                        querystr = ["month==3"]


            :param data:
            :param y:
            :param x:
            If data is None (default) then y,x are interpreted as vector data; else interpreted as names of columns in data.

            :param axis:
            :param ax_props : dict with axis functions and parameters. { funcname : {parameters}}
            :param bins: the number of bins in the histogram, default is 30
            :param levels: the number of height levels in the contour plot, default is 20
            :param max_normalized: whether or not the density will be normalized to make 1 the maximumum value

            :param kwargs:
                    all the rest are used as query (use the and predicate on all the fields).
                    month=1
                    day = 2


            :return:
            """

        # kwargs_used_fields = ["axis", "fig", "bins", "levels"]
        # queryfield = [x for x in kwargs.key() if x not in kwargs_used_fields]
        # querystr = " and ".join(["{fieldname}=={fieldvalue}".format(fieldname=fieldname, fieldvalue=kwargs[fieldname]) for fieldname in queryfield])
        #
        # month_data = data.assign(year=data.index.year) \
        #     .assign(month=data.index.month) \
        #     .assign(day=data.index.day) \
        #     .assign(hour=data.index.hour) \
        #     .query(querystr)

        tmpfig = plt.figure()
        newax = plt.subplot(111)
        if data is None:
            xdata = x; ydata = y
        else:
            xdata = data[x]; ydata = data[y]

        M, x_vals, y_vals, newax = plt.hist2d(xdata, ydata, bins=bins)
        plt.close(tmpfig)
        if axis is None:
            fig, axis = plt.subplots()
        else:
            plt.sca(axis)

        if max_normalized:
            M = M / M.max()
        x_mid = (x_vals[1:] + x_vals[:-1]) / 2
        y_mid = (y_vals[1:] + y_vals[:-1]) / 2
        plt.contourf(x_mid, y_mid, M.transpose(), levels, cmap=plt.cm.jet)

        if not ax_props is None:
            for func in ax_props:
                getattr(axis, func)(**ax_props[func])
