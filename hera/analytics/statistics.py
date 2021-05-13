import matplotlib.pyplot as plt
import numpy as np



def calcDist2d(x, y, data=None, bins=20, normalization="max_normalized"):
    """
    Calculates the distribution of two dimensional data.

    Parameters
    ----------
    data : Pandas
            the data.

    x: str or numpy array or series

    y: str or numpy array or series

    bins: number of bins for hist2d.

    normalization:  the normalization method: density or y_normalized or max_normalized

        max_normalized - normalize the data by the maximal value of the histogram to make 1 the maximum value.
        y_normalized   - normalize the data by group of x values to make the data proportional to the rest of the group values.
        density        - normalize the data by the dXdY of the data. assume the data is equidistant.


    Returns
    -------
    tuple with 3 values:
    (x_mid,y_mid,M.T)

     x_mid: The bin center along the x axis.
     y_mid: The bin center along the y axis.
     M.T:   Transpose of the 2D histogram.
    """

    tmpfig = plt.figure()


    if data is not None:

        xdata = data[x];
        ydata = data[y]

    else:
        xdata = x;
        ydata = y;

    M, x_vals, y_vals, _ = plt.hist2d(xdata, ydata, bins=bins)
    plt.close(tmpfig)

    if normalization == "density":
        square_area = (x_vals[1] - x_vals[0]) * (y_vals[1] - y_vals[0])
        M = M / square_area

    elif normalization == "y_normalized":
        M = M / M.sum(axis=1)[:, None]
        M[np.isnan(M)] = 0

    elif normalization == "max_normalized":
        M = M / M.max()

    else:
        raise ValueError("The normaliztion must be either density,y_normalized or max_normalized")

    x_mid = (x_vals[1:] + x_vals[:-1]) / 2
    y_mid = (y_vals[1:] + y_vals[:-1]) / 2

    return x_mid, y_mid, M.T