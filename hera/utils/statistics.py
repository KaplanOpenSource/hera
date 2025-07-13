import matplotlib.pyplot as plt
import numpy as np

def calcDist2d(x, y, data=None, bins=20, normalization="max_normalized", x_range=None, y_range=None):
    """
    Calculates the distribution of two dimensional data with optional normalization and axis range.

    Parameters
    ----------
    data : DataFrame or None
        Optional DataFrame from which to extract x and y columns.

    x : str or array-like
        x-values or column name in data.

    y : str or array-like
        y-values or column name in data.

    bins : int
        Number of bins along each axis.

    normalization : str
        "density" - normalize by bin area;
        "y_normalized" - normalize each row to sum to 1;
        "max_normalized" - normalize so the max is 1.

    x_range : tuple (optional)
        Lower and upper limits for x-axis bins (e.g. (0, 24)).

    y_range : tuple (optional)
        Lower and upper limits for y-axis bins (e.g. (0, 100)).

    Returns
    -------
    x_mid : np.array
        Midpoints of x bins

    y_mid : np.array
        Midpoints of y bins

    M.T : np.array
        Transposed 2D histogram matrix after normalization
    """

    tmpfig = plt.figure()

    # Extract data
    if data is not None:
        xdata = data[x]
        ydata = data[y]
    else:
        xdata = x
        ydata = y

    # Set bin ranges if specified
    hist_range = None
    if x_range is not None and y_range is not None:
        hist_range = [x_range, y_range]

    # Compute the histogram
    M, x_vals, y_vals, _ = plt.hist2d(xdata, ydata, bins=bins, range=hist_range)
    plt.close(tmpfig)

    # Normalization
    if normalization == "density":
        square_area = (x_vals[1] - x_vals[0]) * (y_vals[1] - y_vals[0])
        M = M / square_area



    elif normalization == "y_normalized":
        M = M.astype(float)
        col_sums = M.sum(axis=0, keepdims=True)
        M = np.divide(M, col_sums, out=np.zeros_like(M), where=(col_sums != 0))



    elif normalization == "max_normalized":
        max_val = M.max()
        if max_val > 0:
            M = M / max_val

    else:
        raise ValueError("The normalization must be one of: density, y_normalized, max_normalized")

    # Midpoints for each bin
    x_mid = (x_vals[1:] + x_vals[:-1]) / 2
    y_mid = (y_vals[1:] + y_vals[:-1]) / 2

    return x_mid, y_mid, M.T
