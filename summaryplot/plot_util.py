import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

colors = {90:'C0', 150:'C1', 220:'C2'}

def plot_timeseries(datetimes, data, ylims, band):
    '''
    Make a timeseries plot with some special features appropriate to the
    summary web page.

    Parameters
    ----------
    datetimes : numpy array
        Datetimes corresponding to each observation to plot
    data : numpy array
        Data for each observation to plot
    ylims : 2-tuple
        [ymin, ymax] containing the min and max y-values to include in the plot
    band : int
        Observing band (90, 150 or 220)

    Returns
    -------
    None
    '''
    # plot data
    datenums = mdates.date2num(datetimes)
    inrange = (data[np.isfinite(data)] >= ylims[0]) & \
              (data[np.isfinite(data)] <= ylims[1])
    plt.plot(datenums[np.isfinite(data)][inrange],
             data[np.isfinite(data)][inrange],
             'o', label='{} GHz'.format(band), color=colors[band])
    
    # plot out-of-range data                                                                                           
    above_range = (data[np.isfinite(data)] > ylims[1])
    below_range = (data[np.isfinite(data)] < ylims[0])
    plt.plot(datenums[np.isfinite(data)][above_range],
             ylims[1]*np.ones(len(datenums[np.isfinite(data)][above_range])),
             '^', color='C3')
    plt.plot(datenums[np.isfinite(data)][below_range],
             ylims[0]*np.ones(len(datenums[np.isfinite(data)][below_range])),
             'v', color='C3')

    # plot light dashed lines for NaNs                                             
    nan_dates = datenums[~np.isfinite(data)]
    for date in nan_dates:
        plt.plot([date, date], ylims, 'r--', linewidth=0.5)

    # set limits
    plot_ymin = ylims[0] - 0.01*(ylims[1]-ylims[0])
    plot_ymax = ylims[1] + 0.01*(ylims[1]-ylims[0])
    plt.ylim([plot_ymin, plot_ymax])
