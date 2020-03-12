import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

colors = {90:'C0', 150:'C1', 220:'C2'}

def plot_timeseries(datetimes, data, band, xlims=None, ylims=None, alpha=0.35):
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
    plt.plot(datenums[np.isfinite(data)],
             data[np.isfinite(data)],
             'o', label='{} GHz'.format(band), color=colors[band],
             alpha=alpha)
    
    # plot out-of-range data and set limits
    if xlims != None:
        plot_xmin = xlims[0] - 0.01*(xlims[1]-xlims[0])
        plot_xmax = xlims[1] + 0.01*(xlims[1]-xlims[0])
        plt.xlim([plot_xmin, plot_xmax])
    if ylims != None:
        above_range = (data[np.isfinite(data)] > ylims[1])
        below_range = (data[np.isfinite(data)] < ylims[0])
        plt.plot(datenums[np.isfinite(data)][above_range],
                 ylims[1]*np.ones(len(datenums[np.isfinite(data)][above_range])),
                 '^', color='C3')
        plt.plot(datenums[np.isfinite(data)][below_range],
                 ylims[0]*np.ones(len(datenums[np.isfinite(data)][below_range])),
                 'v', color='C3')
        plot_ymin = ylims[0] - 0.01*(ylims[1]-ylims[0])
        plot_ymax = ylims[1] + 0.01*(ylims[1]-ylims[0])
        plt.ylim([plot_ymin, plot_ymax])

    # plot light dashed lines for NaNs                                             
    nan_dates = datenums[~np.isfinite(data)]
    ylims_current = plt.gca().get_ylim()
    for date in nan_dates:
        plt.plot([date, date], ylims_current, 'r--', linewidth=0.5)
    plt.ylim(ylims_current)
