import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from spt3g import core, calibration
from spt3g.std_processing import obsid_to_g3time
import datetime as dt

# plot median calibrator S/N for a list of observations; this is a timeseries plot
def MedianCalSN(request):
    bands = [90, 150, 220]
    
    median_calSN = {90: [], 150: [], 220: []}
    
    for obsid in request['observation']:
        try:
            data = [fr for fr in core.G3File('{}/calibration/{}/{}.g3' \
                                             .format(request['caldatapath'],
                                                     request['source'],
                                                     obsid))]
            boloprops = [fr for fr in core.G3File('{}/downsampled/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          obsid))] \
                                                  [0]["NominalBolometerProperties"]
        except RuntimeError:
            return "Could not find data files."

        for band in bands:
            try:
                cal_data = np.array([data[0]['CalibratorResponseSN'][bolo] \
                                for bolo in data[0]['CalibratorResponseSN'].keys() \
                                if boloprops[bolo].band / core.G3Units.GHz == band])
            except KeyError:
                return "CalibratorResponseSN does not exist for this observation."
            median_calSN[band].append(np.median(cal_data[np.isfinite(cal_data)]))

    timestamps = np.sort([obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                        for obsid in request['observation']])
    dts = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
    datenums = mdates.date2num(dts)
    fig, ax = plt.subplots(1)
    for band in bands:
        ax.plot(datenums, median_calSN[band], 'o-', 
                label='{} GHz'.format(band))
    xfmt = mdates.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    plt.xticks(rotation=25)
    plt.legend()
    plt.xlabel('observation time')
    plt.ylabel('median calibration S/N')
    plt.title('Calibrator S/N')
    plt.tight_layout()

    return fig
