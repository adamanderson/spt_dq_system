import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from spt3g import core, calibration
from spt3g.std_processing import obsid_to_g3time

# plot median calibrator S/N for a list of observations; this is a timeseries plot
def MedianCalSN(request):
    bands = [90, 150, 220]
    
    median_calSN = {90: [], 150: [], 220: []}
    
    for obsid in request['observation']:
        try:
            data = [fr for fr in core.G3File('/spt/user/production/calibration/{}/{}.g3' \
                                                 .format(request['source'], obsid))]
            boloprops = [fr for fr in core.G3File('/spt/data/bolodata/downsampled/{}/{}/nominal_online_cal.g3' \
                                                      .format(request['source'], obsid))] \
                                                      [0]["NominalBolometerProperties"]
        except RuntimeError:
            return "Could not find data files."

        for band in bands:
            try:
                cal_data = [data[0]['CalibratorResponseSN'][bolo] \
                                for bolo in data[0]['CalibratorResponseSN'].keys() \
                                if boloprops[bolo].band / core.G3Units.GHz == band]
            except KeyError:
                return "CalibratorResponseSN does not exist for this observation."
            median_calSN[band].append(np.median(cal_data))

    obstimes = [obsid_to_g3time(int(obsid)) for obsid in request['observation']]
    fig, ax = plt.subplots(1)
    for band in bands:
        ax.plot(obstimes, median_calSN[band], label='{} GHz'.format(band))
    ax.format_xdata = mdates.DateFormatter('%m-%d %h:%m')
    plt.legend()
    plt.tight_layout()

    plt.xlabel('observation time')
    plt.ylabel('median calibration S/N')
    plt.title('Calibrator S/N for observation {}'
              .format(request['observation']))

    return fig
