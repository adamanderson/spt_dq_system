import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from spt3g import core, calibration
from spt3g.std_processing import obsid_to_g3time
import datetime as dt
import pickle

# plot median calibrator S/N for a list of observations; this is a timeseries plot
def MedianCalSN(request):
    with open('skims/20180115_skim.pkl', 'rb') as f:
        data = pickle.load(f)

    waferlist = [selection_id for selection_id in \
                 data['calibrator'][list(data['calibrator'].keys())[0]]['MedianCalSN'].keys()
                 if 'w' in str(selection_id)]

    median_calSN = {wafer: [data['calibrator'][obsid]['MedianCalSN'][wafer] \
                                for obsid in request['observation'] \
                                if obsid in data['calibrator'].keys()] \
                        for wafer in waferlist}
    timestamps = np.sort([obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in request['observation'] \
                              if obsid in data['calibrator'].keys()])

    dts = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
    datenums = mdates.date2num(dts)

    fig, ax = plt.subplots(1)
    for wafer in waferlist:
        ax.plot(datenums, median_calSN[wafer], 'o-', label=wafer)

    xfmt = mdates.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    plt.xticks(rotation=25)
    plt.legend()
    plt.xlabel('observation time')
    plt.ylabel('median calibrator S/N')
    plt.title('Calibrator S/N')
    plt.tight_layout()

    return fig
