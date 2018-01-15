import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from spt3g import core, calibration
from spt3g.std_processing import obsid_to_g3time
import datetime as dt
import pickle

# plot median IQ phase angle for a list of observations; this is a timeseries plot
def MedianElnodIQPhaseAngle(request):
    with open('skims/20180115_skim.pkl', 'rb') as f:
        data = pickle.load(f)

    waferlist = [selection_id for selection_id in \
                 data['elnod'][list(data['elnod'].keys())[0]]['MedianElnodIQPhaseAngle'].keys()
                 if 'w' in str(selection_id)]
    median_iq_phase = {wafer: [data['elnod'][obsid]['MedianElnodIQPhaseAngle'][wafer] \
                                   for obsid in request['observation'] \
                                   if obsid in data['elnod'].keys()] \
                           for wafer in waferlist}
    timestamps = np.sort([obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in request['observation'] \
                              if obsid in data['elnod'].keys()])

    dts = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
    datenums = mdates.date2num(dts)

    fig, ax = plt.subplots(1)
    for wafer in waferlist:
        print(median_iq_phase[wafer])
        ax.plot(datenums, median_iq_phase[wafer], 'o-', label=wafer)

    xfmt = mdates.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    plt.xticks(rotation=25)
    plt.legend()
    plt.xlabel('observation time')
    plt.ylabel('phase angle [deg]')
    plt.title('Elnod phase angle')
    plt.tight_layout()

    return fig
