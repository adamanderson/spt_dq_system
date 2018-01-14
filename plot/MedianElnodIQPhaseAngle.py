import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from spt3g import core, calibration
from spt3g.std_processing import obsid_to_g3time
import datetime as dt

# plot median IQ phase angle for a list of observations; this is a timeseries plot
def MedianElnodIQPhaseAngle(request):
    bands = [90, 150, 220]
    
    median_elnodphase = {90: [], 150: [], 220: []}
    
    for obsid in request['observation']:
        try:
            data = [fr for fr in core.G3File('{}/calibration/{}/{}.g3' \
                                             .format(request['caldatapath'],
                                                     request['source'],
                                                     obsid))]
            boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          obsid))] \
                                                  [0]["NominalBolometerProperties"]
        except RuntimeError:
            return "Could not find data files."

        for band in bands:
            try:
                phase_data = np.array([180/np.pi * np.arctan(data[0]['ElnodEigenvalueDominantVectorQ'][bolo] / 
                                                             data[0]['ElnodEigenvalueDominantVectorI'][bolo])
                                for bolo in data[0]['ElnodEigenvalueDominantVectorQ'].keys() \
                                if boloprops[bolo].band / core.G3Units.GHz == band and \
                                       data[0]['ElnodEigenvalueDominantVectorI'][bolo] != 0])
                
                # cut bolometers with very small or zero phase angle 
                phase_data = phase_data[np.isfinite(phase_data)][np.abs(phase_data[np.isfinite(phase_data)])>1e-4]

            except KeyError:
                return "ElnodEigenvalueDominantVectorQ or *I does not exist for this observation."
            median_elnodphase[band].append(np.median(phase_data[np.isfinite(phase_data)]))

    timestamps = np.sort([obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                        for obsid in request['observation']])
    dts = [dt.datetime.fromtimestamp(ts) for ts in timestamps]
    datenums = mdates.date2num(dts)
    fig, ax = plt.subplots(1)
    for band in bands:
        ax.plot(datenums, median_elnodphase[band], 'o-', 
                label='{} GHz'.format(band))
    xfmt = mdates.DateFormatter('%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    plt.xticks(rotation=25)
    plt.legend()
    plt.xlabel('observation time')
    plt.ylabel('phase angle [deg]')
    plt.title('Elnod phase angle')
    plt.tight_layout()

    return fig
