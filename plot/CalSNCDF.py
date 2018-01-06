import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def CalSNCDF(request):
    try:
        data = [fr for fr in core.G3File('{}/calibration/{}/{}.g3' \
                                             .format(request['caldatapath'],
                                                     request['source'],
                                                     request['observation']))]
        boloprops = [fr for fr in core.G3File('{}/downsampled/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          request['observation']))] \
                                                  [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['CalibratorResponseSN'][bolo] \
                                           for bolo in data[0]['CalibratorResponseSN'].keys() \
                                           if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "CalibratorResponseSN does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        histvals, edges = np.histogram(cal_dict[band][np.isfinite(cal_dict[band])],
                                       bins=np.linspace(0,400,50))
        ecdf_unnormed = np.cumsum(histvals)
        ecdf_unnormed = np.hstack([0, ecdf_unnormed])
        inv_ecdf_unnormed = np.max(ecdf_unnormed) - ecdf_unnormed
        plt.step(edges, inv_ecdf_unnormed, label='{} GHz'.format(band))
    #plt.yscale('log')
    plt.legend()

    plt.xlabel('calibrator S/N')
    plt.title('cumulative bolometers above a given calibrator S/N\n'
              'for observation {}\n'
              'chop freq. = {:.1f} Hz'.format(request['observation'],
                                              data[0]["CalibratorResponseFrequency"] / core.G3Units.Hz))
    plt.tight_layout()
    return fig
