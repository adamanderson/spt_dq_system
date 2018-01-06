import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def CalHistogram(request):
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
            cal_dict[band] = np.array([1e15 * data[0]['CalibratorResponse'][bolo] \
                                  for bolo in data[0]['CalibratorResponse'].keys() \
                                  if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "CalibratorResponse does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(-1, 9, 50),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlabel('calibrator response [fW]')
    plt.title('Calibrator response for observation {}\n'
              'chop freq. = {:.1f} Hz'.format(request['observation'],
                                              data[0]["CalibratorResponseFrequency"] / core.G3Units.Hz))
    plt.tight_layout()
    return fig
