import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def NET3to5HzHistogram(request):
    try:
        data = [fr for fr in core.G3File('{}/calibration/{}/{}.g3' \
                                             .format(request['caldatapath'],
                                                     request['source'],
                                                     request['observation']))]
        boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          request['observation']))] \
                                                  [0]["NominalBolometerProperties"]
        boloprops_bolos = list(boloprops.keys())
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['NET_3.0Hz_to_5.0Hz'][bolo] \
                                       for bolo in data[0]['NET_3.0Hz_to_5.0Hz'].keys() \
                                       if bolo in boloprops_bolos and \
                                       boloprops[bolo].band / core.G3Units.GHz == band])
            cal_dict[band] = cal_dict[band] / (1e-6*core.G3Units.K * np.sqrt(core.G3Units.sec))
    except KeyError:
        return "NET_3.0Hz_to_5.0Hz does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(0, 2000, 100),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.grid()
    plt.legend()
    plt.xlim([0,2000])

    plt.xlabel('NET [uK rtsec]')
    plt.title('NET in 3-5Hz for observation {}\n'.format(request['observation']))
    plt.grid()
    plt.tight_layout()
    return fig
