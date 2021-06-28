import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration
from plotutils import get_elevation

# makes a plot of the offline offset given the date                                                                                 
def ElnodSNSlopesHistogram(request):
    try:
        data = [fr for fr in core.G3File('{}/{}/{}.g3' \
                                             .format(request['caldatapath'],
                                                     request['source'],
                                                     request['observation']))]
        boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          request['observation']))] \
                                                  [0]["NominalBolometerProperties"]
        boloprops_bolos = list(boloprops.keys())

        el = get_elevation('{}/{}/{}/0000.g3' \
                   .format(request['bolodatapath'],
                           request['source'],
                           request['observation']))

    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['ElnodSNSlopes'][bolo] \
                                       for bolo in data[0]['ElnodSNSlopes'].keys() \
                                       if bolo in boloprops_bolos and \
                                       boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "ElnodSNSlopes does not exist for this observation."


    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(-100,4000,101),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlim([-100, 4000])
    plt.xlabel('elnod slopes S/N')
    plt.title('Elnod slope S/N for observation {} at el={:.1f} deg'.format(request['observation'], el))
    plt.grid()
    plt.tight_layout()
    return fig
