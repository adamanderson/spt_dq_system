import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def ElnodSlopesHistogram(request):
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
        bolodatafile = core.G3File('{}/{}/{}/0000.g3' \
                                       .format(request['bolodatapath'],
                                               request['source'],
                                               request['observation']))
        for i in range(3):
            pointingframe = bolodatafile.next()
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['ElnodSlopes'][bolo] \
                                       for bolo in data[0]['ElnodSlopes'].keys() \
                                       if bolo in boloprops_bolos and \
                                       boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "ElnodSlopes does not exist for this observation."

    el = np.median([pointingframe["OnlineBoresightEl"][i]*180/np.pi for i in range(len(pointingframe["OnlineBoresightEl"]))])

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(-100000,1000,101),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlim([-100000, 1000])
    plt.xlabel('elnod slope')
    plt.title('Elnod slope for observation {} at el={:.1f} deg'.format(request['observation'], el))
    plt.tight_layout()
    return fig
