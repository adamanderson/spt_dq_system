import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def ElnodSigmaSlopesHistogram(request):
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
            cal_dict[band] = np.array([data[0]['ElnodSigmaSlopes'][bolo] \
                                           for bolo in data[0]['ElnodSigmaSlopes'].keys() \
                                           if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "ElnodSigmaSlopes does not exist for this observation."

    el = np.median([pointingframe["OnlineBoresightEl"][i]*180/np.pi for i in range(len(pointingframe["OnlineBoresightEl"]))])

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(0,500,101),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlim([0, 500])
    plt.xlabel('elnod sigma slopes')
    plt.title('Elnod sigma slope for observation {} at el={:.1f} deg'.format(request['observation'], el))
    plt.tight_layout()
    return fig
