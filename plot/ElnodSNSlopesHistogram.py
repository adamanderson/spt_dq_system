import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def ElnodSNSlopesHistogram(request):
    try:
        data = [fr for fr in core.G3File('/spt/user/production/calibration/{}/{}.g3' \
                                             .format(request['source'], request['observation']))]
        boloprops = [fr for fr in core.G3File('/spt/data/bolodata/downsampled/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['source'], request['observation']))] \
                                                  [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = [data[0]['ElnodSNSlopes'][bolo] \
                                  for bolo in data[0]['ElnodSNSlopes'].keys() \
                                  if boloprops[bolo].band / core.G3Units.GHz == band]
    except KeyError:
        return "ElnodSNSlopes does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band],
                 bins=np.linspace(-100,900,50),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlim([-100, 900])
    plt.xlabel('elnod slopes S/N')
    plt.title('Elnod slope S/N for observation {}'.format(request['observation']))
    plt.tight_layout()
    return fig
