import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date                                                                                 
def ElnodIQPhaseAngle(request):
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
        plot_dict = {}
        for band in bands:
            plot_dict[band] = np.array([180/np.pi * np.arctan(data[0]['ElnodEigenvalueDominantVectorQ'][bolo] / data[0]['ElnodEigenvalueDominantVectorI'][bolo])
                                           for bolo in data[0]['ElnodEigenvalueDominantVectorQ'].keys() 
                                           if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "ElnodEigenvalueDominantVectorQ or *I does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(plot_dict[band][np.isfinite(plot_dict[band])],
                 bins=np.linspace(-45,45,91),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlim([-45, 45])
    plt.xlabel('phase angle [deg]')
    plt.title('Elnod phase angle for observation {}'.format(request['observation']))
    plt.tight_layout()
    return fig
