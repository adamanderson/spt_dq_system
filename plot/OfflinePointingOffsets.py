import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def OfflinePointingOffsets(request):
    try:
        data = [fr for fr in core.G3File('{}/{}/{}.g3' \
                                         .format(request['caldatapath'],
                                                 request['source'],
                                                 request['observation']))][0]
        boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                              .format(request['bolodatapath'],
                                                      request['source'],
                                                      request['observation']))] \
          [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    try:  
        x = np.array([data['PointingOffsetX'][bolo] / core.G3Units.deg
                      for bolo in data['PointingOffsetX'].keys()])
        y = np.array([data['PointingOffsetY'][bolo] / core.G3Units.deg
                      for bolo in data['PointingOffsetY'].keys()])
        bands = np.array([boloprops[bolo].band / core.G3Units.GHz
                          for bolo in data['PointingOffsetX'].keys()])
        pols = np.array([boloprops[bolo].physical_name.split('.')[-1]
                         for bolo in data['PointingOffsetX'].keys()])
    except KeyError:
        return "Offline calibration does not exist for this observation."


    fig = plt.figure(figsize=(15,5))
    for jband, band in enumerate([90, 150, 220]):
        plt.subplot(1,3,jband+1)
        plt.scatter(x[(bands==band) & (pols=='x')],
                    y[(bands==band) & (pols=='x')],
                    s=10, marker='_', linewidth=1, color='b', label='x pol')
        plt.scatter(x[(bands==band) & (pols=='y')],
                    y[(bands==band) & (pols=='y')],
                    s=10, marker='|', linewidth=1, color='r', label='y pol')
        plt.axis(np.array([-0.02, 0.02, -0.02, 0.02]) / core.G3Units.deg)
        plt.grid()
        plt.legend(prop={'size':10})
        plt.title('{} GHz'.format(band))

    plt.subplot(1,3,2)
    plt.xlabel('x Offset [deg]')
    plt.subplot(1,3,1)
    plt.ylabel('y Offset [deg]')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('{}: offline bolometer offsets'.format( request['observation']))
    return fig
