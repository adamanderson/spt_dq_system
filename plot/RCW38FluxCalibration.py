import math
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def RCW38FluxCalibration(request):
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
      cal_dict[band] = np.array([data[0]['RCW38FluxCalibration'][bolo] \
                          for bolo in data[0]['RCW38FluxCalibration'].keys() \
                          if boloprops[bolo].band / core.G3Units.GHz == band])
  except KeyError:
    return "RCW38 flux calibration does not exist for this observation."

  fig = plt.figure()
  for band in bands:
    x = np.array(cal_dict[band])
    x = x[np.isfinite(x)]

    plt.hist(x,
             bins=np.linspace(-90, 10, 50),
             label='{} GHz'.format(band),
             histtype='step')
  plt.legend(loc='upper left')

  plt.xlabel('RCW38FluxCalibration')
  plt.title('RCW38FluxCalibration of ' + request['source'] + ' at time ' + request['observation'])
  return fig
