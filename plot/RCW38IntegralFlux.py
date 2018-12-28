import math
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def RCW38IntegralFlux(request):
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
      cal_dict[band] = np.array([data[0]['RCW38IntegralFlux'][bolo] \
                                 for bolo in data[0]['RCW38IntegralFlux'].keys() \
                                 if bolo in boloprops_bolos and \
                                 boloprops[bolo].band / core.G3Units.GHz == band])
  except KeyError:
    return "RCW38 integral flux does not exist for this observation."

  fig = plt.figure()
  for band in bands:
    x = np.array(cal_dict[band])
    x = x[np.isfinite(x)]

    plt.hist(x,
             bins=np.linspace(2e-7, 6e-7, 50),
             label='{} GHz'.format(band),
             histtype='step')
  plt.legend(loc='upper left')

  plt.xlabel('RCW38IntegralFlux')
  plt.title('RCW38IntegralFlux of ' + request['source'] + ' at time ' + request['observation'])
  return fig
