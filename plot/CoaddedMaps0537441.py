import math
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration, mapmaker

# makes a plot of the offline offset given the date
def CoaddedMaps0537441(request):
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

  fig = plt.figure(figsize=(12,4))
  for jband in range(len(data)):
    plt.subplot(1,3,jband+1)
    mp = data[jband]['T']
    plt.imshow(mp, extent=np.array([ mp.alpha_center - (mp.shape[1]*mp.y_res)/2,
                                     mp.alpha_center + (mp.shape[1]*mp.y_res)/2, 
                                     mp.delta_center - (mp.shape[0]*mp.x_res)/2,
                                     mp.delta_center + (mp.shape[0]*mp.x_res)/2]) / core.G3Units.deg)
    
    plt.xlabel('RA [deg]')
    plt.ylabel('dec [deg]')
    plt.axis([mp.alpha_center / core.G3Units.deg - 0.25,
              mp.alpha_center / core.G3Units.deg + 0.25,
              mp.delta_center / core.G3Units.deg - 0.25,
              mp.delta_center / core.G3Units.deg + 0.25])
    plt.title(data[jband]['Id'])
  plt.tight_layout()

  return fig
