import math
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def RCW38FluxCalibration(request):
  try:
    data = [frame for frame in core.G3File('/spt/user/production/calibration/' + request['source'] + '/' + request['observation'] + '.g3')]
  except RuntimeError:
    return "Could not find data file."

  x = []
  try:
    for key, item in data[0]['RCW38FluxCalibration'].items():
      if (math.isnan(item)):
        continue
      x.append(item)
  except KeyError:
    return "RCW38 flux calibration data does not exist for this observation."

  std = np.std(x)
  mean = np.mean(x)
  y = np.asarray(x)
  for i in range(len(y)):
    if (y[i] > mean + 2.5 * std or y[i] < mean - 2.5 * std):
      x.remove(y[i])
  fig = plt.figure()
  n, bins, patches = plt.hist(x, bins=40)

  plt.xlabel('RCW38FluxCalibration')
  plt.title('RCW38FluxCalibration of ' + request['source'] + ' at time ' + request['observation'])
  return fig
