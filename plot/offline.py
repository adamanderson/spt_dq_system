import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the offline offset given the date
def offline(request):
  try:
    data = [frame for frame in core.G3File('/spt/user/production/calibration/' + request['source'] + '/' + request['observation'] + '.g3')]
  except RuntimeError:
    return "Could not find data file."

  x = []
  y = []
  try:
    for key, item in data[0]['PointingOffsetX'].items():
      x_off = data[0]['PointingOffsetX'][key]
      y_off = data[0]['PointingOffsetY'][key]
      if (x_off < 0.02 and x_off > -0.02 and y_off > -0.02 and y_off < 0.02):
        x.append(x_off)
        y.append(y_off)
  except KeyError:
    return "Offline calibration does not exist for this observation."


  fig = plt.figure()
  plt.scatter(x, y)
  plt.xlabel('x Offset')
  plt.ylabel('y Offset')
  plt.title('Offline Bolometer Positions of ' + request['source'] + ' at time ' + request['observation'])
  return fig
