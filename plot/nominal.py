import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the nominal offset given the date
def nominal(request):
  try:
    data = [frame for frame in core.G3File('/spt/data/bolodata/downsampled/'
      + request['source'] + '/' + request['observation']
      + '/nominal_online_cal.g3')]
  except RuntimeError:
    return "Could not find data file."

  x = []
  y = []
  for key, item in data[0]['NominalBolometerProperties'].items():
    x_off = item.x_offset
    y_off = item.y_offset
    if (x_off < 0.02 and x_off > -0.02 and y_off > -0.02 and y_off < 0.02):
      x.append(x_off)
      y.append(y_off)


  fig = plt.figure()
  plt.scatter(x, y)
  plt.xlabel('x Offset')
  plt.ylabel('y Offset')
  plt.title('Nominal Bolometer Positions of ' + request['source'] + ' at time ' + request['observation'])
  return fig
