import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the nominal offset given the date
def NominalPointingOffsets(request):
  try:
    data = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                  .format(request['bolodatapath'],
                                                          request['source'],
                                                          request['observation']))] \
                                                  [0]["NominalBolometerProperties"]
  except RuntimeError:
    return "Could not find data file."

  x = []
  y = []
  for key, item in data.items():
    x_off = item.x_offset
    y_off = item.y_offset
    if (x_off < 0.02 and x_off > -0.02 and y_off > -0.02 and y_off < 0.02):
      x.append(x_off / core.G3Units.deg)
      y.append(y_off / core.G3Units.deg)


  fig = plt.figure()
  plt.scatter(x, y, s=10)
  plt.xlabel('x Offset [deg]')
  plt.ylabel('y Offset [deg]')
  plt.title('Nominal Bolometer Positions of ' + request['source'] + ' at time ' + request['observation'])
  plt.axis(np.array([-0.02, 0.02, -0.02, 0.02]) / core.G3Units.deg)
  return fig
