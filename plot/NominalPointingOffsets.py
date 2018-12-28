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
  bands = []
  pols = []
  for key, item in data.items():
    x_off = item.x_offset
    y_off = item.y_offset
    if (x_off < 0.02 and x_off > -0.02 and y_off > -0.02 and y_off < 0.02):
      x.append(x_off / core.G3Units.deg)
      y.append(y_off / core.G3Units.deg)
      bands.append(item.band / core.G3Units.GHz)
      pols.append(item.physical_name.split('.')[-1])
  x = np.array(x)
  y = np.array(y)
  bands = np.array(bands)
  pols = np.array(pols)

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
  plt.suptitle('{}: nominal bolometer offsets'.format( request['observation']))
  return fig
