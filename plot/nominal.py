import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration

# makes a plot of the nominal offset given the date
def nominal(request):
    data = [frame for frame in core.G3File(request['path'] + 'nominal_online_cal.g3')]

    x = [item.x_offset for key, item in data[0]['NominalBolometerProperties'].items()]
    y = [item.y_offset for key, item in data[0]['NominalBolometerProperties'].items()]


    fig = plt.figure()
    plt.scatter(x, y)
    plt.xlabel('x Offset')
    plt.ylabel('y Offset')
    plt.title('Nominal Bolometer Positions of ' + request['source'] + ' at time ' + request['observation'])
    return fig
