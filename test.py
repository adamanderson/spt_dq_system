import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from spt3g import core, calibration

# makes a plot of the nominal offset given the date
def nominal_plot(path, target, obs):
    data = [frame for frame in core.G3File(path + 'nominal_online_cal.g3')]

    x = [item.x_offset for key, item in data[0]['NominalBolometerProperties'].items()]
    y = [item.y_offset for key, item in data[0]['NominalBolometerProperties'].items()]


    plt.scatter(x, y)
    plt.xlabel('x Offset')
    plt.ylabel('y Offset')
    plt.title('Nominal Bolometer Positions of ' + target[:-1] + ' at time ' + obs)
    plt.savefig('plot.png')
    plt.close()

# makes a plot of the offline offset given the date
def offline_plot(path, target, obs):
    data = [frame for frame in core.G3File('/spt/user/production/calibration/' + target + '/' + obs + '.g3')]

    x = []
    y = []
    for key, item in data[0]['PointingOffsetX'].items():
        x_off = data[0]['PointingOffsetX'][key]
        y_off = data[0]['PointingOffsetY'][key]
        if (x_off < 0.02 and x_off > -0.02 and y_off > -0.02 and y_off < 0.02):
            x.append(x_off)
            y.append(y_off)


    plt.scatter(x, y)
    plt.xlabel('x Offset')
    plt.ylabel('y Offset')
    plt.title('Offline Bolometer Positions of ' + target[:-1] + ' at time ' + obs)
    plt.savefig('plot.png')
    plt.close()

#Read data from stdin
def read_in():
    lines = sys.stdin.readlines()
    #Since our input would only be having one line, parse our JSON data from that
    return json.loads(lines[0])

def main():
    #get our data as an array from read_in()
    #lines = read_in()

    #make plot and print plot name without any newline
    if (sys.argv[1] == 'nominal'):
      nominal_plot(sys.argv[2], sys.argv[3], sys.argv[4]);
    elif (sys.argv[1] == 'offline'):
      offline_plot(sys.argv[2], sys.argv[3], sys.argv[4]);
    sys.stdout.write('plot.png')

#start process
if __name__ == '__main__':
    main()
