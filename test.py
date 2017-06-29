import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from spt3g import core, calibration

# makes a plot of the nominal offset given the date
def make_plot(path, target, obs):
    data = [frame for frame in core.G3File(path + 'nominal_online_cal.g3')]

    x = [item.x_offset for key, item in data[0]['NominalBolometerProperties'].items()]
    y = [item.y_offset for key, item in data[0]['NominalBolometerProperties'].items()]


    plt.scatter(x, y)
    plt.xlabel('x Offset')
    plt.ylabel('y Offset')
    plt.title('Nominal Bolometer Positions of ' + target[:-1] + ' at time ' + obs)
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
    make_plot(sys.argv[1], sys.argv[2], sys.argv[3])
    sys.stdout.write('plot.png')

#start process
if __name__ == '__main__':
    main()
