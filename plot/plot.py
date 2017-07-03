import sys
import matplotlib.pyplot as plt
from os import listdir

def main():
  # load arguments into a dict to be passed to plotting functions
  request = {'path': sys.argv[1], 'source': sys.argv[2], 'observation': sys.argv[3]}

  for plot_type in sys.argv[4:]:
    request['type'] = plot_type
    # load python file and function
    module = __import__(request['type'])
    func = getattr(module, request['type'] + '_plot')
    # call plotting function print plot name without any newline
    plot = func(request)
    plot.savefig('/tmp/spt_dq/' + request['source'] + request['observation'] + plot_type + '.png')

#start process
if __name__ == '__main__':
    main()
