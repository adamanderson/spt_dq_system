import sys
import matplotlib.pyplot as plt
from os import listdir
from importlib import import_module

def main():
  # load arguments into a dict to be passed to plotting functions
  request = {'path': sys.argv[1], 'source': sys.argv[2], 'observation': sys.argv[3]}

  for plot_type in sys.argv[4:]:
    request['type'] = plot_type
    # load python file and function
    # only import from plot package so that the user passed plot_type
    # cannot import arbitrary modules
    module = import_module('plot.' + plot_type)
    func = getattr(module, request['type'] + '_plot')
    # call plotting function and handle errors
    try:
        plot = func(request)
    except KeyError:
        sys.exit(3)
    except RuntimeError as e:
        # if not a cannot find file error, raise it
        if (str(e).find('Could not find file') == -1):
            raise e
        sys.exit(4)
    finally:
        sys.stdout.write(plot_type)
        sys.stdout.flush()

    plot.savefig('/tmp/spt_dq/' + request['source'] + request['observation'] + plot_type + '.png')

#start process
if __name__ == '__main__':
    main()
