import sys
import matplotlib.pyplot as plt
from os import listdir, path
from importlib import import_module

# send message to user through stdout and send error to server
def err_handler(err, plot_type, msg=None):
  # check if error in plotting file
  if (msg == None):
    msg = 'Error making plot {}. Check the log for details.'.format(plot_type)

  sys.stdout.write(msg)
  sys.stdout.flush()
  # server handles all logging so just raise the error
  raise err

def main():
  # load arguments into a dict to be passed to plotting functions
  request = {'path': sys.argv[1], 'source': sys.argv[2],
      'observation': sys.argv[3]}

  for plot_type in sys.argv[4:]:
    try:
      request['type'] = plot_type
      plot_file = '/tmp/spt_dq/{}{}{}.png'.format(request['source'],
          request['observation'], plot_type)

      # check if plot already exists
      # TODO: remove cached plots if plotting file changes
      if (path.isfile(plot_file)):
          continue
      # load python file and function
      # only import from plot package so that the user passed plot_type
      # cannot import arbitrary modules
      module = import_module('plot.' + plot_type)
    except Exception as e:
      err_handler(e, plot_type)
    try:
      func = getattr(module, plot_type)
    except AttributeError as e:
      err_handler(e, plot_type, 'Error. Plotting file for ' + plot_type
          + ''' formatted incorrectly. There needs to be a
          function with the same name as the file.''')
    # call plotting function
    try:
      plot = func(request)
    except KeyError as e:
      # special case for key error because it likely means a key
      # doesn't exist in a data file.
      err_handler(e, plot_type, 'Error making plot ' + plot_type
          + '''. Likely that requested data doesn't exist in the data file.
          Check the log for more details.''')
    except Exception as e:
      err_handler(e, plot_type)

    # handle if plot is not a figure
    if (isinstance(plot, plt.Figure)):
      plot.savefig(plot_file)
    else:
      err_handler(TypeError('Return object from plot function ' + plot_type
          + ' is not a matplotlib figure.'), plot_type, 'Error making plot '
          + plot_type
          + '. Plot file formatted incorrectly. Check the log for details.')

#start process
if __name__ == '__main__':
    main()
