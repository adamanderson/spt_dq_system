import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from os import path, stat
from importlib import import_module

# send message to user through stdout and send error to server
def err_handler(err, request, msg=None):
  # check if error in plotting file
  if (msg == None):
    msg = 'Error making plot {} {} {}. Check the log for details.'.format(request['source'], request['observation'], request['type'])

  sys.stdout.write(msg)
  sys.stdout.flush()
  # server handles all logging so just raise the error
  raise err

def plot(request):
  try:
    plot_file = '/tmp/spt_dq/{}{}{}.png'.format(request['source'],
        request['observation'], request['plot_type'])

    # check if plot already exists
    # remove cached plots if plotting file has changed
    if (path.isfile(plot_file)):
      plot_time = stat(plot_file).st_mtime
      file_time = stat('./plot/' + request['plot_type'] + '.py').st_mtime
      if (plot_time > file_time):
        return
    # load python file and function
    # only import from plot package so that the user passed plot_type
    # cannot import arbitrary modules
    module = import_module('plot.' + request['plot_type'])
  except Exception as e:
    err_handler(e, request)
  try:
    func = getattr(module, request['plot_type'])
  except AttributeError as e:
    err_handler(e, request, 'Error. Plotting file for ' + request['plot_type']
        + ''' formatted incorrectly. There needs to be a
        function with the same name as the file.''')
  # call plotting function
  try:
    plot = func(request)
  except KeyError as e:
    # special case for key error because it likely means a key
    # doesn't exist in a data file.
    err_handler(e, request, 'Error making plot ' + request['source'] +
        ' ' + requset['observation'] + ' ' + request['plot_type']
        + '''. Likely that requested data doesn't exist in the data file.
        Check the log for more details.''')
  except Exception as e:
    err_handler(e, request)

  # check if an error occurred and the user provided a string
  if (isinstance(plot, str)):
    err_handler(Exception(plot), request,
        'Error making plot ' + request['source'] + ' ' +
        request['observation'] + ' ' + request['plot_type'] + '. ' + plot);

  # handle if plot is not a figure
  if (isinstance(plot, plt.Figure)):
    plot.savefig(plot_file, dpi=200)
  else:
    err_handler(TypeError('Return object from plot function ' + request['plot_type']
        + ' is not a matplotlib figure.'), request, 'Error making plot '
        + request['plot_type']
        + '. Plot file formatted incorrectly. Check the log for details.')

def individual():
  # load arguments into a dict to be passed to plotting functions
  request = {'source': sys.argv[2], 'observation': sys.argv[3]}

  for plot_type in sys.argv[4:]:
      request['plot_type'] = plot_type
      plot(request);

def timeseries():
  # load arguments into a dict to be passed to plotting functions
  request = {'source': sys.argv[2], 'plot_type': sys.argv[3], 'observation': ' '.join(sys.argv[4:])}
  plot(request);

#start process
if __name__ == '__main__':
  if (sys.argv[1] == 'individual'):
    individual()
  elif (sys.argv[1] == 'timeseries'):
    timeseries()
