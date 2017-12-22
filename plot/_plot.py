''' This script is run by the server and it in turn calls the corresponding
python scripts for plotting. This script can be used for testing of plotting
functions. Simply run the script from the main directory (not the plot
directory) using python and pass it the arguments described below.

The plotting scripts communicate with the server using stdout. A three
character code at the beginning of the message tells the server how to handle
it. If a plotting script wants a message to be displayed, prefix your
message with 'msg' followed by the message (no characters inbetween). Every
message is terminated by the first newline character.

Arguments: 
  plotting mode: currently individual or timeseries
  source: The source of the observation
  For individual:
    observation: the observation to be plotted
    plot_type: All the plot types to be made separated by a space
  For timeseries:
    plot_type: The type of timeseries plot to make
    observation: The list of observations included in the timeseries plot
'''

import sys
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from os import path, stat, mkdir
from importlib import import_module
import hashlib

# Argparse structure is kind of messed up.
# Use subparsers for separate arguments in each mode, and deal with variable
# argument lists in a standard way rather than splitting across multiple
# argparse arguments. [AJA]
parser = argparse.ArgumentParser(
    description='Calls plotting functions in other files for a webserver.')
parser.add_argument('mode', type=str,
    help='The plotting mode. Currently individual or timeseries.')
parser.add_argument('source', type=str,
    help='The source of the observation.')
parser.add_argument('plotpath', type=str,
    help='Path to which to save plots.')
parser.add_argument('caldatapath', type=str,
    help='Path to autoprocessed calibration data.')
parser.add_argument('bolodatapath', type=str,
    help='Path to bolometer data.')
parser.add_argument('single', type=str,
    help='''If individual, then the observation id. 
    If timeseries, then the name of the plot.''')
parser.add_argument('multi', type=str, nargs=argparse.REMAINDER,
    help='''If individual, then the name(s) of the plots to be created. 
    If timeseries, then the observations to be used in the plot.''')

args = parser.parse_args()
if (len(args.multi) == 0):
  print('_plot.py: error: ' +
      'the following arguments are required: mode, source, single, multi',
      flush=True)
  exit(1)

# send error message to user through stdout and send error to server
def err_handler(err, request, msg=None):
  # check if error in plotting file
  if (msg == None):
    msg = 'Error making plot {} {} {}. Check the log for details.'.format(
        request['source'], request['observation'], request['plot_type'])

  print('msg' + msg, flush=True)

  # about to leave so this plot is done
  send_plot_done()
  # server handles all logging so just raise the error
  raise err

def send_plot_done():
  print('plt', flush=True)

def send_filename(filename):
  print('fln' + filename, flush=True)

def plot(request):
  try:
    if type(request['observation']) == list:
      obs_string = hashlib.sha512(' '.join(request['observation']).encode('utf-8')).hexdigest()
    elif type(request['observation']) == str:
      obs_string = request['observation']
    plot_file = '{}/{}{}{}.png'.format(request['plotpath'],
                                       request['source'],
                                       obs_string,
                                       request['plot_type'])
    print('msg {}'.format(plot_file))
    plot_basename = path.basename(plot_file)


    # check if plot already exists
    # remove cached plots if plotting file has changed
    if (path.isfile(plot_file)):
      plot_time = stat(plot_file).st_mtime
      file_time = stat('./plot/' + request['plot_type'] + '.py').st_mtime
      if (plot_time > file_time):
        send_filename(plot_basename)
        send_plot_done()
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
        ' ' + request['observation'] + ' ' + request['plot_type']
        + '''. Likely that requested data doesn't exist in the data file.
        Check the log for more details.''')
  except Exception as e:
    err_handler(e, request)

  # check if an error occurred and the user provided a string
  if (isinstance(plot, str)):
    err_handler(Exception(plot), request,
        'Error making plot {} {} {}. {}'.format(request['source'],
                                                obs_string,
                                                request['plot_type'],
                                                plot))

  # handle if plot is not a figure
  if (isinstance(plot, plt.Figure)):
    # check if /tmp/spt_dq directory exists and make it if not
    if (not path.isdir(args.plotpath)):
      mkdir(args.plotpath)
    plot.savefig(plot_file, dpi=200)
  else:
    err_handler(TypeError('Return object from plot function '
        + request['plot_type']
        + ' is not a matplotlib figure.'), request, 'Error making plot '
        + request['plot_type']
        + '. Plot file formatted incorrectly. Check the log for details.')

  # write plot file name to stdout so that nodejs server can send it to client
  # we wait until the end to send the filename so that if an error occurs
  # the website doesn't display a blank image
  send_filename(plot_basename)

  # signal the loading screen that the plot is done
  send_plot_done();

def individual():
  print('plotting!')
  # load arguments into a dict to be passed to plotting functions
  request = {'source': args.source,
             'observation': args.single,
             'plotpath': args.plotpath,
             'caldatapath': args.caldatapath,
             'bolodatapath': args.bolodatapath}

  for plot_type in args.multi:
    try:
      request['plot_type'] = plot_type
      plot(request);
    # errors are handled by the plot method so just print their message
    # to stderr and ignore them
    except Exception as e:
      print(e, file=sys.stderr)

def timeseries():
  # load arguments into a dict to be passed to plotting functions
  request = {'source': args.source,
             'plot_type': args.single,
             'observation': args.multi,
             'plotpath': args.plotpath,
             'caldatapath': args.caldatapath, 
             'bolodatapath': args.bolodatapath}
  print(request)
  plot(request);

#start process
if __name__ == '__main__':
  if (args.mode == 'individual'):
    individual()
  elif (args.mode == 'timeseries'):
    timeseries()
