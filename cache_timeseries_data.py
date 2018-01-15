from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
import numpy as np
import pickle as pickle
import argparse as ap
import glob
import datetime
import os.path

P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                          help='Function to perform. For help, call: '
                          '%(prog)s %(metavar)s -h')

S0 = S.add_parser('rebuild', help='Rebuild the pickle files from scratch.',
                  formatter_class=ap.ArgumentDefaultsHelpFormatter)
S0.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
S0.add_argument('bolodatapath', action='store', default=None,
                help='Path to bolometer data (for bolometer properties.')
S0.add_argument('outfilename', action='store', default=None,
                help='name of output data file to write.')
S0.add_argument('--min-time', action='store', default='20180101',
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S0.add_argument('--max-time', action='store',
                default=datetime.datetime.now().strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')

S1 = S.add_parser('update', help='Updates the pickle file data skim when the '
                  'autoprocessing file timestamp is newer than the stored '
                  'timestamp or when there is no stored data for an '
                  'autoprocessing data file.',
                  formatter_class=ap.ArgumentDefaultsHelpFormatter)
S1.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
S1.add_argument('bolodatapath', action='store', default=None,
                help='Path to bolometer data (for bolometer properties.')
S1.add_argument('infilename', action='store', default=None,
                help='Name of input file to update.')
S1.add_argument('outfilename', action='store', default=None,
                help='Name of output data file to write.')
S1.add_argument('--min-time', action='store', default='20180101',
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S1.add_argument('--max-time', action='store',
                default=datetime.datetime.now().strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')

args = P0.parse_args()


# convert min/max times in observation IDs that we can compare with filenames
min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(args.min_time)))
max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(args.max_time)))


# functions that define splits
def select_band(boloprops, bolo, band):
    return boloprops[bolo].band / core.G3Units.GHz == band

def select_wafer(boloprops, bolo, wafer):
    return boloprops[bolo].wafer_id == wafer

selector_dict = {90: select_band,
                 150: select_band,
                 220: select_band,
                 'w172': select_wafer,
                 'w174': select_wafer,
                 'w176': select_wafer,
                 'w177': select_wafer,
                 'w180': select_wafer,
                 'w181': select_wafer,
                 'w187': select_wafer,
                 'w188': select_wafer,
                 'w201': select_wafer,
                 'w203': select_wafer}


# functions that define quantities to be saved
def median_cal_sn(frame, selector_dict):
    calSN = np.array([frame['CalibratorResponseSN'][bolo] 
                      for bolo in frame['CalibratorResponseSN'].keys()])
    calSN_on_selection = {}
    for select_value, f_select in selector_dict.items():
        selection = np.array([f_select(boloprops, bolo, select_value)
                              for bolo in frame['CalibratorResponseSN'].keys()])
        calSN_on_selection[select_value] = np.median(calSN[selection]
                                                     [np.isfinite(calSN[selection])])
    return calSN_on_selection

def median_elnod_iq_phase_angle(frame, selector_dict):
    phase_data = np.array([180/np.pi * np.arctan(frame['ElnodEigenvalueDominantVectorQ'][bolo] / \
                                                     frame['ElnodEigenvalueDominantVectorI'][bolo]) \
                               for bolo in frame['ElnodEigenvalueDominantVectorQ'].keys() \
                               if frame['ElnodEigenvalueDominantVectorI'][bolo] != 0])
    phase_data_on_selection = {}
    for select_value, f_select in selector_dict.items():
        selection = np.array([f_select(boloprops, bolo, select_value)
                              for bolo in frame['ElnodEigenvalueDominantVectorQ'].keys()
                              if frame['ElnodEigenvalueDominantVectorI'][bolo] != 0])
        phase_data_on_selection[select_value] = np.median(phase_data[selection]
                                                          [np.isfinite(phase_data[selection])])
    return phase_data_on_selection
    
function_dict = {'calibrator': {'MedianCalSN': median_cal_sn},
                 'elnod': {'MedianElnodIQPhaseAngle': median_elnod_iq_phase_angle}}


# create the output data dictionary
print(args.mode)
if args.mode == 'update':
    with open(args.infilename, 'rb') as f:
        data = pickle.load(f)
else:
    data = {}


for source, quantities in function_dict.items():
    calfiles = glob.glob('{}/{}/*g3'.format(args.caldatapath, source))
    files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                          int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]
    
    print('Analyzing source: {}'.format(source))

    if source not in data.keys():
        data[source] = {}
    for fname in files_to_parse:
        obsid = os.path.splitext(os.path.basename(fname))[0]
        print('observation: {}'.format(obsid))

        d = [fr for fr in core.G3File(fname)]
        boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                  .format(args.bolodatapath, \
                                                          source, \
                                                          obsid))] \
                                                  [0]["NominalBolometerProperties"]
        if obsid not in data[source].keys() or \
                data[source][obsid]['timestamp'] != os.path.getctime(fname):
            data[source][obsid] = {'timestamp': os.path.getctime(fname)}
            for quantity_name in function_dict[source]:
                data[source][obsid][quantity_name] = \
                    function_dict[source][quantity_name](d[0], selector_dict)
            
with open(args.outfilename, 'wb') as f:
    pickle.dump(data, f)
