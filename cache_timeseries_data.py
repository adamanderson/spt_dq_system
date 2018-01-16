from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
from spt3g.std_processing import obsid_to_g3time
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
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

timenow = datetime.datetime.now()
default_mintime = datetime.datetime(timenow.year, timenow.month, timenow.day - (timenow.weekday()+3))

S0 = S.add_parser('rebuild', help='Rebuild the pickle files from scratch.',
                  formatter_class=ap.ArgumentDefaultsHelpFormatter)
S0.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
S0.add_argument('bolodatapath', action='store', default=None,
                help='Path to bolometer data (for bolometer properties.')
S0.add_argument('outdir', action='store', default=None,
                help='Path in which to store output data.')
S0.add_argument('--min-time', action='store', default=default_mintime.strftime('%Y%m%d'),
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S0.add_argument('--max-time', action='store',
                default=timenow.strftime('%Y%m%d'),
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
S1.add_argument('outdir', action='store', default=None,
                help='Path in which to store output data.')
S1.add_argument('--min-time', action='store', default=default_mintime.strftime('%Y%m%d'),
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S1.add_argument('--max-time', action='store',
                default=timenow.strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')
args = P0.parse_args()
print(args.min_time)
print(args.max_time)

outdir = '{}/{}_cached_dq_plots'.format(args.outdir,args.min_time)

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
wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180', 'w181', 'w187', 'w188', 'w201', 'w203']


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
# update the directory if in update mode and the directory already exists
if args.mode == 'update' and os.path.exists(outdir):
    with open('{}/data_cache.pkl'.format(outdir), 'rb') as f:
        data = pickle.load(f)
# otherwise, build a new directory
else:
    os.mkdir('{}'.format(outdir))
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
        
        print(data[source].keys())
        if obsid not in data[source].keys() or \
                data[source][obsid]['timestamp'] != os.path.getctime(fname):
            data[source][obsid] = {'timestamp': os.path.getctime(fname)}
            d = [fr for fr in core.G3File(fname)]
            boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                      .format(args.bolodatapath, \
                                                                  source, \
                                                                  obsid))] \
                                                                  [0]["NominalBolometerProperties"]

            for quantity_name in function_dict[source]:
                data[source][obsid][quantity_name] = \
                    function_dict[source][quantity_name](d[0], selector_dict)
            
with open('{}/data_cache.pkl'.format(outdir), 'wb') as f:
    pickle.dump(data, f)




def median_cal_sn(data):
    for wafer in wafer_list:
        obsids = [obsid for obsid in data['calibrator']]
        median_calSN = [data['calibrator'][obsid]['MedianCalSN'][wafer] for obsid in data['calibrator']]
        print(obsids)
        print(median_calSN)

        f = plt.figure()
        timestamps = np.sort([obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids])
        dts = [datetime.datetime.fromtimestamp(ts) for ts in timestamps]
        datenums = mdates.date2num(dts)

        plt.plot(datenums, median_calSN, 'o', label=wafer)

        xfmt = mdates.DateFormatter('%m-%d %H:%M')
        plt.gca().xaxis.set_major_formatter(xfmt)
        plt.xticks(rotation=25)
        plt.legend()
        plt.xlabel('observation time')
        plt.ylabel('median calibrator S/N')
        plt.title('Calibrator S/N')
        plt.tight_layout()
        plt.savefig('{}/median_cal_sn_{}.png'.format(outdir, wafer))
        plt.close()

# create the plots
median_cal_sn(data)
