import matplotlib
matplotlib.use('Agg')

from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
import numpy as np
import pickle as pickle
import argparse as ap
import glob
import datetime
import os.path
import shutil
from calibrator import *
from elnod import *
from mat5a import *
from noise import *
from rcw38 import *


P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                          help='Function to perform. For help, call: '
                          '%(prog)s %(metavar)s -h')

timenow = datetime.datetime.now()
dt = datetime.timedelta(-1*(timenow.weekday()+1))
default_mintime = timenow + dt
cache_dir_stub = 'cached_dq_plots'

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
S1.add_argument('--new-plots', action='store_true', default=False,
                help='Force regeneration of new plots, even if data has not changed.')
args = P0.parse_args()


# functions that define splits
def select_band(boloprops, bolo, band):
    try:
        return boloprops[bolo].band / core.G3Units.GHz == band
    except:
        return False

def select_wafer(boloprops, bolo, wafer):
    if wafer == 'all':
        return True
    else:
        try:
            return boloprops[bolo].wafer_id == wafer
        except:
            return False

selector_dict = {('w172', 90): (select_wafer, select_band),
                 ('w172', 150): (select_wafer, select_band),
                 ('w172', 220): (select_wafer, select_band),
                 ('w174', 90): (select_wafer, select_band),
                 ('w174', 150): (select_wafer, select_band),
                 ('w174', 220): (select_wafer, select_band),
                 ('w176', 90): (select_wafer, select_band),
                 ('w176', 150): (select_wafer, select_band),
                 ('w176', 220): (select_wafer, select_band),
                 ('w177', 90): (select_wafer, select_band),
                 ('w177', 150): (select_wafer, select_band),
                 ('w177', 220): (select_wafer, select_band),
                 ('w180', 90): (select_wafer, select_band),
                 ('w180', 150): (select_wafer, select_band),
                 ('w180', 220): (select_wafer, select_band),
                 ('w181', 90): (select_wafer, select_band),
                 ('w181', 150): (select_wafer, select_band),
                 ('w181', 220): (select_wafer, select_band),
                 ('w187', 90): (select_wafer, select_band),
                 ('w187', 150): (select_wafer, select_band),
                 ('w187', 220): (select_wafer, select_band),
                 ('w188', 90): (select_wafer, select_band),
                 ('w188', 150): (select_wafer, select_band),
                 ('w188', 220): (select_wafer, select_band),
                 ('w201', 90): (select_wafer, select_band),
                 ('w201', 150): (select_wafer, select_band),
                 ('w201', 220): (select_wafer, select_band),
                 ('w203', 90): (select_wafer, select_band),
                 ('w203', 150): (select_wafer, select_band),
                 ('w203', 220): (select_wafer, select_band),
                 ('w204', 90): (select_wafer, select_band),
                 ('w204', 150): (select_wafer, select_band),
                 ('w204', 220): (select_wafer, select_band),
                 ('w206', 90): (select_wafer, select_band),
                 ('w206', 150): (select_wafer, select_band),
                 ('w206', 220): (select_wafer, select_band),
                 ('all', 90): (select_wafer, select_band),
                 ('all', 150): (select_wafer, select_band),
                 ('all', 220): (select_wafer, select_band)}
wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180',
              'w181', 'w187', 'w188', 'w201', 'w203',
              'w204', 'w206', 'all']
    
function_dict = {'RCW38':             {'RCW38SkyTransmission': rcw38_sky_transmission},
                 'RCW38-pixelraster': {'MedianRCW38FluxCalibration': median_rcw38_fluxcal,
                                       'MedianRCW38IntegralFlux': median_rcw38_intflux},
                 'MAT5A':             {'MAT5ASkyTransmission': mat5a_sky_transmission},
                 'MAT5A-pixelraster': {'MedianMAT5AFluxCalibration': median_mat5a_fluxcal,
                                       'MedianMAT5AIntegralFlux': median_mat5a_intflux},
                 'calibrator':        {'MedianCalSN_4Hz': median_cal_sn_4Hz,
                                       'MedianCalResponse_4Hz': median_cal_response_4Hz,
                                       'AliveBolosCal_4Hz': alive_bolos_cal_4Hz},
                 'elnod':             {'MedianElnodIQPhaseAngle': median_elnod_iq_phase_angle,
                                       'MedianElnodSNSlopes': median_elnod_sn_slope,
                                       'AliveBolosElnod': alive_bolos_elnod},
                 'noise':             {'NEI_0.1Hz_to_0.5Hz': median_nei_01Hz_to_05Hz,
                                       'NEI_1.0Hz_to_2.0Hz': median_nei_1Hz_to_2Hz,
                                       'NEI_3.0Hz_to_5.0Hz': median_nei_3Hz_to_5Hz,
                                       'NEI_10.0Hz_to_15.0Hz': median_nei_10Hz_to_15Hz,
                                       'NEP_0.1Hz_to_0.5Hz': median_nep_01Hz_to_05Hz,
                                       'NEP_1.0Hz_to_2.0Hz': median_nep_1Hz_to_2Hz,
                                       'NEP_3.0Hz_to_5.0Hz': median_nep_3Hz_to_5Hz,
                                       'NEP_10.0Hz_to_15.0Hz': median_nep_10Hz_to_15Hz,
                                       'NET_0.1Hz_to_0.5Hz': median_net_01Hz_to_05Hz,
                                       'NET_1.0Hz_to_2.0Hz': median_net_1Hz_to_2Hz,
                                       'NET_3.0Hz_to_5.0Hz': median_net_3Hz_to_5Hz,
                                       'NET_10.0Hz_to_15.0Hz': median_net_10Hz_to_15Hz}}

                 
# loop over data by week
dt_mintime = datetime.datetime(year=int(args.min_time[:4]),
                         month=int(args.min_time[4:6]),
                         day=int(args.min_time[6:8]))
dt_maxtime = datetime.datetime(year=int(args.max_time[:4]),
                         month=int(args.max_time[4:6]),
                         day=int(args.max_time[6:8]))
d_to_next_week = datetime.timedelta(days = 7 - dt_mintime.weekday())
date_boundaries = [dt_mintime]
next_day = dt_mintime + d_to_next_week
while next_day < dt_maxtime:
    date_boundaries.append(next_day)
    next_day = next_day + datetime.timedelta(days=7)
date_boundaries.append(dt_maxtime)

# delete the full output directory tree if we are in rebuild mode
if args.mode == 'rebuild' and os.path.exists(args.outdir):
    shutil.rmtree('{}'.format(args.outdir))
    os.mkdir('{}'.format(args.outdir))

for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
    print(mindate)

    # convert min/max times in observation IDs that we can compare with filenames
    min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
    max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))
    
    # make the subdirectory based on time
    outdir = '{}/{}_{}'.format(args.outdir, mindate.strftime('%Y%m%d'), cache_dir_stub)
    if args.mode == 'rebuild':
        os.mkdir(outdir)
        data = {}
    if args.mode == 'update':
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if os.path.exists(os.path.join(outdir, 'data_cache.pkl')):
            data = pickle.load(open(os.path.join(outdir, 'data_cache.pkl'), 'rb'))
        else:
            data = {}

    was_data_updated = False

    for source, quantities in function_dict.items():
        calfiles = glob.glob('{}/calibration/{}/*g3'.format(args.caldatapath, source))
        files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                              int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]

        print('Analyzing source: {}'.format(source))

        if source not in data.keys():
            data[source] = {}
        for fname in files_to_parse:
            obsid = os.path.splitext(os.path.basename(fname))[0]
            print('observation: {}'.format(obsid))

            if obsid not in data[source].keys() or \
                    data[source][obsid]['timestamp'] != os.path.getctime(fname):
                was_data_updated = True
                data[source][obsid] = {'timestamp': os.path.getctime(fname)}
                d = [fr for fr in core.G3File(fname)]

                boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                          .format(args.bolodatapath, \
                                                                      source, \
                                                                      obsid))] \
                                                                      [0]["NominalBolometerProperties"]

                for quantity_name in function_dict[source]:
                    func_result = function_dict[source][quantity_name](d[0], boloprops, selector_dict)
                    if func_result:
                        data[source][obsid][quantity_name] = func_result

    with open('{}/data_cache.pkl'.format(outdir), 'wb') as f:
        pickle.dump(data, f)

    # only update figures if the underlying data actually changed.
    if was_data_updated or (args.mode == 'update' and args.new_plots):
        # create the plots
        plot_median_cal_sn_4Hz(data, wafer_list, outdir)
        plot_median_cal_response_4Hz(data, wafer_list, outdir)
        plot_alive_bolos_cal_4Hz(data, wafer_list, outdir)
        plot_median_elnod_sn(data, wafer_list, outdir)
        plot_median_elnod_iq_phase(data, wafer_list, outdir)
        plot_alive_bolos_elnod(data, wafer_list, outdir)
        plot_median_rcw38_fluxcal(data, wafer_list, outdir)
        plot_median_rcw38_intflux(data, wafer_list, outdir)
        plot_rcw38_sky_transmission(data, wafer_list, outdir)
        plot_median_mat5a_fluxcal(data, wafer_list, outdir)
        plot_median_mat5a_intflux(data, wafer_list, outdir)
        plot_mat5a_sky_transmission(data, wafer_list, outdir)
        plot_median_noise(data, 'NEI_0.1Hz_to_0.5Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEI_1.0Hz_to_2.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEI_3.0Hz_to_5.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEI_10.0Hz_to_15.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEP_0.1Hz_to_0.5Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEP_1.0Hz_to_2.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEP_3.0Hz_to_5.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NEP_10.0Hz_to_15.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NET_0.1Hz_to_0.5Hz', wafer_list, outdir)
        plot_median_noise(data, 'NET_1.0Hz_to_2.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NET_3.0Hz_to_5.0Hz', wafer_list, outdir)
        plot_median_noise(data, 'NET_10.0Hz_to_15.0Hz', wafer_list, outdir)

        # Create symlink from latest data directory to current
        # The latest data directory is not necessarily reliably identified by
        # the modification timestamp; instead, it is defined by having the
        # latest value of timestamp in its filename. Based on the naming scheme
        # for directories, an alphanumeric sort should also produce
        # chronological ordering, so we'll rely on this assumption.
        symlinkname = '{}/current'.format(args.outdir)
        dirnames = np.sort(glob.glob('{}/*{}'.format(args.outdir, cache_dir_stub)))
        latest_dirname = dirnames[-1]
        os.symlink(latest_dirname, '{}/temp'.format(args.outdir))
        os.rename('{}/temp'.format(args.outdir), '{}/current'.format(args.outdir))
