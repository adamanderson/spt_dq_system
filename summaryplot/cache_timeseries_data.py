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

timenow = datetime.datetime.now()
dt = datetime.timedelta(-1*(timenow.weekday()+1))
default_mintime = timenow + dt
cache_dir_stub = 'cached_dq_plots'


P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', title='subcommands',
                      help='Function to perform. For help, call: '
                      '%(prog)s %(metavar)s -h')

# skim data mode
parser_skim = S.add_parser('skim',
                           help='Skim summary data from autoprocessing outputs '
                           'for timeseries plotting.')
parser_skim.add_argument('caldatapath', default=None,
                         help='Path to calibration data to skim.')
parser_skim.add_argument('bolodatapath', default=None,
                         help='Path to bolometer data (for bolometer properties.')

# plot data mode
parser_plot = S.add_parser('plot',
                           help='Plot data from skims of autoprocessing data.')
parser_plot.add_argument('timeinterval', default=None,
                         help='Time interval at which to generate plots. '
                         'Available options are "weekly", "monthly" or "N", '
                         'where N is generates a plot containing data from '
                         'only the most recent N days.')

# arguments common to each mode
for parser in [parser_skim, parser_plot]:
    parser.add_argument('outdir', action='store', default=None,
                        help='Path containing skimmed data and plots.')
    parser.add_argument('--min-time', action='store',
                        default=default_mintime.strftime('%Y%m%d'),
                        help='Minimum time of observations to skim. Format: '
                        'YYYYMMDD')
    parser.add_argument('--max-time', action='store',
                        default=timenow.strftime('%Y%m%d'),
                        help='Maximum time of observations to skim. Format: '
                        'YYYYMMDD')
    parser.add_argument('action', choices=['update', 'rebuild'], default=None,
                        help='Update or rebuild the data skims or plots.')

args = P0.parse_args()

# check timeinterval argument
if args.mode == 'plot' and \
   args.timeinterval != 'monthly' and \
   args.timeinterval != 'weekly':
    try:
        float(args.timeinterval)
    except:
        raise ValueError('Argument `timeinterval` is none of `monthly`, '
                         '`weekly` or a number of days.')


wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180',
              'w181', 'w187', 'w188', 'w201', 'w203',
              'w204', 'w206', 'all']

# parse times to loop over
dt_mintime = datetime.datetime(year=int(args.min_time[:4]),
                         month=int(args.min_time[4:6]),
                         day=int(args.min_time[6:8]))
dt_maxtime = datetime.datetime(year=int(args.max_time[:4]),
                         month=int(args.max_time[4:6]),
                         day=int(args.max_time[6:8]))

date_boundaries = []
next_day = dt_mintime
while next_day < dt_maxtime:
    date_boundaries.append(next_day)

    if args.mode == 'skim' or args.timeinterval == 'weekly': # weekly mode
        next_day = next_day + datetime.timedelta(days = 7 - next_day.weekday())

    elif args.timeinterval == 'monthly': # monthly mode
        try:
            next_day = datetime.datetime(year=next_day.year,
                                         month=next_day.month+1,
                                         day=next_day.day)
        except ValueError:
            next_day = datetime.datetime(year=next_day.year+1,
                                         month=1,
                                         day=next_day.day)
date_boundaries.append(dt_maxtime)


# SKIM MODE
if args.mode == 'skim':
    # manage directory structure
    os.makedirs(os.path.join(args.outdir, 'data'),
                     exist_ok=True)

    # delete the full output directory tree if we are in rebuild mode
    if args.mode == 'rebuild' and os.path.exists(datadir):
        shutil.rmtree(plotsdir)        
    datadir = os.path.join(args.outdir, 'data')
    os.makedirs(datadir, exist_ok=True)

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

    # loop over weeks
    for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
        # convert min/max times to observation IDs that we can compare with filenames
        min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
        max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))

        datafile = os.path.join(datadir, '{}_data_cache.pkl'.format(mindate.strftime('%Y%m%d')))
        if os.path.exists(datafile):
            with open(datafile, 'rb') as f:
                data = pickle.load(f)
        else:
            data = {}

        # update the data skim
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

        with open(datafile, 'wb') as f:
            pickle.dump(data, f)


# PLOT MODE
elif args.mode == 'plot':
    plotsdir = os.path.join(args.outdir, 'plots')
    datadir = os.path.join(args.outdir, 'data')

    # delete the full output directory tree if we are in rebuild mode
    if args.mode == 'rebuild' and os.path.exists(plotsdir):
        shutil.rmtree(plotsdir)        

    if args.timeinterval == 'monthly' or args.timeinterval == 'weekly':
        outdir = os.path.join(args.outdir, 'plots', args.timeinterval)
        os.makedirs(outdir, exist_ok=True)
    else: # checked arguments at top so we know this is castable to float
        outdir = os.path.join(args.outdir, 'plots',
                              'last_{}'.format(float(args.timeinterval)))
        os.makedirs(outdir, exist_ok=True)

    # check if plots are stale
    weekly_filenames = glob(os.join(datadir, '*pkl'))
    weekly_datetimes = {datetime.strptime(os.path.basename(dname).split('_')[0]): fname
                        for fname in weekly_filenames}

    for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
        # convert min/max time for this interval to obsids that we can compare
        # to data obsids
        min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
        max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))

        # load data from this date range
        data = {}
        for dt, fname in weekly_datetimes.iteritems():
            if dt > mindate and dt < maxdate:
                with open(fname, 'rb') as f:
                 nextdata = pickle.load(f)
                 
                 # restrict data to time range
                 for source in nextdata:
                     for obsid in nextdata[source]:
                         if obsid < min_obsid or obsid > max_obsid:
                             nextdata[source].pop(obsid)
                 data = {**data, **nextdata}

        # get maximum obsid that went into existing plots
        if os.path.exists(os.path.join(outdir, 'max_obsid.dat')):
            with open(os.path.join(outdir, 'max_obsid.dat'), 'r') as f:
                plot_obsid = int(f.readline())
        else:
            plot_obsid = 0

        # check whether data contains newer obsids than the plot, and if so,
        # remake the plots
        data_obsids = [obsid for source in data.keys() for obsid in data[source].keys()
                       if obsid > min_obsid and obsid < max_obsid]
        if any([plot_obsid < oid for oid in data_obsids]):
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

            # Get maximum obsid of data and save it to a file so that we can
            # detect when plots need to be updated.
            max_obsid = np.max([np.max(list(data[source].keys())) for source in data])
            with open(os.path.join(outdir, 'max_obsid.dat'), 'w') as f:
                f.write(max_obsid)

            # # Create symlink from latest data directory to current
            # # The latest data directory is not necessarily reliably identified by
            # # the modification timestamp; instead, it is defined by having the
            # # latest value of timestamp in its filename. Based on the naming scheme
            # # for directories, an alphanumeric sort should also produce
            # # chronological ordering, so we'll rely on this assumption.
            # symlinkname = '{}/current'.format(args.outdir)
            # dirnames = np.sort(glob.glob('{}/*{}'.format(args.outdir, cache_dir_stub)))
            # latest_dirname = dirnames[-1]
            # os.symlink(latest_dirname, '{}/temp'.format(args.outdir))
            # os.rename('{}/temp'.format(args.outdir), '{}/current'.format(args.outdir))
