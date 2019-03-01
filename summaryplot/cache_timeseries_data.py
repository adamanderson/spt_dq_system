import matplotlib
matplotlib.use('Agg')

from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
import numpy as np
import pickle as pickle
import argparse as ap
from glob import glob
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


P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', title='subcommands',
                      help='Function to perform. For help, call: '
                      '%(prog)s %(metavar)s -h')

# skim data mode
parser_skim = S.add_parser('skim',
                           help='Skim summary data from autoprocessing outputs '
                           'for timeseries plotting.')
parser_skim.add_argument('action', choices=['update', 'rebuild'], default=None,
                         help='Update or rebuild the data skims or plots.')
parser_skim.add_argument('outdir', action='store', default=None,
                         help='Path containing skimmed data and plots.')
parser_skim.add_argument('caldatapath', default=None,
                         help='Path to calibration data to skim.')
parser_skim.add_argument('bolodatapath', default=None,
                         help='Path to bolometer data (for bolometer properties.')
parser_skim.add_argument('--min-time', action='store',
                         default=default_mintime.strftime('%Y%m%d'),
                         help='Minimum time of observations to skim. Format: '
                         'YYYYMMDD')
parser_skim.add_argument('--max-time', action='store',
                         default=timenow.strftime('%Y%m%d'),
                         help='Maximum time of observations to skim. Format: '
                         'YYYYMMDD')

# plot data mode
parser_plot = S.add_parser('plot',
                           help='Plot data from skims of autoprocessing data.')
parser_plot.add_argument('action', choices=['update', 'rebuild'], default=None,
                         help='Update or rebuild the data skims or plots.')
parser_plot.add_argument('timeinterval', default=None,
                         help='Time interval at which to generate plots. '
                         'Available options are "weekly", "monthly" or "N", '
                         'where N is generates a plot containing data from '
                         'only the most recent N days.')
parser_plot.add_argument('outdir', action='store', default=None,
                         help='Path containing skimmed data and plots.')
parser_plot.add_argument('--min-time', action='store',
                         default=default_mintime.strftime('%Y%m%d'),
                         help='Minimum time of observations to skim. Format: '
                         'YYYYMMDD')
parser_plot.add_argument('--max-time', action='store',
                         default=timenow.strftime('%Y%m%d'),
                         help='Maximum time of observations to skim. Format: '
                         'YYYYMMDD')
args = P0.parse_args()

# check timeinterval argument
if args.mode == 'plot' and \
   args.timeinterval != 'monthly' and \
   args.timeinterval != 'weekly':
    try:
        int(args.timeinterval)
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

if args.mode == 'skim' or (args.timeinterval == 'weekly' or
                           args.timeinterval == 'monthly'):
    while next_day < dt_maxtime:
        date_boundaries.append(next_day)

        if args.mode == 'skim' or args.timeinterval == 'weekly': # weekly mode
            next_day = next_day + datetime.timedelta(days = 7 - next_day.weekday())

        elif args.timeinterval == 'monthly': # monthly mode
            try:
                next_day = datetime.datetime(year=next_day.year,
                                             month=next_day.month+1,
                                             day=1)
            except ValueError:
                next_day = datetime.datetime(year=next_day.year+1,
                                             month=1,
                                             day=1)
    date_boundaries.append(dt_maxtime)
else:
    date_boundaries = [dt_maxtime - datetime.timedelta(days=int(args.timeinterval)),
                       dt_maxtime]


# SKIM MODE
if args.mode == 'skim':
    # manage directory structure
    os.makedirs(os.path.join(args.outdir, 'data'),
                     exist_ok=True)

    # delete the full output directory tree if we are in rebuild mode
    datadir = os.path.join(args.outdir, 'data')
    if args.action == 'rebuild' and os.path.exists(datadir):
        shutil.rmtree(datadir)        
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
                                           'MedianRCW38IntegralFlux': median_rcw38_intflux,
                                           'AliveBolosRCW38': alive_bolos_rcw38_fluxcal},
                     'MAT5A':             {'MAT5ASkyTransmission': mat5a_sky_transmission},
                     'MAT5A-pixelraster': {'MedianMAT5AFluxCalibration': median_mat5a_fluxcal,
                                           'MedianMAT5AIntegralFlux': median_mat5a_intflux,
                                           'AliveBolosMAT5A': alive_bolos_mat5a_fluxcal},
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
    function_dict_raw = {'calibrator':    {'elevation': mean_cal_elevation}}

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
            calfiles = glob('{}/calibration/{}/*g3'.format(args.caldatapath, source))
            files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                              int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]

            print('Analyzing source: {}'.format(source))

            if source not in data.keys():
                data[source] = {}
            for fname in files_to_parse:
                obsid = os.path.splitext(os.path.basename(fname))[0]
                cal_fname = '{}/{}/{}/nominal_online_cal.g3' \
                    .format(args.bolodatapath, source, obsid)                
                print('observation: {}'.format(obsid))

                if (obsid not in data[source].keys() or \
                    data[source][obsid]['timestamp'] != os.path.getctime(fname)) and \
                    os.path.exists(fname) and os.path.exists(cal_fname):
                    data[source][obsid] = {'timestamp': os.path.getctime(fname)}
                    try:
                        d = [fr for fr in core.G3File(fname)]

                        boloprops = [fr for fr in core.G3File(cal_fname)][0]["NominalBolometerProperties"]

                        for quantity_name in function_dict[source]:
                            func_result = function_dict[source][quantity_name](d[0], boloprops, selector_dict)
                            if func_result:
                                data[source][obsid][quantity_name] = func_result

                        if source in function_dict_raw:
                            rawpath = os.path.join(args.bolodatapath, source,
                                                   obsid, '0000.g3')
                            for quantity_name in function_dict_raw[source]:
                                func_result = function_dict_raw[source][quantity_name](rawpath, boloprops)
                                if func_result:
                                    data[source][obsid][quantity_name] = func_result
                    except Exception as ex:
                        print(ex)

        with open(datafile, 'wb') as f:
            pickle.dump(data, f)


# PLOT MODE
elif args.mode == 'plot':
    if args.timeinterval == 'monthly' or args.timeinterval == 'weekly':
        timeinterval_stub = args.timeinterval
    else: # checked arguments at top so we know this is castable to float
        timeinterval_stub = 'last_{}days'.format(int(args.timeinterval))

    plotsdir = os.path.join(args.outdir, 'plots')
    plotstimedir = os.path.join(args.outdir, 'plots', timeinterval_stub)
    datadir = os.path.join(args.outdir, 'data')

    # delete the full output directory tree if we are in rebuild mode
    if args.action == 'rebuild' and os.path.exists(plotstimedir):
        shutil.rmtree(plotstimedir)        

    weekly_filenames = np.array(glob(os.path.join(datadir, '*pkl')))
    weekly_datetimes = np.array([datetime.datetime.strptime(os.path.basename(fname).split('_')[0],
                                                            '%Y%m%d')
                                 for fname in weekly_filenames])
    ind_sort = np.argsort(weekly_datetimes)
    weekly_filenames = weekly_filenames[ind_sort]
    weekly_datetimes = weekly_datetimes[ind_sort]
    for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
        # convert min/max time for this interval to obsids that we can compare
        # to data obsids
        min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
        max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))

        if args.timeinterval == 'weekly':
            outdir = os.path.join(args.outdir, 'plots', timeinterval_stub,
                                  mindate.strftime('%Y%m%d'))
        elif args.timeinterval == 'monthly':
            outdir = os.path.join(args.outdir, 'plots', timeinterval_stub,
                                  mindate.strftime('%Y%m'))
        else:
            outdir = os.path.join(args.outdir, 'plots', timeinterval_stub)
        os.makedirs(outdir, exist_ok=True)

        # load data from this date range
        dt_ind = np.arange(len(weekly_datetimes))
        # case 1: plot time range is newer than all data timestamps; then load
        # last data file only
        if np.all(weekly_datetimes <= maxdate) and \
           np.all(weekly_datetimes <= mindate):
            dt_ind_inrange = np.array([len(weekly_datetimes) - 1])
        # case 2: plot time range is older than all data timestamps; then load
        # no data files
        elif np.all(weekly_datetimes >= maxdate) and \
             np.all(weekly_datetimes >= mindate):
            dt_ind_inrange = np.array([])
        # case 3: all others; then load all data files that are between the
        # endpoints of the plot time range, plus the data file with timestamp
        # that falls just before the beginning of the plot time range, if such
        # a file exists
        else:
            dt_ind_inrange = dt_ind[(weekly_datetimes <= maxdate) &
                                    (weekly_datetimes >= mindate)]
            if np.min(dt_ind_inrange) > 0:
                dt_ind_inrange = np.append(dt_ind_inrange,
                                           np.min(dt_ind_inrange) - 1)
        data = {}
        for ind in dt_ind_inrange:
            fname = weekly_filenames[ind]
            with open(fname, 'rb') as f:
                nextdata = pickle.load(f)
            for source in nextdata:
                if source in data:
                    data[source] = {**data[source], **nextdata[source]}
                else:
                    data[source] = nextdata[source]

        # restrict data to time range
        sourcelist = list(data.keys())
        for source in sourcelist:
            obsidlist = list(data[source].keys())
            for obsid in obsidlist:
                if int(obsid) <= min_obsid or int(obsid) >= max_obsid:
                    data[source].pop(obsid)

        # get maximum obsid that went into existing plots
        if os.path.exists(os.path.join(outdir, 'max_obsid.dat')):
            with open(os.path.join(outdir, 'max_obsid.dat'), 'r') as f:
                plot_obsid = int(f.readline())
        else:
            plot_obsid = 0

        # check whether data contains newer obsids than the plot, and if so,
        # remake the plots
        data_obsids = [int(obsid) for source in data.keys()
                       for obsid in data[source].keys()
                       if int(obsid) >= min_obsid and int(obsid) <= max_obsid]
        if any([plot_obsid < oid for oid in data_obsids]):
            # create the plots
            plot_median_cal_sn_4Hz(data, wafer_list, outdir, 'low')
            plot_median_cal_response_4Hz(data, wafer_list, outdir, 'low')
            plot_alive_bolos_cal_4Hz(data, wafer_list, outdir, 'low')
            plot_median_cal_sn_4Hz(data, wafer_list, outdir, 'high')
            plot_median_cal_response_4Hz(data, wafer_list, outdir, 'high')
            plot_alive_bolos_cal_4Hz(data, wafer_list, outdir, 'high')
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
            max_obsids = np.max([int(obsid) for source in data.keys()
                                 for obsid in data[source].keys()])
            with open(os.path.join(outdir, 'max_obsid.dat'), 'w') as f:
                f.write(str(max_obsid))

        # Create symlink from latest data directory to current
        # The latest data directory is not necessarily reliably identified by
        # the modification timestamp; instead, it is defined by having the
        # latest value of timestamp in its filename. Based on the naming scheme
        # for directories, an alphanumeric sort should also produce
        # chronological ordering, so we'll rely on this assumption.
        symlinkname = '{}/current'.format(plotstimedir)
        if os.path.exists(symlinkname):
            os.unlink(symlinkname)
        dirnames = np.sort(glob('{}/*'.format(plotstimedir)))
        latest_dirname = dirnames[-1]
        if args.timeinterval == 'monthly' or args.timeinterval == 'weekly':
            os.symlink(latest_dirname, '{}/temp'.format(plotstimedir))
        else:
            os.symlink(plotstimedir, '{}/temp'.format(plotstimedir))
        os.rename('{}/temp'.format(plotstimedir), '{}/current'.format(plotstimedir))
