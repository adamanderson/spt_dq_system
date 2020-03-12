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
from summaryplot.calibrator import *
from summaryplot.elnod import *
from summaryplot.mat5a import *
from summaryplot.noise import *
from summaryplot.rcw38 import *
from summaryplot.pmnj0210_5101 import *
import hashlib


def file_hash(filelist):
    '''
    Utility function for computing md5 checksum hash of a list of files. Useful
    for keeping track of whether data files have changed when generating plots.

    Parameters
    ----------
    filelist : list
        List of files for which to compute md5 sum hash

    Returns
    -------
    hash : str
        Hex hash
    '''
    md5 = hashlib.md5()
    for fname in filelist:
        with open(fname, 'rb') as f:
            while True:
                data = f.read(65536) # read data in 64kB chunks
                if not data:
                    break
                md5.update(data)
    return md5.hexdigest()


def update(mode, action, outdir, caldatapath=None, bolodatapath=None,
           timeinterval=None, min_time=None, max_time=None):
    timenow = datetime.datetime.utcnow()
    dt = datetime.timedelta(-1*(timenow.weekday()+1))
    default_mintime = timenow + dt

    # check timeinterval argument
    if mode == 'plot' and timeinterval not in ['weekly', 'monthly', 'yearly']:
        try:
            int(timeinterval)
        except:
            raise ValueError('Argument `timeinterval` is none of `monthly`, '
                             '`weekly` or a number of days.')

    # check --max-time argument
    if max_time == None:
        if timeinterval == 'yearly':
            max_time = timenow.replace(year=timenow.year+1, month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
        if timeinterval == 'monthly':
            if timenow.month != 12:
                max_time = timenow.replace(month=timenow.month+1, day=1, hour=0, minute=0, second=0, microsecond=0)
            else:
                max_time = timenow.replace(year=timenow.year+1, month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
        elif timeinterval == 'weekly':
            max_time = (timenow + datetime.timedelta(days=7-timenow.weekday())).replace(hour=0, minute=0, second=0, microsecond=0)
        else:
            max_time = (timenow + datetime.timedelta(days=1))
        max_time = max_time.strftime('%Y%m%d')

    wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180',
                  'w181', 'w187', 'w188', 'w201', 'w203',
                  'w204', 'w206', 'all']

    # parse times to loop over
    dt_mintime = datetime.datetime(year=int(min_time[:4]),
                                   month=int(min_time[4:6]),
                                   day=int(min_time[6:8]))
    dt_maxtime = datetime.datetime(year=int(max_time[:4]),
                                   month=int(max_time[4:6]),
                                   day=int(max_time[6:8]))

    date_boundaries = []
    next_day = dt_mintime

    if mode == 'skim' or (timeinterval in ['weekly', 'monthly', 'yearly']):
        while next_day < dt_maxtime:
            date_boundaries.append(next_day)

            if mode == 'skim' or timeinterval == 'weekly': # weekly mode
                next_day = next_day + datetime.timedelta(days = 7 - next_day.weekday())
            elif timeinterval == 'monthly': # monthly mode
                try:
                    next_day = datetime.datetime(year=next_day.year,
                                                 month=next_day.month+1,
                                                 day=1)
                except ValueError:
                    next_day = datetime.datetime(year=next_day.year+1,
                                                 month=1,
                                                 day=1)
            elif timeinterval == 'yearly':
                next_day = datetime.datetime(year=next_day.year+1,
                                             month=1,
                                             day=1)
        date_boundaries.append(dt_maxtime)
    else:
        date_boundaries = [dt_maxtime - datetime.timedelta(days=int(timeinterval)), dt_maxtime]


    # SKIM MODE
    if mode == 'skim':
        # manage directory structure
        os.makedirs(os.path.join(outdir, 'data'),
                         exist_ok=True)

        # delete the full output directory tree if we are in rebuild mode
        datadir = os.path.join(outdir, 'data')
        if action == 'rebuild' and os.path.exists(datadir):
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
                         'PMNJ0210-5101':     {'BenchPosAndFittingResults': benchpos_min_fwhm_ellip},
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
        function_dict_raw = {'calibrator':    {'el': mean_cal_elevation},
                             'elnod'     :    {'ppa': median_pelec_per_airmass}}
        key_dict_raw = {'el' : ['elevation'],
                        'ppa': ['MedianPelecPerAirmass', 'MedianElnodOpacity']}


        # loop over weeks
        for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
            # convert min/max times to observation IDs that we can compare with filenames
            min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
            max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))

            datafile = os.path.join(datadir, '{}_data_cache.pkl'.format(mindate.strftime('%Y%m%d')))
            updated = False
            if os.path.exists(datafile):
                with open(datafile, 'rb') as f:
                    data = pickle.load(f)
            else:
                data = {}
                updated = True

            # update the data skim
            for source, quantities in function_dict.items():
                calfiles = glob('{}/calibration/{}/*g3'.format(caldatapath, source))
                files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                                  int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]

                print('Analyzing source: {}'.format(source))

                if source not in data.keys():
                    data[source] = {}
                for fname in files_to_parse:
                    obsid = os.path.splitext(os.path.basename(fname))[0]
                    cal_fname = '{}/{}/{}/nominal_online_cal.g3' \
                        .format(bolodatapath, source, obsid)                
                    print('observation: {}'.format(obsid))

                    if (obsid not in data[source].keys() or \
                        data[source][obsid]['timestamp'] != os.path.getctime(fname)) and \
                        os.path.exists(fname) and os.path.exists(cal_fname):
                        updated = True

                        data[source][obsid] = {'timestamp': os.path.getctime(fname)}
                        d = [fr for fr in core.G3File(fname)]
                        boloprops = [fr for fr in core.G3File(cal_fname)][0]["NominalBolometerProperties"]
                        if source not in ['PMNJ0210-5101']:
                            d = d[0]

                        for quantity_name in function_dict[source]:
                            func_result = function_dict[source][quantity_name](d, boloprops, selector_dict)
                            if func_result:
                                data[source][obsid][quantity_name] = func_result

                        if source in function_dict_raw:
                            rawpath = os.path.join(bolodatapath, source,
                                                   obsid, '0000.g3')
                            for quantity_name in function_dict_raw[source]:
                                func_result = function_dict_raw[source][quantity_name](rawpath, cal_fname, fname, boloprops, selector_dict)
                                if func_result:
                                    for actual_name, result in zip(key_dict_raw[quantity_name], func_result):
                                        data[source][obsid][actual_name] = result

            if updated:
                with open(datafile, 'wb') as f:
                    pickle.dump(data, f)


    # PLOT MODE
    elif mode == 'plot':
        if timeinterval in ['weekly', 'monthly', 'yearly']:
            timeinterval_stub = timeinterval
        else: # checked arguments at top so we know this is castable to float
            timeinterval_stub = 'last_n'

        plotstimedir = os.path.join(outdir, 'plots', timeinterval_stub)
        datadir = os.path.join(outdir, 'data')

        # delete the full output directory tree if we are in rebuild mode
        if action == 'rebuild' and os.path.exists(plotstimedir):
            shutil.rmtree(plotstimedir)        

        weekly_filenames = np.array(glob(os.path.join(datadir, '*pkl')))
        weekly_datetimes = np.array([datetime.datetime.strptime(os.path.basename(fname).split('_')[0], '%Y%m%d') for fname in weekly_filenames])
        ind_sort = np.argsort(weekly_datetimes)
        weekly_filenames = weekly_filenames[ind_sort]
        weekly_datetimes = weekly_datetimes[ind_sort]
        for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
            # convert min/max time for this interval to obsids that we can compare
            # to data obsids
            min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
            max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))

            if timeinterval == 'weekly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub,
                                        mindate.strftime('%Y%m%d'))
            elif timeinterval == 'monthly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub,
                                        mindate.strftime('%Y%m'))
            elif timeinterval == 'yearly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub,
                                        mindate.strftime('%Y'))
            else:
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub,
                                        'last_{}'.format(timeinterval))
            os.makedirs(plotsdir, exist_ok=True)

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

            # compute hash of data files to check whether plots need updating
            filelist = np.sort([weekly_filenames[ind] for ind in dt_ind_inrange])
            new_pkl_hash = file_hash(filelist)
            hash_filename = os.path.join(plotsdir, 'data_hash.dat')
            if os.path.exists(hash_filename):
                with open(hash_filename, 'r') as f:
                    old_pkl_hash = f.readline()
            else:
                # if hash file doesn't exist, set to empty string to force rebuild
                old_pkl_hash = ''

            # rebuild the data if the new and old hashes differ or if a rebuild
            # is requested
            if new_pkl_hash != old_pkl_hash or action == 'rebuild':
                # load all the data
                data = {}
                for fname in filelist:
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

                # create the plots
                if timeinterval == 'yearly':
                    plot_median_rcw38_fluxcal(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_median_rcw38_intflux(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_alive_bolos_rcw38(data, wafer_list, plotsdir,
                                           xlims=[mindate, maxdate])
                    plot_median_mat5a_fluxcal(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_median_mat5a_intflux(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_pmnj0210_5101_fitting_results(data, plotsdir,
                                                       xlims=[mindate, maxdate])
                    plot_alive_bolos_mat5a(data, wafer_list, plotsdir,
                                           xlims=[mindate, maxdate])
                else:
                    plot_median_cal_sn_4Hz(data, wafer_list, plotsdir, 'low',
                                           xlims=[mindate, maxdate])
                    plot_median_cal_response_4Hz(data, wafer_list, plotsdir, 'low',
                                                 xlims=[mindate, maxdate])
                    plot_alive_bolos_cal_4Hz(data, wafer_list, plotsdir, 'low',
                                             xlims=[mindate, maxdate])
                    plot_median_cal_sn_4Hz(data, wafer_list, plotsdir, 'high',
                                           xlims=[mindate, maxdate])
                    plot_median_cal_response_4Hz(data, wafer_list, plotsdir, 'high',
                                                 xlims=[mindate, maxdate])
                    plot_alive_bolos_cal_4Hz(data, wafer_list, plotsdir, 'high',
                                             xlims=[mindate, maxdate])
                    plot_median_elnod_response(data, wafer_list, plotsdir,
                                               xlims=[mindate, maxdate])
                    plot_median_elnod_sn(data, wafer_list, plotsdir,
                                         xlims=[mindate, maxdate])
                    plot_median_elnod_opacity(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_alive_bolos_elnod(data, wafer_list, plotsdir,
                                           xlims=[mindate, maxdate])
                    plot_median_elnod_iq_phase(data, wafer_list, plotsdir,
                                               xlims=[mindate, maxdate])
                    plot_median_rcw38_fluxcal(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_median_rcw38_intflux(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_alive_bolos_rcw38(data, wafer_list, plotsdir,
                                           xlims=[mindate, maxdate])
                    plot_rcw38_sky_transmission(data, wafer_list, plotsdir,
                                                xlims=[mindate, maxdate])
                    plot_median_mat5a_fluxcal(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_median_mat5a_intflux(data, wafer_list, plotsdir,
                                              xlims=[mindate, maxdate])
                    plot_alive_bolos_mat5a(data, wafer_list, plotsdir,
                                           xlims=[mindate, maxdate])
                    plot_mat5a_sky_transmission(data, wafer_list, plotsdir,
                                                xlims=[mindate, maxdate])
                    plot_pmnj0210_5101_fitting_results(data, plotsdir,
                                                       xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEI_0.1Hz_to_0.5Hz', wafer_list, plotsdir, 
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEI_1.0Hz_to_2.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEI_3.0Hz_to_5.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEI_10.0Hz_to_15.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEP_0.1Hz_to_0.5Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEP_1.0Hz_to_2.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEP_3.0Hz_to_5.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NEP_10.0Hz_to_15.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NET_0.1Hz_to_0.5Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NET_1.0Hz_to_2.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NET_3.0Hz_to_5.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])
                    plot_median_noise(data, 'NET_10.0Hz_to_15.0Hz', wafer_list, plotsdir,
                                      xlims=[mindate, maxdate])

                # update the hash
                with open(os.path.join(plotsdir, 'data_hash.dat'), 'w') as f:
                    f.write(new_pkl_hash)


if __name__ == '__main__':

    timenow = datetime.datetime.utcnow()
    dt = datetime.timedelta(-1*(timenow.weekday()+1))
    default_mintime = timenow + dt
    default_maxtime = timenow

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
                             'YYYYMMDD (starts at beginning of day)')
    parser_skim.add_argument('--max-time', action='store',
                             default=default_maxtime.strftime('%Y%m%d'),
                             help='Maximum time of observations to skim. Format: '
                             'YYYYMMDD (ends at end of day)')

    # plot data mode
    parser_plot = S.add_parser('plot',
                               help='Plot data from skims of autoprocessing data.')
    parser_plot.add_argument('action', choices=['update', 'rebuild'], default=None,
                             help='Update or rebuild the data skims or plots.')
    parser_plot.add_argument('timeinterval', default=None,
                             help='Time interval at which to generate plots. '
                             'Available options are '
                             '"weekly", "monthly", "yearly" or "N", '
                             'where N is generates a plot containing data from '
                             'only the most recent N days.')
    parser_plot.add_argument('outdir', action='store', default=None,
                             help='Path containing skimmed data and plots.')
    parser_plot.add_argument('--min-time', action='store',
                             default=default_mintime.strftime('%Y%m%d'),
                             help='Minimum time of observations to skim. Format: '
                             'YYYYMMDD (starts at beginning of day)')
    parser_plot.add_argument('--max-time', action='store',
                             default=default_maxtime.strftime('%Y%m%d'),
                             help='Maximum time of observations to skim. Format: '
                             'YYYYMMDD (ends at end of day)')
    args = P0.parse_args()
    
    
    if args.mode == 'plot':
        update(mode=args.mode,
               action=args.action,
               timeinterval=args.timeinterval,
               outdir=args.outdir,
               min_time=args.min_time,
               max_time=args.max_time)
    elif args.mode == 'skim':
        update(mode=args.mode,
               action=args.action,
               outdir=args.outdir,
               caldatapath=args.caldatapath,
               bolodatapath=args.bolodatapath,
               min_time=args.min_time,
               max_time=args.max_time)
