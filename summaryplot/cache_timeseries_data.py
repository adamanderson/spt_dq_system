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
import re
from summaryplot.calibrator import *
from summaryplot.elnod import *
from summaryplot.noise import *
from summaryplot.htwo import *
from summaryplot.focus import *
import hashlib

from spt3g.autoprocessing.status_db import AutoprocDatabase, ScanifyDatabase, field_regex

db = None
sdb = None


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

    if mode == 'plot':
        if timeinterval == 'weekly':
            next_day = next_day - datetime.timedelta(days=next_day.weekday())
        elif timeinterval == 'monthly':
            next_day = datetime.datetime(year=next_day.year, month=next_day.month, day=1)
        elif timeinterval == 'yearly':
            next_day = datetime.datetime(year=next_day.year, month=1, day=1)

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
            elif timeinterval == 'yearly': # yearly mode
                next_day = datetime.datetime(year=next_day.year+1,
                                             month=1,
                                             day=1)
        date_boundaries.append(next_day)
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

        function_dict = {'RCW38':             {'RCW38SkyTransmission': htwo_sky_transmission},
                         'RCW38-pixelraster': {'MedianRCW38FluxCalibration': median_htwo_fluxcal,
                                               'MedianRCW38IntegralFlux': median_htwo_intflux,
                                               'AliveBolosRCW38': alive_bolos_htwo_fluxcal},
                         'MAT5A':             {'MAT5ASkyTransmission': htwo_sky_transmission},
                         'MAT5A-pixelraster': {'MedianMAT5AFluxCalibration': median_htwo_fluxcal,
                                               'MedianMAT5AIntegralFlux': median_htwo_intflux,
                                               'AliveBolosMAT5A': alive_bolos_htwo_fluxcal},
                         'W28A2':             {'W28A2SkyTransmission': htwo_sky_transmission},
                         'W28A2-pixelraster': {'MedianW28A2FluxCalibration': median_htwo_fluxcal,
                                               'MedianW28A2IntegralFlux': median_htwo_intflux,
                                               'AliveBolosW28A2': alive_bolos_htwo_fluxcal},
                         'IRAS17258':             {'IRAS17258SkyTransmission': htwo_sky_transmission},
                         'IRAS17258-pixelraster': {'MedianIRAS17258FluxCalibration': median_htwo_fluxcal,
                                                   'MedianIRAS17258IntegralFlux': median_htwo_intflux,
                                                   'AliveBolosIRAS17258': alive_bolos_htwo_fluxcal},
                         'RCW122A':             {'RCW122ASkyTransmission': htwo_sky_transmission},
                         'RCW122A-pixelraster': {'MedianRCW122AFluxCalibration': median_htwo_fluxcal,
                                                 'MedianRCW122AIntegralFlux': median_htwo_intflux,
                                                 'AliveBolosRCW122A': alive_bolos_htwo_fluxcal},
                         'PMNJ0210-5101':     {'BenchPosAndFittingResults': benchpos_min_fwhm_ellip},
                         'PMNJ0522-3628':     {'BenchPosAndFittingResults': benchpos_min_fwhm_ellip},
                         'ra0hdec-52.25':     {'NominalBenchPosition': extract_benchpos},
                         'ra1h40dec-33.25':   {'NominalBenchPosition': extract_benchpos},
                         'ra5hdec-33.25':     {'NominalBenchPosition': extract_benchpos},
                         'ra12h30dec-33.25':  {'NominalBenchPosition': extract_benchpos},
                         'calibrator':        {'MedianCalSN_4Hz': median_cal_sn_4Hz,
                                               'MedianCalResponse_4Hz': median_cal_response_4Hz,
                                               'AliveBolosCal_4Hz': alive_bolos_cal_4Hz,
                                               'elevation': mean_cal_elevation},
                         'elnod':             {'MedianElnodIQPhaseAngle': median_elnod_iq_phase_angle,
                                               'MedianElnodSNSlopes': median_elnod_sn_slope,
                                               'AliveBolosElnod': alive_bolos_elnod,
                                               'MedianPelecPerAirmass': median_pelec_per_airmass,
                                               'MedianElnodOpacity': median_elnod_opacity},
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
                                               'NET_10.0Hz_to_15.0Hz': median_net_10Hz_to_15Hz,
                                               'NumberOfLinesInMedianPSDs': number_of_lines_in_median_psds}}
        function_dict_raw = {}
        key_dict_raw = {}


        global db
        if db is None:
            db = AutoprocDatabase(read_only=True)

        global sdb
        if sdb is None:
            sdb = ScanifyDatabase(read_only=True)

        # loop over weeks
        for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
            smindate = mindate.strftime('%Y%m%d')
            smaxdate = maxdate.strftime('%Y%m%d')
            print('Updating cache files between {} and {}'.format(smindate, smaxdate))
            # convert min/max times to observation IDs that we can compare with filenames
            min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(smindate)))
            max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(smaxdate)))

            datafile = os.path.join(datadir, '{}_data_cache.pkl'.format(smindate))
            updated = False
            if os.path.exists(datafile):
                try:
                    with open(datafile, 'rb') as f:
                        data = pickle.load(f)
                except Exception as e:
                    raise e.__class__('Error loading {}: {}'.format(datafile, e))
            else:
                data = {}
                updated = True

            data['updated_sources'] = set()

            # update the data skim
            for source, quantities in function_dict.items():
                # incremental update after every source
                if updated:
                    with open(datafile, 'wb') as f:
                        pickle.dump(data, f)
                    updated = False

                print('Analyzing source: {}'.format(source))
                if source not in data.keys():
                    data[source] = {}
                    updated = True
                    data['updated_sources'].add(source)

                isfieldobs = re.match(field_regex, source) is not None
                if isfieldobs and bolodatapath.startswith('/sptgrid'):
                    bolodatapath = bolodatapath.replace('fullrate', 'downsampled')
                    data_rate = 'downsampled'
                else:
                    bolodatapath = bolodatapath.replace('downsampled', 'fullrate')
                    data_rate = 'fullrate'

                entries = db.match(
                    '{}/calframe'.format(source) if isfieldobs else source,
                    '{}:{}'.format(min_obsid, max_obsid + 1),
                    status='complete',
                    return_df=True,
                )

                for _, entry in entries.iterrows():
                    obsid = int(entry['observation'])

                    if source in focus_sources:
                        if not sdb.check_ready(source, obsid, rate=data_rate):
                            # need data on disk to grab observation frame
                            continue

                    if isfieldobs:
                        fname = os.path.join(bolodatapath, source, str(obsid))
                    else:
                        fname = os.path.join(caldatapath, source, "{}.g3".format(obsid))
                    cal_fname = os.path.join(bolodatapath, source, str(obsid), 'nominal_online_cal.g3')
                    rawpath = os.path.join(bolodatapath, source, str(obsid), '0000.g3')
                    print('observation: {}'.format(obsid))

                    # make sure all obsidss are ints
                    for k in list(data[source]):
                        if isinstance(k, str):
                            data[source][int(k)] = data[source].pop(k)

                    if obsid in data[source]:
                        # replace ctime with database modified time
                        tstamp = data[source][obsid]['timestamp']
                        if isinstance(tstamp, float):
                            ctime = os.path.getctime(fname)
                            if ctime == tstamp:
                                data[source][obsid]['timestamp'] = entry['modified']
                                updated = True
                                data['updated_sources'].add(source)
                            else:
                                del data[source][obsid]

                    if (obsid not in data[source] or \
                        data[source][obsid]['timestamp'] != entry['modified']):
                        updated = True
                        data[source][obsid] = {'timestamp': entry['modified']}
                        data['updated_sources'].add(source)

                        """
                        if caldatapath.startswith("/sptgrid"):
                            def get_gfal_copied_file(path):
                                origin = "gsiftp://osg-gridftp.grid.uchicago.edu:2811" + path
                                destin = os.path.join(os.getcwd(), datadir, os.path.basename(path))
                                os.system("gfal-copy {} file://{} -t 86400".format(origin, destin))
                                return destin
                            if not isfieldobs:
                                fname = get_gfal_copied_file(fname)
                            cal_fname = get_gfal_copied_file(cal_fname)
                            if isfieldobs or (source in function_dict_raw):
                                rawpath = get_gfal_copied_file(rawpath)
                        """
                        ### We can just read the files directly if gfal-copy is cumbersome or too slow.
                        ### Remember to also take care of the two lines below that delete files.

                        if source in focus_sources:
                            iterator = core.G3File(rawpath)
                            while True:
                                frame = iterator.next()
                                if frame.type == core.G3FrameType.Observation:
                                    obf = frame
                                    break
                            d = [obf]
                            if not isfieldobs:
                                d += [fr for fr in core.G3File(fname)]
                        else:
                            d = [fr for fr in core.G3File(fname)][0]
                            if 'ObservationID' not in d:
                                d['ObservationID'] = obsid
                        boloprops = [fr for fr in core.G3File(cal_fname)][0]["NominalBolometerProperties"]
                        if len(d) == 0:
                            data[source].pop(obsid)
                            continue

                        for quantity_name in function_dict[source]:
                            func_result = function_dict[source][quantity_name](d, boloprops, selector_dict)
                            if func_result:
                                data[source][obsid][quantity_name] = func_result

                        if source in function_dict_raw:
                            for quantity_name in function_dict_raw[source]:
                                func_result = function_dict_raw[source][quantity_name](
                                    rawpath, cal_fname, fname, boloprops, selector_dict)
                                if func_result:
                                    for actual_name, result in zip(key_dict_raw[quantity_name], func_result):
                                        data[source][obsid][actual_name] = result

                        """
                        if caldatapath.startswith("/sptgrid"):
                            for f in [fname, cal_fname, rawpath]:
                                if not f.startswith("/sptgrid"):
                                    os.system("rm {}".format(f))
                        """

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
            smindate = mindate.strftime('%Y%m%d')
            smaxdate = maxdate.strftime('%Y%m%d')
            print('Updating {} plots between {} and {}'.format(timeinterval, smindate, smaxdate))
            # convert min/max time for this interval to obsids that we can compare
            # to data obsids
            min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(smindate)))
            max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(smaxdate)))

            if timeinterval == 'weekly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub, smindate)
            elif timeinterval == 'monthly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub, smindate[:6])
            elif timeinterval == 'yearly':
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub, smindate[:4])
            else:
                plotsdir = os.path.join(outdir, 'plots', timeinterval_stub,
                                        'last_{:02d}'.format(int(timeinterval)))
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
                    old_pkl_hash = f.readline().strip()
            else:
                # if hash file doesn't exist, set to empty string to force rebuild
                old_pkl_hash = ''

            # compute hashes for each file individually
            new_pkl_hash_dict = {os.path.basename(f): file_hash([f]) for f in filelist}
            hash_dict_filename = os.path.join(plotsdir, 'data_hash_dict.dat')
            old_pkl_hash_dict = {os.path.basename(k): '' for k in filelist}
            if os.path.exists(hash_dict_filename):
                with open(hash_dict_filename, 'r') as f:
                    for line in f:
                        k, v = line.strip().split()
                        old_pkl_hash_dict[os.path.basename(k)] = v

            # rebuild the data if the new and old hashes differ or if a rebuild
            # is requested
            if new_pkl_hash != old_pkl_hash or action == 'rebuild':
                # load all the data
                data = {}
                updated_sources = set()
                for fname in filelist:
                    fname1 = os.path.basename(fname)
                    with open(fname, 'rb') as f:
                        nextdata = pickle.load(f)
                        # make sure all obsidss are ints
                        for source, dd in nextdata.items():
                            if source == 'updated_sources':
                                continue
                            for k in list(dd):
                                if isinstance(k, str):
                                    dd[int(k)] = dd.pop(k)

                    next_updated_sources = nextdata.pop('updated_sources', set())
                    if action == 'rebuild' or bolodatapath.startswith('/sptgrid'):
                        # rebuild all plots in chicago because the rsync'ed cache files
                        # are not updated at the same rate as at pole
                        updated_sources |= set(nextdata.keys())
                    elif new_pkl_hash_dict[fname1] != old_pkl_hash_dict[fname1]:
                        updated_sources |= next_updated_sources
                    for source in nextdata:
                        if source in data:
                            data[source].update(nextdata[source])
                        else:
                            data[source] = nextdata[source]

                # restrict data to time range
                sourcelist = list(data.keys())
                for source in sourcelist:
                    obsidlist = list(data[source].keys())
                    for obsid in obsidlist:
                        # convert obsids to ints
                        if type(obsid) is str:
                            data[source][int(obsid)] = data[source][obsid]
                            data[source].pop(obsid)

                    obsidlist = list(data[source].keys())
                    for obsid in obsidlist:
                        if obsid <= min_obsid or obsid >= max_obsid:
                            data[source].pop(obsid)

                opts = dict(xlims=[mindate, maxdate])

                print('Found updated sources: {}'.format(list(updated_sources)))

                # create the plots
                if timeinterval == 'yearly':
                    for src in set(htwo_sources) & updated_sources:
                        plot_median_htwo_fluxcal(src, data, wafer_list, plotsdir, **opts)
                        plot_median_htwo_intflux(src, data, wafer_list, plotsdir, **opts)
                        plot_alive_bolos_htwo(src, data, wafer_list, plotsdir, **opts)

                    if len(set(focus_sources) & updated_sources) > 0:
                        plot_focus_quasar_fitting_results(data, plotsdir, **opts)
                else:
                    for src in set(htwo_sources) & updated_sources:
                        plot_median_htwo_fluxcal(src, data, wafer_list, plotsdir, **opts)
                        plot_median_htwo_intflux(src, data, wafer_list, plotsdir, **opts)
                        plot_alive_bolos_htwo(src, data, wafer_list, plotsdir, **opts)
                        plot_htwo_skytrans(src, data, wafer_list, plotsdir, **opts)

                    if len(set(focus_sources) & updated_sources) > 0:
                        plot_focus_quasar_fitting_results(data, plotsdir, **opts)

                    if 'calibrator' in updated_sources:
                        plot_median_cal_sn_4Hz(data, wafer_list, plotsdir, 'low', **opts)
                        plot_median_cal_response_4Hz(data, wafer_list, plotsdir, 'low', **opts)
                        plot_alive_bolos_cal_4Hz(data, wafer_list, plotsdir, 'low', **opts)
                        plot_median_cal_sn_4Hz(data, wafer_list, plotsdir, 'high', **opts)
                        plot_median_cal_response_4Hz(data, wafer_list, plotsdir, 'high', **opts)
                        plot_alive_bolos_cal_4Hz(data, wafer_list, plotsdir, 'high', **opts)

                    if 'elnod' in updated_sources:
                        plot_median_elnod_response(data, wafer_list, plotsdir, **opts)
                        plot_median_elnod_sn(data, wafer_list, plotsdir, **opts)
                        plot_median_elnod_opacity(data, wafer_list, plotsdir, **opts)
                        plot_alive_bolos_elnod(data, wafer_list, plotsdir, **opts)
                        plot_median_elnod_iq_phase(data, wafer_list, plotsdir, **opts)

                    if 'noise' in updated_sources:
                        for quant in ['NEI', 'NEP', 'NET']:
                            for rng in [
                                '0.1Hz_to_0.5Hz',
                                '1.0Hz_to_2.0Hz',
                                '3.0Hz_to_5.0Hz',
                                '10.0Hz_to_15.0Hz',
                            ]:
                                key = '{}_{}'.format(quant, rng)
                                plot_median_noise(data, key, wafer_list, plotsdir, **opts)
                        plot_number_of_lines(data, wafer_list, plotsdir, **opts)

                # update the hash
                with open(hash_filename, 'w') as f:
                    f.write('{}\n'.format(new_pkl_hash))
                with open(hash_dict_filename, 'w') as f:
                    for k, v in sorted(new_pkl_hash_dict.items()):
                        f.write('{} {}\n'.format(os.path.basename(k), v))

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
    parser_plot.add_argument('bolodatapath', default=None,
                             help='Path to bolometer data.')
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
               bolodatapath=args.bolodatapath,
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
