# ============================================================================ #
#  This script is intended to be called by db_server.js, which periodically    #
#  calls this script to update the data products and figures used by the       #
#  summary webpage.                                                            #
#                                                                              #
# ============================================================================ #

import os
import sys
import time
import datetime
import argparse
import logging

from matplotlib import use
use('agg')

from summaryplot import cache_timeseries_data
from summaryplot import cache_field_maps
from summaryplot import cache_archive_data



# ==============================================================================
# Parse command line arguments
# ------------------------------------------------------------------------------


parser = argparse.ArgumentParser(
             description='Wrapper for updating data products and figures used '
                         'by the summary webpage.',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

S = parser.add_subparsers(
        dest='mode', title='subcommands',
        help='Whether to update the data products and figures used for the '
             '"Calibration Observations", "Field Observations", or '
             '"Weather etc." section of the data quality summary webpage.'
             'For help, call: %(prog)s %(metavar)s -h')


parser_summarystats = \
    S.add_parser('summarystats',
                 help='Updates skimmed data and plots that appear on the '
                      '"Calibration Observation" tab of the webpage.')

parser_summarystats.add_argument('staticplotdir',
                                 help='Path to static plots and data skims.')
parser_summarystats.add_argument('caldatadir',
                                 help='Path to calibration data.')
parser_summarystats.add_argument('bolodatadir',
                                 help='Path to raw bolometer data.')
parser_summarystats.add_argument('mintime',
                                 help='Time from which to start plots.')
parser_summarystats.add_argument('--no-data-update', action='store_true',
                                 help='Do not update the skimmed data.')


parser_maps = \
    S.add_parser('maps',
                 help='Updates the coadded maps and the figures showing '
                      'map-related quantities that appear in the '
                      '"Field Observations" tabs of the webpage.')

parser_maps.add_argument('season',
                         help='Season from which the field maps are '
                              'to be monitored. Available choices are '
                              '"winter", "summer", "summerb", and "summerc". '
                              'If "winter" is chosen, '
                              'then the relevant fields are '
                              'the four ra0hdec* fields. '
                              'If "summer*" is chosen, '
                              'then the relevant fields are '
                              'the fourteen ra5hdec* fields, '
                              'four ra1h40dec* fields, '
                              'or four ra12h40dec* fields.')
parser_maps.add_argument('mapsdatadir',
                         help='Path to the directory containing auto-processed '
                              'maps of individual observations.')
parser_maps.add_argument('coaddsdatadir',
                         help='Path to the directory where coadded maps from '
                              'different time intervals will be saved.')
parser_maps.add_argument('coaddsfigsdir',
                         help='Path to the directory where figures showing '
                              'coadded maps and map-related quantities such as '
                              'noise levels will be saved.')
parser_maps.add_argument('coaddslogsdir',
                         help='Path to the directory where the logs generated '
                              'during the coadding and plotting processes '
                              'will be saved.')
parser_maps.add_argument('caldatadir',
                         help='Path to the directory where results of the '
                              'auto-processing of the calibration observation '
                              'are saved.')
parser_maps.add_argument('bolodatadir',
                         help='Path to the directory where raw timestreams '
                              'are saved.')
parser_maps.add_argument('mintime',
                         help='Any maps from observations taken before this '
                              'time will be ignored.')
parser_maps.add_argument('-m', '--maxtime',
                         action='store', default=None, type=str,
                         help='Any maps from observations taken after time '
                              'time will be ignored.')
parser_maps.add_argument('-c', '--nocoadding',
                         action='store_true', default=False,
                         help='Skip the coadding part.')
parser_maps.add_argument('-p', '--noplotting',
                         action='store_true', default=False,
                         help='Skip the plotting part.')
parser_maps.add_argument('-o', '--onlymonthly',
                         action='store_true', default=False,
                         help='Only invoke the monthly mode of updating.')
parser_maps.add_argument('-e', '--onlyyearly',
                         action='store_true', default=False,
                         help='Only invoke the yearly mode of updating.')
parser_maps.add_argument('-r', '--rebuildinsteadofupdate',
                         action='store_true', default=False,
                         help='Remake coadded maps and relevant figures '
                              'instead of updating the existing files.')
parser_maps.add_argument('-l', '--nologfiles',
                         action='store_true', default=False,
                         help='Log files are not to be generated.')
parser_maps.add_argument('-j', '--justseecommands',
                         action='store_true', default=False,
                         help='Just see the functions that will be called and '
                              'the arguments that will be supplied to them '
                              'without actually invoking the functions.')


parser_arc = \
    S.add_parser('arcs',
                 help='Updates the pickle files and figures relevant to the '
                      '"Weather etc." tab of the summary webpage, which shows '
                      'weather-related data recorded in the ARC files.')
parser_arc.add_argument('arcdir',
                        help='Path to the directory that contains '
                             'all the ARC files.')
parser_arc.add_argument('arcpklsdir',
                        help='Path to the directory where pickle files '
                             'storing reduced data from the ARC files '
                             'will be saved.')
parser_arc.add_argument('arcfigsdir',
                        help='Path to the directory where figures showing '
                             'the reduced data will be saved.')
parser_arc.add_argument('arclogsdir',
                        help='Path to the directory where logs generated '
                             'during the data reduction and plotting processes '
                             'will be saved.')
parser_arc.add_argument('mintime',
                        help='ARC file data taken from before this day '
                             'will not be used for reduction and plotting. '
                             'The format of the argument is YYYYMMDD.')
parser_arc.add_argument('-m', '--maxtime',
                        action='store', default=None, type=str,
                        help='ARC file data taken from after this day '
                             'will not be used for reduction and plotting. '
                             'The format of the argument is YYYYMMDD.')
parser_arc.add_argument('-r', '--noreduction',
                        action='store_true', default=False,
                        help='Skip the data reduction part.')
parser_arc.add_argument('-f', '--nofigures',
                        action='store_true', default=False,
                        help='Skip the figure generation part.')
parser_arc.add_argument('-b', '--rebuild',
                        action='store_true', default=False,
                        help='Remake all the pickle files and/or figures '
                             'instead of updating them by just incorporating '
                             'newer data.')
parser_arc.add_argument('-d', '--dryrun',
                        action='store_true', default=False,
                        help='Just see what data reduction and plotting '
                             'functions will be called with what arguments '
                             'without actually invoking those functions.')


args = parser.parse_args()


# ==============================================================================



# ==============================================================================
# Update the data products and figures used by
# the "Calibration Observations" tab
# ------------------------------------------------------------------------------

if args.mode == 'summarystats':
    
    # update data skims
    if args.no_data_update == False:
        cache_timeseries_data.update(mode='skim', action='update',
                                     outdir=args.staticplotdir,
                                     caldatapath=args.caldatadir,
                                     bolodatapath=args.bolodatadir,
                                     min_time=args.mintime)

    # update plots
    for interval in ['weekly', 'monthly', 'yearly']:
        cache_timeseries_data.update(mode='plot', action='update',
                                     timeinterval=interval,
                                     outdir=args.staticplotdir,
                                     bolodatapath=args.bolodatadir,
                                     min_time=args.mintime)

    for n in [1, 3, 7, 30]:
        cache_timeseries_data.update(mode='plot', action='update',
                                     timeinterval=n,
                                     outdir=args.staticplotdir,
                                     bolodatapath=args.bolodatadir,
                                     min_time=args.mintime)



# Update the data products and figures used by
# the "Field Observations (Winter)" and "Field Observations (Summer)" tabs
# ------------------------------------------------------------------------------

if args.mode == 'maps':
    
    # update coadded maps and relevant figures
    
    for d in [args.coaddsdatadir, args.coaddsfigsdir, args.coaddslogsdir]:
        if not os.path.isdir(d):
            os.mkdir(d)
    
    current_time = datetime.datetime.utcnow()
    
    logger_name = "overall_progress"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter('%(message)s')
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(log_format)
    logger.addHandler(stream_handler)
    
    if args.nologfiles:
        pass
    else:
        time_string = current_time.strftime('20%y%m%d_%H%M%S')
        log_file = os.path.join(args.coaddslogsdir,
                                '{}_{}'.format(time_string, logger_name))
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(log_format)
        logger.addHandler(file_handler)
    
    logger.info('\n')
    logger.info('Starting to update coadded maps and figures at')
    logger.info('%s (UTC) ...', current_time)
    logger.info('')
    time.sleep(1)
    
    modes = ['coadding', 'plotting']
    if args.nocoadding:
        modes.remove('coadding')
    if args.noplotting:
        modes.remove('plotting')
        
    if args.rebuildinsteadofupdate:
        action = 'rebuild'
    else:
        action = 'update'
    
    current_day = args.maxtime
    
    for mode in modes:
        if args.onlymonthly or args.onlyyearly:
            ns = []
        else:
            ns = [7, 30]
        for n in ns:
            current_time = datetime.datetime.utcnow()
            logger.info('Running commands for last_n and %s with n = %s',
                        mode, n)
            logger.info('starting from %s (UTC) ...', current_time)
            logger.info('')
            
            logger_name = 'last_{}_{}'.format(n, mode)
            time_string = current_time.strftime('20%y%m%d_%H%M%S')
            if args.nologfiles:
                log_file = None
            else:
                log_file = \
                    os.path.join(args.coaddslogsdir,
                                 '{}_{}'.format(time_string, logger_name))
            
            try:
                cache_field_maps.update(args.season,
                                        mode, action, time_interval='last_n',
                                        last_how_many_days=n,
                                        oldest_time_to_consider=args.mintime,
                                        current_time=current_day,
                                        original_maps_dir=args.mapsdatadir,
                                        calibration_data_dir=args.caldatadir,
                                        bolo_timestreams_dir=args.bolodatadir,
                                        coadds_dir=args.coaddsdatadir,
                                        figs_dir=args.coaddsfigsdir,
                                        just_see_commands=args.justseecommands,
                                        logger_name=logger_name,
                                        log_file=log_file)
            except Exception:
                logger.exception('Something did not go well '
                                 'while updating maps/figures ...')
                sys.exit()
            time.sleep(1)
        
        
        if args.onlymonthly:
            intervals = ['monthly']
        elif args.onlyyearly:
            intervals = ['yearly']
        else:
            intervals = ['weekly', 'yearly', 'monthly']
        for interval in intervals:
            current_time = datetime.datetime.utcnow()
            logger.info('Running commands for %s and %s', interval, mode)
            logger.info('starting from %s (UTC) ...', current_time)
            logger.info('')
            
            logger_name = '{}_{}'.format(interval, mode)
            time_string = current_time.strftime('20%y%m%d_%H%M%S')
            if args.nologfiles:
                log_file = None
            else:
                log_file = \
                    os.path.join(args.coaddslogsdir,
                                 '{}_{}'.format(time_string, logger_name))
            
            try:
                cache_field_maps.update(args.season,
                                        mode, action, time_interval=interval,
                                        oldest_time_to_consider=args.mintime,
                                        current_time=current_day,
                                        original_maps_dir=args.mapsdatadir,
                                        calibration_data_dir=args.caldatadir,
                                        bolo_timestreams_dir=args.bolodatadir,
                                        coadds_dir=args.coaddsdatadir,
                                        figs_dir=args.coaddsfigsdir,
                                        just_see_commands=args.justseecommands,
                                        logger_name=logger_name,
                                        log_file=log_file)
            except Exception:
                logger.exception('Something did not go well '
                                 'while updating maps/figures ...')
                sys.exit()
            time.sleep(1)
    
    current_time = datetime.datetime.utcnow()
    logger.info('All done at %s (UTC)!', current_time)
    logger.info('')



# Update the data products and figures used by the "Weather etc." tab
# ------------------------------------------------------------------------------

if args.mode == 'arcs':

    # Update pickle files containing reduced arc file data and relevant figures
    
    for d in [args.arcpklsdir, args.arcfigsdir, args.arclogsdir]:
        if not os.path.isdir(d):
            os.mkdir(d)
    
    current_time = datetime.datetime.utcnow()
    
    logger_name = 'overall_progress'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter('%(message)s')
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(log_format)
    logger.addHandler(stream_handler)
    
    time_string = current_time.strftime('20%y%m%d_%H%M%S')
    log_file = os.path.join(args.arclogsdir,
                            '{}_{}'.format(time_string, logger_name))
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)
    
    logger.info('\n')
    logger.info('Starting to update arc files-related data products at')
    logger.info('%s (UTC) ...', current_time)
    logger.info('')
    time.sleep(1)
    
    modes = ['reducing', 'plotting']
    if args.noreduction:
        modes.remove('reducing')
    if args.nofigures:
        modes.remove('plotting')
    
    if args.rebuild:
        action = 'rebuild'
    else:
        action = 'update'
    
    for mode in modes:        
        if mode == 'reducing':
            current_time = datetime.datetime.utcnow()
            logger.info('Reducing the ARC data and pickling the results...')
            logger.info('(Started at %s (UTC))', current_time)
            
            logger_name = 'reducing'
            time_string = current_time.strftime('20%y%m%d_%H%M%S')
            log_file    = os.path.join(args.arclogsdir,
                                       '{}_{}'.format(time_string, logger_name))
            
            try:
                cache_archive_data.run(mode=mode,
                                       action=action,
                                       oldest_time_to_consider=args.mintime,
                                       newest_time_to_consider=args.maxtime,
                                       arc_files_dir=args.arcdir,
                                       pickles_dir=args.arcpklsdir,
                                       just_see_commands=args.dryrun,
                                       logger_name=logger_name,
                                       log_file=log_file)
            except:
                logger.exception('\nSomething did not go well...')
                sys.exit()
            logger.info('')
            time.sleep(1)
    
        elif mode == 'plotting':
            for interval in ['last_01', 'last_03', 'last_07', 'last_30',
                             'weekly',  'monthly', 'fridge']:
                current_time = datetime.datetime.utcnow()
                logger.info('Plotting the pickled results '
                            'for the time interval "%s"...', interval)
                logger.info('(Started at %s (UTC))', current_time)
                
                logger_name = 'plotting_{:s}'.format(interval)
                time_string = current_time.strftime('20%y%m%d_%H%M%S')
                log_file    = os.path.join(args.arclogsdir,
                                           '{}_{}'.format(
                                               time_string, logger_name))
                
                try:
                    calfigsdir = args.arcfigsdir.replace('arcs/figs', 'plots')
                    cache_archive_data.run(mode=mode,
                                           action=action,
                                           interval=interval,
                                           oldest_time_to_consider=args.mintime,
                                           newest_time_to_consider=args.maxtime,
                                           pickles_dir=args.arcpklsdir,
                                           figures_dir=args.arcfigsdir,
                                           calibra_dir=calfigsdir,
                                           just_see_commands=args.dryrun,
                                           logger_name=logger_name,
                                           log_file=log_file)
                except:
                    logger.exception('\nSomething did not go well...')
                    sys.exit()
                logger.info('')
                time.sleep(1)
                
    current_time = datetime.datetime.utcnow()
    logger.info('All done at %s (UTC)!', current_time)
    logger.info('')
                



# ==============================================================================
