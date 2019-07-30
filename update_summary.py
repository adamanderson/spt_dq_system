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

from summaryplot import cache_timeseries_data
from summaryplot import cache_field_maps



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
             '"summarystats" section or maps "section". For help, call: '
             '%(prog)s %(metavar)s -h')


parser_summarystats = \
    S.add_parser('summarystats',
                 help='Updates skimmed data and plots that appear on the '
                      '"Summary Stats" tab of the webpage.')

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
                      'map-related quantities that appear on the "Maps" tab
                      'of the webpage.')

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

args = parser.parse_args()


# ==============================================================================



# ==============================================================================
# Update the data products and figures used by the "Summary Stats" section
# ------------------------------------------------------------------------------

if args.mode == 'summarystats':
    
    # update data skims
    if args.no_data_update == False:
        cache_timeseries_data.update(mode='skim', action='update',
                                     outdir=args.staticplotdir,
                                     caldatapath=args.caldatadir,
                                     bolodatapath=args.bolodatadir,
                                     min_time=args.mintime)

    # update weekly plots
    cache_timeseries_data.update(mode='plot', action='update',
                                 timeinterval='weekly',
                                 outdir=args.staticplotdir,
                                 min_time=args.mintime)

    # update monthly plots
    cache_timeseries_data.update(mode='plot', action='update',
                                 timeinterval='monthly',
                                 outdir=args.staticplotdir,
                                 min_time=args.mintime)

    # update "last 3 days" plots
    cache_timeseries_data.update(mode='plot', action='update',
                                 timeinterval='3',
                                 outdir=args.staticplotdir,
                                 min_time=args.mintime)



# Update the data products and figures used by the "Maps" section
# ------------------------------------------------------------------------------

if args.mode == 'maps':
    
    # update coadded maps
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
                cache_field_maps.update(mode, action, time_interval='last_n',
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
                cache_field_maps.update(mode, action, time_interval=interval,
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


# ==============================================================================
