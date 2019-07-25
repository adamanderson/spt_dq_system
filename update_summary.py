import os
import sys
import argparse as ap
import logging
import time
import datetime
from summaryplot import cache_timeseries_data
from summaryplot import cache_field_maps

parser = ap.ArgumentParser(description='Wrapper for updating data and plots.',
                           formatter_class=ap.ArgumentDefaultsHelpFormatter)

S = parser.add_subparsers(dest='mode', title='subcommands',
                          help='Set of plots to update. For help, call: '
                          '%(prog)s %(metavar)s -h')
parser_summarystats = S.add_parser('summarystats', help='Updates skimmed data '
                                   'and plots that appear on the `Summary '
                                   'Stats` tab of the webpage')
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

parser_maps = S.add_parser('maps', help='Updates the coadd data and map '
                           'monitoring plots.')
parser_maps.add_argument('mapsdatadir',
                         help='Path to map data.')
parser_maps.add_argument('coaddsdatadir',
                         help='Path to coadded maps to be created.')
parser_maps.add_argument('coaddsfigsdir',
                         help='Path to figures related to coadded maps.')
parser_maps.add_argument('coaddslogsdir',
                         help='Path to log files to be generated.')
parser_maps.add_argument('mintime',
                         help='Time from which to care about.')
parser_maps.add_argument('-m', '--maxtime', action='store',
                         type=str, default=None,
                         help='Time up to which to care about.')
parser_maps.add_argument('-c', '--nocoadding',
                         action='store_true', default=False,
                         help='Skip the coadding part.')
parser_maps.add_argument('-p', '--noplotting',
                         action='store_true', default=False,
                         help='Skip the plotting part.')
parser_maps.add_argument('-r', '--rebuildinsteadofupdate',
                         action='store_true', default=False,
                         help='Rebuild maps and figures instead of update them.')
parser_maps.add_argument('-l', '--nologfiles',
                         action='store_true', default=False,
                         help='Whether to redirect stdout/err to txt files.')
parser_maps.add_argument('-j', '--justseecommands',
                         action='store_true', default=False,
                         help='Whether to just see commands to be run.')

args = parser.parse_args()


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
        for n in [7, 30]:
            current_time = datetime.datetime.utcnow()
            logger.info('Running commands for last_n and %s with n = %s', mode, n)
            logger.info('starting from %s (UTC) ...', current_time)
            logger.info('')
            
            logger_name = 'last_{}_{}'.format(n, mode)
            time_string = current_time.strftime('20%y%m%d_%H%M%S')
            if args.nologfiles:
                log_file = None
            else:
                log_file = os.path.join(args.coaddslogsdir,
                                        '{}_{}'.format(time_string, logger_name))
            
            try:
                cache_field_maps.update(mode, action, time_interval='last_n',
                                        last_how_many_days=n,
                                        oldest_time_to_consider=args.mintime,
                                        current_time=current_day,
                                        original_maps_dir=args.mapsdatadir,
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
        
        
        for interval in ['yearly', 'weekly', 'monthly']:
            current_time = datetime.datetime.utcnow()
            logger.info('Running commands for %s and %s', interval, mode)
            logger.info('starting from %s (UTC) ...', current_time)
            logger.info('')
            
            logger_name = '{}_{}'.format(interval, mode)
            time_string = current_time.strftime('20%y%m%d_%H%M%S')
            if args.nologfiles:
                log_file = None
            else:
                log_file = os.path.join(args.coaddslogsdir,
                                        '{}_{}'.format(time_string, logger_name))
            
            try:
                cache_field_maps.update(mode, action, time_interval=interval,
                                        oldest_time_to_consider=args.mintime,
                                        current_time=current_day,
                                        original_maps_dir=args.mapsdatadir,
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

