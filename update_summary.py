import os
import argparse as ap
import datetime

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
        os.system('python summaryplot/cache_timeseries_data.py skim update '
                  '{} {} {} --min-time {}'.format(args.staticplotdir,
                                                  args.caldatadir,
                                                  args.bolodatadir,
                                                  args.mintime))

    # update weekly plots
    os.system('python summaryplot/cache_timeseries_data.py plot update weekly '
              '{} --min-time {}'.format(args.staticplotdir,
                                        args.mintime))

    # update monthly plots
    os.system('python summaryplot/cache_timeseries_data.py plot update monthly '
              '{} --min-time {}'.format(args.staticplotdir,
                                        args.mintime))

    # update "last 3 days" plots
    os.system('python summaryplot/cache_timeseries_data.py plot update 3 '
              '{} --min-time {}'.format(args.staticplotdir,
                                        args.mintime))

if args.mode == 'maps':
    # update coadded maps
    
    modes = ['coadding', 'plotting']
    if args.nocoadding:
        modes.remove('coadding')
    if args.noplotting:
        modes.remove('plotting')
    
    current_day = args.maxtime
    if current_day == None:
        current_day_flag = ""
    else:
        current_day_flag = "-c " + current_day
    
    if args.justseecommands:
        cmds_exec_flag = '-j'
    else:
        cmds_exec_flag = ''
    
    print()
    for n in ['7', '30']:
        for mode in modes:
            current_time = datetime.datetime.utcnow()
            print('Running commands for last_n and', mode, 'with n =', n)
            print('starting from', current_time, '(UTC) ...')
            print()
            if args.nologfiles:
                log_file_flag = ''
            else:
                log_file = args.coaddslogsdir + \
                           current_time.strftime('20%y%m%d_%H%M%S') + \
                           '_last_' + n + '_' + mode
                log_file_flag = '&>' + ' ' + log_file
            
            os.system('python summaryplot/cache_field_maps.py '
                      '-m {} -a update -t last_n -n {} '
                      '-o {} {} '
                      '-d {} -D {} '
                      '-F {} {} {}'.\
                      format(mode, n,
                             args.mintime, current_day_flag,
                             args.mapsdatadir, args.coaddsdatadir,
                             args.coaddsfigsdir, log_file_flag, cmds_exec_flag))
    
    for interval in ['yearly', 'monthly', 'weekly']:
        for mode in modes:
            current_time = datetime.datetime.utcnow()
            print("Running commands for", interval, "and", mode)
            print("starting from", current_time, "(UTC) ...")
            print()
            if args.nologfiles:
                log_file_flag = ''
            else:
                log_file = args.coaddslogsdir + \
                           current_time.strftime('20%y%m%d_%H%M%S') + \
                           "_" + interval + "_" + mode
                log_file_flag = '&>' + ' ' + log_file
            
            os.system('python summaryplot/cache_field_maps.py '
                      '-m {} -a update -t {} '
                      '-o {} {} '
                      '-d {} -D {} '
                      '-F {} {} {}'.\
                      format(mode, interval,
                             args.mintime, current_day_flag,
                             args.mapsdatadir, args.coaddsdatadir,
                             args.coaddsfigsdir, log_file_flag, cmds_exec_flag))

