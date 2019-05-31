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
                         help='Path to coadded maps.')
parser_maps.add_argument('coaddsfigsdir',
                         help='Path to figures related to coadded maps.')
parser_maps.add_argument('coaddslogsdir',
                         help='Path to the log files.')
parser_maps.add_argument('mintime',
                         help='Time from which to start plots.')



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
    current_time = datetime.datetime.utcnow()
    logfile = args.coaddslogsdir + current_time.strftime('20%y%m%d_%H%M%S')

    for mode in ['coadding', 'plotting']:
        for interval in ['last_n', 'weekly', 'monthly', 'yearly']:
            os.system('python summaryplot/cache_field_maps.py '
                      '-m {} -a update -t {} '
                      '-o {} -d {} -D {} -F {} &> {}'.format(mode,
                                                            interval,
                                                            args.mintime,
                                                            args.mapsdatadir,
                                                            args.coaddsdatadir,
                                                            args.coaddsfigsdir,
                                                            logfile))


    """
    os.system('python summaryplot/cache_field_maps.py -m coadding -a update -t last_n -n 3'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m coadding -a update -t weekly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m coadding -a update -t monthly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m coadding -a update -t yearly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    # update figures related to maps



    os.system('python summaryplot/cache_field_maps.py -m plotting -a update -t last_n -n 3'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m plotting -a update -t weekly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m plotting -a update -t monthly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))

    os.system('python summaryplot/cache_field_maps.py -m plotting -a update -t yearly'
              '-o {} -d {} -D {} -F {}'.format(args.mintime,
                                               args.mapsdatadir
                                               args.coaddsdatadir
                                               args.coaddsfigsdir))
    """


