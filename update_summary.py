import os
import argparse as ap
import datetime

parser = ap.ArgumentParser(description='Wrapper for updating data and plots.',
                           formatter_class=ap.ArgumentDefaultsHelpFormatter)
parser.add_argument('staticplotdir',
                    help='Path to static plots and data skims.')
parser.add_argument('caldatadir',
                    help='Path to calibration data.')
parser.add_argument('bolodatadir',
                    help='Path to raw bolometer data.')
parser.add_argument('mintime',
                    help='Time from which to start plots.')
parser.add_argument('--no-data-update', action='store_true',
                    help='Do not update the skimmed data.')

"""parser.add_argument('mapsdatadir',
                    help='Path to map data.')
parser.add_argument('coaddsdatadir',
                    help='Path to coadded maps.')
parser.add_argument('coaddsfigsdir',
                    help='Path to figures related to coadded maps.')
parser.add_argument('coaddslogsdir',
                    help='Path to the log files.')"""



args = parser.parse_args()



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


# update coadded maps

"""current_day = "20190519"

print()
for n in ['7']:
    for mode in ['coadding', 'plotting']:
        current_time = datetime.datetime.utcnow()
        logfile = args.coaddslogsdir + current_time.strftime('20%y%m%d_%H%M%S') + "_last_" + n + "_" + mode
        print("Running commands for last_n and", mode, "with n =", n)
        print("starting from", current_time, "...")
        print()
        os.system('python summaryplot/cache_field_maps.py '
                  '-m {} -a update -t last_n -n {} '
                  '-o {} -c {} -d {} -D {} -F {} &> {}'.\
                  format(mode, n,
                         args.mintime, current_day,
                         args.mapsdatadir, args.coaddsdatadir, args.coaddsfigsdir, logfile))
        

for interval in ['weekly', 'monthly', 'yearly']:
    for mode in ['coadding', 'plotting']:
        current_time = datetime.datetime.utcnow()
        logfile = args.coaddslogsdir + current_time.strftime('20%y%m%d_%H%M%S') + "_" + interval + "_" + mode
        print("Running commands for", interval, "and", mode)
        print("starting from", current_time, "...")
        print()
        os.system('python summaryplot/cache_field_maps.py '
                  '-m {} -a update -t {} '
                  '-o {} -c {} -d {} -D {} -F {} &> {}'.\
                  format(mode, interval,
                         args.mintime, current_day,
                         args.mapsdatadir, args.coaddsdatadir, args.coaddsfigsdir, logfile))"""

