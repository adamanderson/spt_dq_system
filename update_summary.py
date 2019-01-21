import os
import argparse as ap

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
args = parser.parse_args()

# update data skims
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
