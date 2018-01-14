from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
import numpy as np
import pickle as pickle
import argparse as ap
import glob
import datetime
import os.path

P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
P0.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
P0.add_argument('outfilename', action='store', default=None,
                help='name of output data file to write.')
P0.add_argument('--min-time', action='store', default='20180101',
                help='Minimum time of observations to skim. Format: YYYYMMDD')
P0.add_argument('--max-time', action='store',
                default=datetime.datetime.now().strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')
args = P0.parse_args()

# convert min/max times in observation IDs that we can compare with filenames
min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(args.min_time)))
max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(args.max_time)))

# functions that define quantities to be saved
def median_cal_sn(frame):
    return np.median([frame['CalibratorResponseSN'][bolo]
                      for bolo in frame['CalibratorResponseSN'].keys()])

function_dict = {'calibrator': {'MedianCalSN': median_cal_sn}}

data = {}

for source, quantities in function_dict.items():
    calfiles = glob.glob('{}/{}/*g3'.format(args.caldatapath, source))
    files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                          int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]
    
    print('Analyzing source: {}'.format(source))

    data[source] = {}
    for fname in files_to_parse:
        obsid = os.path.splitext(os.path.basename(fname))[0]
        print('observation: {}'.format(obsid))

        d = [fr for fr in core.G3File(fname)]
        for quantity_name in function_dict[source]:
            data[source][obsid] = function_dict[source][quantity_name](d[0])
            
with open(args.outfilename, 'wb') as f:
    pickle.dump(data, f)
