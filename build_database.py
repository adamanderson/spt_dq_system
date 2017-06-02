import sqlite3
import os
import glob
import pandas as pd
import numpy as np
from spt3g import core, cluster
import argparse as ap
import datetime
import tarfile

P0 = ap.ArgumentParser(description='Tool for building SPT3G data quality and '
                       'transfer databases.', formatter_class=ap.RawTextHelpFormatter)
P0.add_argument('mode', choices=['rebuild', 'update'], action='store',
                default=None, help='Action to perform on file database.')
P0.add_argument('plotdir', default=None, help='Directory where plots should be '
                'stored and updated.')
P0.add_argument('logdir', default=None, help='Directory where log files should '
                'be stored and updated.')
P0.add_argument('proxycert', default=None, help='Proxy certificate for plot '
                'generation on the grid.')
P0.add_argument('--dbroot', default='/scratch/arahlin/public_html/spt_dq_system',
                help='Path to directory where the database is stored')
# P0.add_argument('--overwrite-plots', default=True,
#                 help='Regenerate plots that already exist for any entry in the '
#                 'database that script updates. Note that if run in \'rebuild\' '
#                 'mode, this will regenerate all plots for the entire dataset.')
args = P0.parse_args()


dbname = 'master_database'
dbfile = os.path.join(args.dbroot, '{}.db'.format(dbname))

if not os.path.exists(dbfile) and args.mode is None:
    args.mode = 'rebuild'

names = ['source', 'observation', 'status_fullrate', 'status_downsampled',
         'transfer_fullrate', 'transfer_downsampled', 'modified',
         'plots_complete', 'plots_inprogress',
         'path', 'autoproc_path', 'start_time', 'date', 'time']
modnames = ['status_fullrate', 'status_downsampled', 'transfer_fullrate',
            'transfer_downsampled', 'modified']
typestr = {'source': 'text',
           'observation': 'int',
           'status_fullrate': 'text',
           'status_downsampled': 'text',
           'transfer_fullrate': 'int',
           'transfer_downsampled': 'int',
           'modified': 'text',
           'plots_complete': 'int',
           'plots_inprogress': 'int',
           'path': 'text',
           'autoproc_path': 'text',
           'start_time': 'float',
           'date': 'text',
           'time': 'text'}

# transfer database files
transfer_path = '/spt/data/transfer_database/'
transfer_file = '{}/transfer.txt'.format(transfer_path)

# fill the database from Sasha's transfer stuff to start
df_transfer = pd.read_csv(
    transfer_file, sep='\t',
    dtype={'source': np.str,
           'observation': np.int,
           'status_fullrate': np.str,
           'status_downsampled': np.str,
           'transfer_fullrate': np.int,
           'transfer_downsampled': np.int,
           'modified': np.str})

if args.mode == 'rebuild':
    # delete existing database file
    if os.path.isfile(dbfile):
        os.remove(dbfile)
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()
    column_schema = ', '.join(['{} {}'.format(col, typestr[col]) for col in names])
    c.execute('CREATE TABLE {} ({})'.format(dbname, column_schema))
    conn.commit()
    conn.close()

n_files_per_worker = 100

conn = sqlite3.connect(dbfile)
conn.row_factory = sqlite3.Row
c = conn.cursor()
for index, row in df_transfer.iterrows():
    result = c.execute('SELECT * FROM {} WHERE observation == {}'.format(
            dbname, row['observation'])).fetchall()

    # add entry to the sqlite DB if it doesn't exist already
    if not result:
        start_time = np.nan
        date = 'nan'
        time = 'nan'

        # get the path of the observation
        path = '/spt/data/bolodata/downsampled/{}/{}/' \
            .format(row['source'], row['observation'])
        if not os.path.isdir(path):
            print('ERROR processing {}: Missing path'.format(path))
            path = 'none'
        autoproc_path = '/spt/user/production/calibration/{}/{}.g3' \
            .format(row['source'], row['observation'])
        if not os.path.exists(autoproc_path):
            autoproc_path = 'none'
        else:
            # get start time by skimming files
            first_file = path + '0000.g3'
            try:
                # f = core.G3File(first_file)
                print('Processing {}'.format(first_file))
                # frame = f.next()
                # start_time = frame['ObservationStart'].time
                # dt = datetime.datetime.fromtimestamp(
                #     int(start_time/core.G3Units.second))
                # date = dt.strftime('%d/%m/%Y')
                # time = dt.strftime('%H:%M')
                date = 0
                time = 0
            except Exception as e:
                print('ERROR processing {}: {}'.format(first_file, e))
                path = 'none'

        newrow = {key: row[key] for key in row.keys()}
        newrow['path'] = path
        newrow['autoproc_path'] = autoproc_path
        newrow['start_time'] = start_time
        newrow['date'] = date
        newrow['time'] = time
        newrow['plots_complete'] = int(os.path.exists('{}/{}/'.format(args.plotdir, newrow['observation'])))
        newrow['plots_inprogress'] = 0
        # replace the line below with a call to plotter.check_plots when implemented
        keylist = ', '.join(['\'{}\''.format(key) for key in newrow.keys()])
        datalist = ', '.join(['\'{}\''.format(newrow[key]) for key in newrow.keys()])
        c.execute("INSERT INTO {} ({}) VALUES ({})".format(
                dbname, keylist, datalist))

    # if the entry already exists, update any modified rows
    else:
        rdict = dict(zip(names, result[0]))
        for name in names:
            if name not in row:
                rdict.pop(name)
            elif row[name] == rdict[name]:
                rdict.pop(name)
            elif np.isnan(row[name]) and rdict[name] == 'nan':
                rdict.pop(name)

        if rdict:
            print('Updating observation {}/{}'.format(
                    row['source'], row['observation']))
            arglist = ', '.join(['{} = {}'.format(k, v) for k, v in rdict.items()])
            c.execute('UPDATE {} SET {} WHERE observation == {}'.format(
                    dbname, arglist, row['observation']))

    # get the newly added row for further processing
    result = c.execute('SELECT * FROM {} WHERE observation == {}'.format(
                       dbname, row['observation'])).fetchone()

# handle plot generation
# collect files that need to be processed
result = c.execute('SELECT * FROM {} WHERE plots_complete == 0'.format(dbname)).fetchall()
plot_queue = [] # list of observations to process
for row in result:
    if row['plots_inprogress'] == 0 and row['autoproc_path'] != 'none':
        plot_queue.append((os.path.dirname(row['autoproc_path']).split('/')[-1],
                           row['autoproc_path'],
                           '{}/nominal_online_cal.g3'.format(row['path'])))
        c.execute('UPDATE {} SET {} WHERE observation == {}'.format(
                    dbname, ('plots_inprogress = 1'), row['observation']))

conn.commit()
conn.close()

jbatch = 0
timestamp = '{:%Y%m%d_%H%M%S}'.format(datetime.datetime.now())
while jbatch*n_files_per_worker < len(plot_queue):
    # write the pathlist file
    startfile = jbatch*n_files_per_worker
    stopfile = np.min([(jbatch+1)*n_files_per_worker, 
                       len(plot_queue)])
    pathlist_file = '{}_{}.txt'.format(timestamp, jbatch)
    output_file = '{}_{}_plots.tar.gz'.format(timestamp, jbatch)

    # tar up relevant g3 files
    tar_input_fname = '{}_{}_input_data.tar'.format(timestamp, jbatch)
    tar_cal_fname = '{}_{}_online_cal.tar'.format(timestamp, jbatch)
    tar_input = tarfile.open(tar_input_fname, 'w')
    tar_cal = tarfile.open(tar_cal_fname, 'w')
    for p in plot_queue[startfile:stopfile]:
        tar_input.add(p[1])
        tar_cal.add(p[2])
    tar_input.close()
    tar_cal.close()

    input_files = []
    with open(pathlist_file, 'w') as f:
        for p in plot_queue[startfile:stopfile]:
            f.write('\t'.join(p) + '\n')
            input_files.append(p[1])
            input_files.append(p[2])
            
    condor_configs = {
        'script': 'plotter.py',
        'args': ['{}'.format(pathlist_file), 'plots'],
        'log_root': args.logdir,
        'request_cpus': 1,
        'request_memory': 2*core.G3Units.GB,
        'request_disk': 10*core.G3Units.GB,
        'requirements': '(OpSysAndVer == "CentOS6" || OpSysAndVer == "RedHat6" || OpSysAndVer == "SL6")',
        'grid_proxy': '{}'.format(args.proxycert),
        'when_to_transfer_output': 'ON_EXIT',
        'aux_input_files': [pathlist_file],
        'input_files': [tar_input_fname, tar_cal_fname],
        'output_root': args.plotdir,
        'output_files': [output_file],
        'globus_uri': 'gsiftp://gridftp.grid.uchicago.edu:2811/cephfs',
        'user_code': '\n'.join(['tar -xvf {}'.format(tar_input_fname),
                                'tar -xvf {}'.format(tar_cal_fname)])
    }
    print tar_input_fname
    cluster.condor_submit(**condor_configs)
    jbatch = jbatch+1
