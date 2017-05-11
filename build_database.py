import sqlite3
import os
import glob
import pandas as pd
import numpy as np
from spt3g import core
import argparse as ap

P0 = ap.ArgumentParser(description='Tool for building SPT3G data quality and '
                       'transfer databases.', formatter_class=ap.RawTextHelpFormatter)
P0.add_argument('mode', choices=['rebuild', 'update'], action='store',
                default=None, help='Action to perform on file database.')
args = P0.parse_args()


dbname = 'master_database'
dbfile = '{}.db'.format(dbname)

names = ['source', 'observation', 'status_fullrate', 'status_downsampled',
         'transfer_fullrate', 'transfer_downsampled', 'path', 'start_time']
typestr = {'source': 'text',
           'observation': 'int',
           'status_fullrate': 'text',
           'status_downsampled': 'text',
           'transfer_fullrate': 'int',
           'transfer_downsampled': 'int',
           'path': 'text',
           'start_time': 'float'}

# transfer database files
transfer_path = '/spt/data/transfer_database/'
transfer_file = '{}/transfer.txt'.format(transfer_path)

# fill the database from Sasha's transfer stuff to start
df_transfer = pd.read_csv(transfer_file, sep='\t',
                            dtype={'source': np.str,
                                    'observation': np.int,
                                    'status_fullrate': np.str,
                                    'status_downsampled': np.str,
                                    'transfer_fullrate': np.int,
                                    'transfer_downsampled': np.int})
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


conn = sqlite3.connect(dbfile)
c = conn.cursor()
for index, row in df_transfer.iterrows():
    result = c.execute('SELECT * FROM master_database WHERE observation == {}'.format(row['observation'])).fetchall()
    if not result:
        # get the path of the observation
        path = '/spt/data/bolodata/downsampled/{}/{}/' \
            .format(row['source'], row['observation'])
        if not os.path.isdir(path):
            path = 'none'

        # get start time by skimming files
        first_file = path + '0000.g3'
        try:
            f = core.G3File(first_file)
            print('Processing {}'.format(first_file))
            frame = f.next()
            time = frame['ObservationStart'].time
        except:
            time = np.nan

        newrow = [row[key] for key in row.keys()] + [path, time]
        print("INSERT INTO {} VALUES ({})" \
                      .format(dbname, ', '.join(['\'{}\''.format(var) for var in newrow])))
        c.execute("INSERT INTO {} VALUES ({})" \
                      .format(dbname, ', '.join(['\'{}\''.format(var) for var in newrow])))

conn.commit()
conn.close()
