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
P0.add_argument('--dbroot', default='/scratch/arahlin/public_html/spt_dq_system',
                help='Path to directory where the database is stored')
args = P0.parse_args()


dbname = 'aux_database'
dbfile = os.path.join(args.dbroot, '{}.db'.format(dbname))

if not os.path.exists(dbfile) and args.mode is None:
    args.mode = 'rebuild'

names = ['filename', 'type', 'size', 'status', 'modified', 'path']
modnames = ['type', 'size', 'status', 'modified']
typestr = {'filename': 'text',
           'type': 'text',
           'size': 'int',
           'status': 'text',
           'modified': 'text',
           'path': 'text'}

# transfer database files
transfer_path = '/spt/data/transfer_database/'
transfer_file = '{}/aux_transfer.txt'.format(transfer_path)

# fill the database from Sasha's transfer stuff to start
df_transfer = pd.read_csv(
    transfer_file, sep='\t',
    dtype={'filename': np.str,
           'type': np.str,
           'size': np.int,
           'status': np.str,
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


conn = sqlite3.connect(dbfile)
c = conn.cursor()
for index, row in df_transfer.iterrows():
    tp = row['type']
    fname = row['filename']

    result = c.execute('SELECT * FROM {} WHERE filename == \'{}\''.format(
            dbname, row['filename'])).fetchall()

    path = 'none'
    if result:
        rdict = dict(zip(names, result[0]))
        path = rdict['path']

    if path == 'none' and row['status'] == 'verified':
        if tp == 'pydfmux':
            root = '/scratch/pydfmux_output'
        elif tp == 'pydfmuxlog':
            root = '/spt/data/pydfmux_output'
            fname = 'aux-pydfmuxlog_{}.tar.gz'.format(row['filename'])
        elif tp in ['arc', 'eht', 'tar', 'rsync']:
            root = os.path.join('/spt/data', tp)
            if tp == 'rsync':
                fname = os.path.basename(row['filename'])
            elif tp == 'tar':
                root = '/scratch' + os.path.dirname(row['filename'])
                fname = os.path.basename(row['filename'])

        path = os.path.join(root, fname)

    if not result:

        if not os.path.exists(path):
            path = 'none'

        print('Processing aux entry {}'.format(row['filename']))
        newrow = [row[key] for key in row.keys()] + [path]
        c.execute("INSERT INTO {} VALUES ({})".format(
                dbname, ', '.join(['\'{}\''.format(var) for var in newrow])))

    else:

        for name in names:
            if name == 'path':
                if rdict['path'] == path:
                    rdict.pop('path')
                else:
                    if path != 'none' and os.path.exists(path):
                        rdict['path'] = path
            elif name not in row:
                rdict.pop(name)
            elif row[name] == rdict[name]:
                rdict.pop(name)
            elif rdict[name] == 'nan':
                try:
                    isnan = np.isnan(row[name])
                except TypeError:
                    pass
                else:
                    if isnan:
                        rdict.pop(name)
                    else:
                        rdict[name] = row[name]
            else:
                rdict[name] = row[name]

        if rdict:
            print('Updating aux entry {}'.format(row['filename']))
            arglist = ', '.join(['{} = \'{}\''.format(k, v) for k, v in rdict.items()])
            c.execute('UPDATE {} SET {} WHERE filename == \'{}\''.format(
                    dbname, arglist, row['filename']))

conn.commit()
conn.close()
