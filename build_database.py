import sqlite3
import os
import glob
import pandas as pd
import numpy as np

dbname = 'master_database'

# transfer database files
transfer_path = '/spt/data/transfer_database/'
transfer_file = '{}/transfer.txt'.format(transfer_path)

dbfile = '{}.db'.format(dbname)

# delete existing database file
if os.path.isfile(dbfile):
    os.remove(dbfile)

typestr = {np.dtype('O'): 'text',
            np.dtype('int64'): 'int'}

# fill the database
df_transfer = pd.read_csv(transfer_file, sep='\t',
                            dtype={'source': np.str,
                                    'observation': np.int,
                                    'status_fullrate': np.str,
                                    'status_downsampled': np.str,
                                    'transfer_fullrate': np.int,
                                    'transfer_downsampled': np.int})
df_transfer.fillna('NULL', inplace=True)
def generate_path(row):
    path = '/spt/data/bolodata/downsampled/{}/{}/' \
            .format(row['source'], row['observation'])
    if os.path.isdir(path):
        return path
    else:
        return 'none'
df_transfer['path'] = df_transfer.apply(lambda row: generate_path(row), axis=1)

# create database and table
conn = sqlite3.connect(dbfile)
c = conn.cursor()
column_schema = ', '.join(['{} {}'.format(var, typestr[df_transfer[var].dtype]) for var in df_transfer.columns])
c.execute('CREATE TABLE {} ({})'.format(dbname, column_schema))

for index, row in df_transfer.iterrows():
    c.execute("INSERT INTO {} VALUES ({})" \
                  .format(dbname, ', '.join(['\'{}\''.format(row[var]) if type(row[var]) is str 
                                             else '{}'.format(row[var]) for var in df_transfer.columns])))

conn.commit()
conn.close()
