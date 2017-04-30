import sqlite3
import os
import glob

data_dir = '/spt/data/bolodata/downsampled/'
source_types = ['calibrator', 'CenA', 'elnod', 'RCW38', 'saturn']
dbname = 'obsfiles'
dbfile = '{}.db'.format(dbname)

# delete existing database file
if os.path.isfile(dbfile):
    os.remove(dbfile)

# create database and table
conn = sqlite3.connect(dbfile)
c = conn.cursor()
c.execute('CREATE TABLE {} (path text, type text, id int)'.format(dbname))

# fill the database
for source in source_types:
    files = glob.glob('{}/{}/*'.format(data_dir, source))
    for fname in files:
        c.execute("INSERT INTO {} VALUES ('{}', '{}', {})" \
                  .format(dbname, os.path.abspath(fname), source, int(os.path.basename(fname))))
conn.commit()
conn.close()

