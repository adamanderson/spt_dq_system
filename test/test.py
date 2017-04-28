import sqlite3
import os

files = os.listdir(os.getcwd())
print files

conn = sqlite3.connect('filelist.db')
c = conn.cursor()
c.execute('create table filetable(id integer primary key, one varchar(100), two info);')

for fname in files:
    c.execute('insert into filetable values(null, \'{}\', 10);'.format(fname))
# c.execute()
conn.commit()
c.execute('select * from filetable;')
print c.fetchall()
