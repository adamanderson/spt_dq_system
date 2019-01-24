from spt3g.std_processing import status_db
db = status_db.AutoprocDatabase(filename='/home/centos/spt_dq_system/autoproc.txt', user='centos', host='spt3g-dq.grid.uchicago.edu', read_only=True)
db.write_sql()
