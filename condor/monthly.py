from glob import glob
import tarfile
import os
import datetime
import time
import numpy as np
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils,cluster
import argparse as ap

'''
Generate monthly summary plots for a given source, this submits condor jobs
'''

P = ap.ArgumentParser(description='Taring files for condor',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)


P.add_argument('source', action='store',
           default=None, help='name of source')
P.add_argument('year',action='store', default=None,
           help='the date range of data to be plotted')

args = P.parse_args()
source=args.source
year=args.year

#directories to be modified
plotpath= '/scratch/dyang18/dq_plots/'
script='/home/dyang18/spt3g/spt3g_software/scratch/dyang18/'+source+'PlotCondor.py'
log_root='/scratch/dyang18'
datatarpath='/spt/user/dyang18/tar/'+str(year)+'/'+str(year)+'_'+source+'.tar.gz'
bolotarpath='/spt/user/dyang18/tar/'+str(year)+'/'+str(year)+'_'+source+'BolometerProperties.tar.gz'
grid_proxy='/home/dyang18/dyang'
outputpath='/spt/user/dyang18/output/'


currentDT = datetime.datetime.now()
ctime= currentDT.strftime("%Y-%m")
a=[os.path.basename(f) for f in glob(plotpath+'*')]
t=[str(x)+'-'+y for x in [2018] for y in ['01','02','03','04','05','06','07','08','09','10','11','12']]


#final=[date for date in t if date not in a]

final=[]
for ti in t:
    a=[os.path.basename(f) for f in glob(plotpath+ ti+'/*')]
    ck= ''.join(a)
    if source not in ck:
        final.append(ti)



final.append(ctime)
final=list(np.unique(final))
#make updated plots for the current month
print('making plots for '+ ' '.join(final))

#Submission section for different months of data



input_files=[datatarpath,bolotarpath]
output_files=[source+".tar.gz"]

args=[ '-d' +' ' +' '.join(final)]
jobname='condor_{}'.format(time.strftime('%Y%m%d_%H%M%S'))
output_root=outputpath+jobname
cluster.condor_submit(script=script, args=args, caller='python', jobname=jobname,
                      log_root=log_root, output_root=output_root, input_files=input_files,
                      output_files=output_files, grid_proxy=grid_proxy, aux_input_files=[],
                      aux_output_files=[], spt3g_env=True, python3=False,
                      user_code='', requirements=None, request_cpus=1,
                      request_memory=4*core.G3Units.GB,
                      request_disk=6*core.G3Units.GB, queue=1,
                      when_to_transfer_output='ON_EXIT', create_only=False,
                      globus_uri='gsiftp://gridftp.grid.uchicago.edu:2811/cephfs')
