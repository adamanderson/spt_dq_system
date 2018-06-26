from glob import glob
import tarfile
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils,cluster
import argparse as ap
import time

'''
submitting plots for all the individual observation plots. This submits condor jobs


'''


P = ap.ArgumentParser(description='making individual observation plots on condor',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)


P.add_argument('source', action='store',
           default=None, help='name of source')
P.add_argument('year',action='store', default=None,
           help='the year to gather data for')

args = P.parse_args()

source=args.source
year=args.year

#directories to be modified
plotpath='/scratch/dyang18/dq_plots/'
datatarpath='/spt/user/dyang18/tar/'+str(year)+'/'+str(year)+'_'+source+'.tar.gz'
bolotarpath='/spt/user/dyang18/tar/'+str(year)+'/'+str(year)+'_'+source+'BolometerProperties.tar.gz'
script='/home/dyang18/spt3g/spt3g_software/scratch/dyang18/'+source+'PlotCondor.py'
log_root='/scratch/dyang18'
grid_proxy='/home/dyang18/dyang'
outputpath='/spt/user/dyang18/output/'

datatar = tarfile.open(datatarpath)
bolotar=  tarfile.open(bolotarpath)

dat=[member.name for member in datatar.getmembers()]
bolo=[member.name for member in bolotar.getmembers()]
bolo=[x.split('Bolo')[0] for x in bolo]
#storing the available data
adata=[x for x in bolo if x+'.g3' in dat]

#aplots stores a list of availble observation numbers and date range
aplots= [x.split('/')[-1] for x in glob(plotpath+'*')]
aplots= [ x.split('_')[0] for x in aplots]

submit = [x for x in adata if x not in aplots]



#section for developing the date range to make plots for

def split(x,n):
    return [x[i:i+n] for i in range(0, len(x), n)]
# splitting up the list of files and then submit them

final=split(submit, 50)

#Submission section for individual


input_files=[datatarpath,bolotarpath]
output_files=[source+".tar.gz"]


for i in range(len(final)):
    f=final[i]

    args=[ '-o' +' ' +' '.join(f)]
    jobname='condor_{}_batch_{}'.format(time.strftime('%Y%m%d_%H%M%S'),str(i))
    output_root=outputpath+jobname
    i+=1
    cluster.condor_submit(script=script, args=args, caller='python', jobname=jobname,
                      log_root=log_root, output_root=output_root, input_files=input_files,
                      output_files=output_files, grid_proxy=grid_proxy, aux_input_files=[],
                      aux_output_files=[], spt3g_env=True, python3=False,
                      user_code='', requirements=None, request_cpus=1,
                      request_memory=4*core.G3Units.GB,
                      request_disk=6*core.G3Units.GB, queue=1,
                      when_to_transfer_output='ON_EXIT', create_only=False,
                      globus_uri='gsiftp://gridftp.grid.uchicago.edu:2811/cephfs')
