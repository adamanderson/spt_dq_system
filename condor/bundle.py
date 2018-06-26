import tarfile
import numpy
from glob import glob
import argparse as ap
import os.path
from astropy.time import Time
from spt3g.std_processing import obsid_to_g3time
import sys

# this script creates tarballs for data for a specific source in a given date (2018 or 2018-01)
outputroot="/spt/user/dyang18/tar/"

P = ap.ArgumentParser(description='Taring files for condor',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)


P.add_argument('source', action='store',
           default=None, help='name of source')
P.add_argument('bp', action='store',
           default=None, help='True False, whether to also process bolometer property files')
P.add_argument('-d','--daterange', nargs='+', default=None,
           help='the daterange to gather data for')

args = P.parse_args()
obstype = str(args.source)
root_folder ='/spt/user/production/calibration/'
input_file_lst = glob(root_folder +obstype+'/'+'*.g3')

if args.daterange:
    directory= outputroot +''.join(args.daterange)+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    tar = tarfile.open(outputroot +''.join(args.daterange) +'/'+ ''.join(args.daterange)+'_'+obstype+".tar.gz", "w:gz")
    for fil in input_file_lst:


        # processing the time format
        f= fil.split('/')[-1]
        obsid= f.split('.')[0]

        t= obsid_to_g3time(obsid)
        t= Time(t.mjd,format='mjd').iso

        for d in args.daterange:

            if str(d) in t:

                tar.add(fil,arcname=fil.split('/')[-1],recursive =False)


    tar.close()

    if args.bp:
        print("Taring the corresponding bolometer properties")
        tar = tarfile.open(outputroot +''.join(args.daterange) +'/'+ ''.join(args.daterange)+'_'+obstype+ "BolometerProperties" +".tar.gz", "w:gz")
        root_folder ='/spt/data/bolodata/downsampled/'+obstype+'/'
        input_file_lst = glob(root_folder+'*')
        for fil in input_file_lst:
            obsid= int((fil.split('/'))[-1])
            t= obsid_to_g3time(obsid)
            t= Time(t.mjd,format='mjd').iso
            for d in args.daterange:

                if str(d) in t and os.path.exists(fil+'/nominal_online_cal.g3'):
                    tar.add(fil+'/nominal_online_cal.g3',arcname=str(obsid)+'BolometerProperties.g3')



        print("finished tarring the " + obstype + "files")

        tar.close()
