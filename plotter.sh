eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`
alias 3g=/home/adama/SPT/spt3g_software/build/env-shell.sh
export X509_USER_PROXY=$PWD/my_proxy

pathlistfile = "$1"

mkdir data
mkdir plots

while read -r line
do
    path="$line"
    globus-url-copy -vb gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/$path file://$PWD/data
done < "$pathlistfile"

# NOTE: 'pathlistfile' below does not conform to the format expected by plotter.py,
# plotter.py expects a format that is essentially inconsistent with portable/grid usage
3g python plotter.py pathlistfile $PWD/plots

tar cvfz plots.tar.gz plots
globus-url-copy -vb file://$PWD/plots.tar.gz gsiftp://gridftp.grid.uchicago.edu:2811/cephfs/spt/user/adama
