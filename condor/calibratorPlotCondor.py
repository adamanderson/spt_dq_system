import matplotlib
matplotlib.use('Agg')
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
import scipy.stats
import os, numpy
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import cm
import pandas as pd
from spt3g.std_processing import obsid_to_g3time
from glob import glob
import tarfile
from astropy.time import Time
# Script used on OSG to generate plots for calibrator
# Pass in the needed parameters using argparse

P = ap.ArgumentParser(description='Various Calibrator Plots',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('-o','--obsid', nargs='+', 
               help='list of observation number', default=None)
P.add_argument('-d','--daterange', nargs='+', default=None,
           help='the daterange to make plots for')
args = P.parse_args()
year='2018'
#extract all the observation files which were packed into tar files
dpi_setting=200
tar = tarfile.open(year+'_'+'calibrator.tar.gz', "r:gz")
tar.extractall()
tar.close()
tar = tarfile.open(year+'_'+'calibratorBolometerProperties.tar.gz', "r:gz")
tar.extractall()
tar.close()

# obsp keeps track of a list of parameter that we want to plot for
obsp=['CalibratorResponse', 'CalibratorResponseSN']
obstype='calibrator'
# Helper functions for stats things
def finiteStats(a, remove_percent=3):
    '''
    Throw out top and bottom remove_percent% of outliers (unless boolean)
    return (min, max, mean, sigma) of the finite
    '''
    # Throw out nan
    finite = np.isfinite(a)
    af = a[finite]  # Also flattens
    if len(af)==0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    #amed  = np.median(af)
    if af.dtype == bool:
        remove_percent = 0  # Don't remove boolean outliers
    amid  = af[ (af >= np.percentile(af,remove_percent)) * 
                (af <= np.percentile(af,100-remove_percent)) ]
    # Note it must be >= and <= or else an array of 0s and 1s would come back empty
    mina   = np.min(af)
    maxa   = np.max(af)
    mean   = np.mean(amid)  # Note empty arrays return nan
    sigma  = np.std(amid)
    median = np.median(amid)
    return (mina, maxa, mean, sigma, median)

def autoRange(a, nSigma=4.1, remove_percent=3):
    '''
    Throw out top and bottom remove_percent% of outliers
    Return (vmin, vmax, ticks) were
    (vmin, vmax) ~= (mean-nSigma*sigma, mean+nSigma*sigma)
    unless either side of this exceeds the actual min/max,
    in which case return that.
    The approximatly-equal uses matplotlib's ticker.AutoLocator to round to nice values, 
    since the histogram will extend the limits to these anyway
    NOTE: Unfortunately, vmin and vmax are sometimes slightly too big: 1.00000000000000008882
          Which is why clip() and the Histograming part of makePlots() must shrink this range by epsilon
    for bool or and int types with a small range, lock the min and max exactly to take full advantage of the colors.
    '''
    def isbool(a):
        return type(a)==np.ndarray and a.dtype.kind=='b'
    def issmallint(a):
        return type(a)==np.ndarray and a.dtype.kind=='i' and (maxa-mina)<=20
    if isbool(a):  # a type of bool
        low, high = (0, 1)
    else:
        (mina, maxa, mean, sigma, median) = finiteStats(a, remove_percent)  # returns all nan if there are no finite
        if not np.isfinite([mina, maxa, mean, sigma, median]).all():
            low, high = (0, 1)
        elif issmallint(a): # a type of int with low values
            low, high = (mina, maxa)  # We want the full range in this case.
        else:
            low  = max(mina, mean-nSigma*sigma)
            high = min(maxa, mean+nSigma*sigma)
    # Use matplotlib's ticker.AutoLocator to round
    al = matplotlib.ticker.AutoLocator()
    ticks = al.bin_boundaries(low, high)
    if issmallint(a):
        (vmin, vmax) = (low, high)  # Really take exact integer range.
    else:
        (vmin, vmax) = ticks[[0,-1]]   # pick out AutoLocator's first and last tick
    return [vmin, vmax]
def StackedHistogram(request):
    '''
    Makes a stacked histogram for a parameter given in request

    '''
    obsid=request['observation']
    varname=request['parameter']
    # setting up the appopriate variables
    for fr in core.G3File(request['input']):
        data_frame=fr
        break
        
    for fr in core.G3File(request['boloprop']):
        bolo_props=fr['NominalBolometerProperties']

    # Grabbing all the correct data files
    bands=[900.0,1500.0,2200.0]

    #temporarily set 3 frequencies because some bolos are giving negative
    #frequency which is causing issues. Commented out is the original code

    #bands = np.unique([bolo_props[bname].band for bname in bolo_props.keys()])
    #bands.sort()

    wafers = np.unique([bolo_props[bname].wafer_id for bname in bolo_props.keys()])
    wafers.sort()
    #Check what wafers are available and make a list of them
    try:
        data={k: v for k, v in data_frame[varname].items() if k in bolo_props.keys()}
               #{k:v for k,v in data_frame[varname].iteritems() if k in bolo_props.keys()}
    except:
        return
    # converting the data into dictionary in which the bolo ID is the key and the data of interest is
    # the value
    colors = np.linspace(0,1,len(bands)*len(wafers))
    color_list = cm.rainbow(colors).tolist()
    color=[]
    for i in color_list:
        color.append(tuple(i))
    color.reverse()
    #establishing color conventions for the band and wafer combinations
    
    plotrange=autoRange(np.array((data_frame[varname]).values()))
    # generate the numerical ranges for the plot


    if np.isfinite(plotrange[0]) and np.isfinite(plotrange[1]):
        
        
        # check that we don't get nan in plot ranges

        f = plt.figure()
        bindata=[]
        bindatakey=[]
        histbins = np.linspace(plotrange[0], plotrange[1], 30)
        for band in bands:
            for wafer in wafers:
                bindatakey.append(str(band/core.G3Units.GHz)+'/'+str(wafer))
                bindata.append([data[bname] for bname in data.keys() if bolo_props[bname].wafer_id == wafer and bolo_props[bname].band==band])
        
        plt.hist([np.clip(bindata[i],plotrange[0],plotrange[1]) for i in range(len(bindata))], bins=histbins,range=(plotrange[0],plotrange[1]), color=color, stacked =True, label=bindatakey)
        
        plt.legend(bbox_to_anchor=(0.9,1.05),loc=2,prop={'size':7})
        plt.xlim(plotrange[0],plotrange[1])
        plt.title(obstype+' ' + obsid  +' ' + varname)
        #f.savefig(obsid+'_'+'calibrator'+'StackedHistogram'+'.png',bbox_inches='tight')
        plt.savefig(obsid+'_'+'StackedHistogram'+'_'+obstype+'_' + str(varname)+'.png',bbox_inches='tight',dpi=dpi_setting)
        plt.close()
        

def FocalPlane(request):
    '''
    This should be able to be replaced by the existing version of focal plane plot.
    '''
    obsid=request['observation']
    varname=request['parameter'] #CalibratorResponse'
    #this grabs the right calframe
    for fr in core.G3File(request['input']):
        data_frame=fr
        break
    # need to be changed according to how I treat bolometer property file.
    for fr in core.G3File(request['boloprop']):
        bolo_props=fr['NominalBolometerProperties']
    bands = np.unique([bolo_props[bname].band for bname in bolo_props.keys()])
    bands.sort()
    try:
        data={k: v for k, v in data_frame[varname].items() if k in bolo_props.keys()}
    except: 
        return


    #color map
    cmap = plt.get_cmap('plasma_r')
    radius=0.00015
    linewidth=1.3

    for i in range(len(bands)):

        bindata=np.array([np.abs(data[bname]) for bname in data.keys() if bolo_props[bname].band==bands[i]])
        binangle=np.array([bolo_props[bname].pol_angle for bname in data.keys() if bolo_props[bname].band==bands[i]])
        xdata=np.array([bolo_props[bname].x_offset for bname in data.keys() if bolo_props[bname].band==bands[i]])
        ydata=np.array([bolo_props[bname].y_offset for bname in data.keys() if bolo_props[bname].band==bands[i]])

        if len(xdata)>0 and len(ydata)>0 and len(bindata)>0:

            f, ax1 = plt.subplots(2,gridspec_kw = {'height_ratios':[10, 1]},figsize=(8,8))
            plt.title(obstype+' ' + obsid  +' ' + varname)
            pltrange=autoRange(bindata)
            if np.isfinite(pltrange[0]) and np.isfinite(pltrange[1]):

                marginx=(np.nanmax(xdata)- np.nanmin(xdata))*0.2
                marginy=(np.nanmax(ydata)- np.nanmin(ydata))*0.2

                for j in range(len(xdata)):

                    if np.isfinite(bindata[j]):
                        c=(bindata[j]-pltrange[0])/(pltrange[1]-pltrange[0])
                        x_ll, x_ur = xdata[j]-radius*np.sin(binangle[j]), xdata[j]+radius*np.sin(binangle[j])
                        y_ll, y_ur = ydata[j]-radius*np.cos(binangle[j]), ydata[j]+radius*np.cos(binangle[j])
                        if binangle[j]:
                            ax1[0].add_line( matplotlib.lines.Line2D([x_ll, x_ur], [y_ll, y_ur], linewidth=linewidth, color=cmap(c), alpha=1.0) )

                ax1[0].set_xlim([np.nanmin(xdata)-marginx,np.nanmax(xdata)+marginx])
                ax1[0].set_ylim([np.nanmin(ydata)-marginy,np.nanmax(ydata)+marginy])
                norm = matplotlib.colors.Normalize(vmin=pltrange[0], vmax=pltrange[1])
                cb = matplotlib.colorbar.ColorbarBase(ax1[1], cmap=cmap,norm=norm,orientation='horizontal')
                
                plt.savefig(str(obsid)+'_'+str(bands[i]/core.G3Units.GHz)+'GHz_'+'focalplane_'+obstype+'_' + str(varname)+'.png',dpi=dpi_setting)
                
                plt.close()
                
def HistogramOverTime(request):
    fsize=8
    bins=100
    obstype='calibrator'
    obsparameter=request['parameter']
    filelst=request['input']
    obslist=request['observation']
    dater=request['daterange']

    #end of required parameters


    if len(filelst)!=0:
        # grab the right files and observation numbers


        obslist.sort()
        filelst.sort()

        for fr in core.G3File(filelst):
            if obsparameter in fr:
                bolos=fr[obsparameter].keys()
            break

        #pulling out the available bolometers

        frs=[]
        newobslist=[]
        #make an updated list of frames
        for i in range(len(filelst)):
            for fr in core.G3File(filelst[i]):
                if obsparameter in fr:
                    frs.append(fr[obsparameter].values())
                    newobslist.append(obslist[i])

        frs=pd.DataFrame(frs)
        frs=frs[frs!=0]
        frs=frs.as_matrix()
        frs=np.absolute(frs)

        flatfrs=frs.flatten()
        flatfrs=flatfrs[np.logical_not(np.isnan(flatfrs))]
        vmin,vmax=autoRange(flatfrs)
        if np.isfinite(vmin) and np.isfinite(vmax):
            histo_edges=[np.histogram(np.clip(i,vmin,vmax),bins=bins,range=(vmin,vmax)) for i in frs]
            histo_lines = [he[0] for he in histo_edges]
            histomax = np.percentile([max(h) for h in histo_lines], 50)

            step=float(vmax-vmin)/float(bins)


            f=plt.figure(figsize=(fsize,fsize+fsize*len(frs)/bins))
            plt.imshow(histo_lines,interpolation='none',vmin=vmin, vmax=histomax)
            plt.yticks([int(i) for i in np.linspace(0,len(frs)-1,num=len(frs)/10)],[str(obsid_to_g3time(newobslist[j]))[:17] for j in [int(i) for i in np.linspace(0,len(frs)-1,num=len(frs)/10)]],fontsize=8)
            plt.xticks(np.linspace(0,bins,5),[step*i for i in np.linspace(0,bins,5)],rotation='vertical',fontsize=8)
            plt.colorbar(fraction=0.02, pad=0.02, aspect=40)
            plt.xlabel(obsparameter)
            plt.title(obstype+' ' + dater +' ' + obsparameter)
            #f.savefig(str(dater)+'_'+'calibrator'+str(obsparameter)+'HistogramOverTime.png',bbox_inches='tight',dpi=100)
            plt.savefig(str(dater)+'_'+'HistogramOverTime'+'_'+str(obsparameter)+'_'+obstype+'.png',bbox_inches='tight',dpi=dpi_setting)
            
            plt.close()


def SingleWaferArray(request):

    mydpi=dpi_setting
    fsize=3
    obstype='calibrator'

    obsparameter=request['parameter'] #CalibratorResponse'

    filelst=request['input']
    obslist=request['observation']
    dater=request['daterange']

    for fr in core.G3File(request['boloprop']):
        bp=fr['NominalBolometerProperties']
        

    wafers = np.unique([bp[bname].wafer_id for bname in bp.keys()])
    wafers.sort()
    

    if len(filelst)!=0:
        obslist.sort()
        filelst.sort()

        frs=[]
        newobslist=[]
        #make an updated list of frames
        for i in range(len(filelst)):
            d=core.G3File(filelst[i])
            for fr in d:
                if obsparameter in fr:
                    frs.append(fr[obsparameter])
                    newobslist.append(obslist[i])
        for wafer in wafers:
            for fr in core.G3File(request['boloprop']):
                bp=fr['NominalBolometerProperties']

                bolos=[key for key in bp.keys() if bp[key].wafer_id==wafer]
                break


            if len(bolos)<15:
                #print( 'skipped '+ wafer + ' ' + dater)
                continue
        # make and fill in data
            arr=np.empty((len(frs), len(bolos),))
            arr[:] = np.nan
            
            
            for i in range(len(arr)):
                for j in range(len(arr[0])):
                    if bolos[j] in frs[i]:
                        arr[i][j]=(frs[i])[bolos[j]]
                    #if not np.isfinite(arr[i][j]):
                        #arr[i][j]=0
                        
            if len(arr)!=0:
                
                vmin,vmax=autoRange(arr.flatten())
                
                if np.isfinite(vmin) and np.isfinite(vmax):
                    arr= np.clip(arr, vmin,vmax)
                    f=plt.figure(figsize=(fsize*len(bolos)/mydpi,2*fsize+fsize*len(frs)/mydpi),dpi=dpi_setting)
                    plt.imshow(arr,interpolation='none')
                    plt.yticks([int(i) for i in np.linspace(0,len(frs)-1,num=len(frs)/10)],
                    [str(obsid_to_g3time(newobslist[j]))[:-13] for j in 
                    [int(i) for i in np.linspace(0,len(frs)-1,num=len(frs)/10)]],fontsize=7)

                    plt.xticks([int(i) for i in np.linspace(0,len(bolos)-1,num=len(bolos)/10)],
                               [bolos[j] for j in 
                               [int(i) for i in np.linspace(0,len(bolos)-1,num=len(bolos)/10)]],rotation='vertical',fontsize=5)
                    plt.colorbar(fraction=0.03)
                    plt.title(obstype+' ' + dater +' ' + obsparameter+ ' '+ wafer)
                    #plt.savefig(str(dater)+'_'+'calibrator'+wafer+str(obsparameter)+'Array.png',bbox_inches='tight',dpi=100)
                    plt.savefig(str(dater)+'_'+'array'+'_'+str(obsparameter)+'_'+obstype+'_'+wafer+'.png',bbox_inches='tight',dpi=dpi_setting)
                    plt.close()

#new functions
def CalSNCDF(request):
    obsid = request['observation']
    try:
        data = [fr for fr in core.G3File(request['input'])]
        boloprops = [fr for fr in core.G3File(request['boloprop'])] \
                                                  [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['CalibratorResponseSN'][bolo] \
                                           for bolo in data[0]['CalibratorResponseSN'].keys() \
                                           if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "CalibratorResponseSN does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        histvals, edges = np.histogram(cal_dict[band][np.isfinite(cal_dict[band])],
                                       bins=np.linspace(0,400,50))
        ecdf_unnormed = np.cumsum(histvals)
        ecdf_unnormed = np.hstack([0, ecdf_unnormed])
        inv_ecdf_unnormed = np.max(ecdf_unnormed) - ecdf_unnormed
        plt.step(edges, inv_ecdf_unnormed, label='{} GHz'.format(band))
    #plt.yscale('log')
    plt.legend()

    plt.xlabel('calibrator S/N')
    plt.title('cumulative bolometers above a given calibrator S/N\n'
              'for observation {}\n'
              'chop freq. = {:.1f} Hz'.format(request['observation'],
                                              data[0]["CalibratorResponseFrequency"] / core.G3Units.Hz))
    plt.tight_layout()
    plt.savefig(str(obsid)+'_'+'CalSNCDF_'+obstype+'.png', bbox_inches='tight')
    plt.close()

def CalHistogram(request):
    try:
        data = [fr for fr in core.G3File(request['input'])]
        boloprops = [fr for fr in core.G3File(request['boloprop'])] \
                                                  [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([1e15 * data[0]['CalibratorResponse'][bolo] \
                                  for bolo in data[0]['CalibratorResponse'].keys() \
                                  if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "CalibratorResponse does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(-1, 6, 50),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()
    plt.xlim([-1,6])

    plt.xlabel('calibrator response [fW]')
    plt.title('Calibrator response for observation {}\n'
              'chop freq. = {:.1f} Hz'.format(request['observation'],
                                              data[0]["CalibratorResponseFrequency"] / core.G3Units.Hz))
    plt.tight_layout()
    plt.savefig(str(obsid)+'_'+'CalHistogram_'+obstype+'.png', bbox_inches='tight')
    plt.close()




def CalSNHistogram(request):
    try:
        data = [fr for fr in core.G3File(request['input'])]
        boloprops = [fr for fr in core.G3File(request['boloprop'])] \
                                                  [0]["NominalBolometerProperties"]
    except RuntimeError:
        return "Could not find data file."

    bands = [90, 150, 220]
    try:
        cal_dict = {}
        for band in bands:
            cal_dict[band] = np.array([data[0]['CalibratorResponseSN'][bolo] \
                                  for bolo in data[0]['CalibratorResponseSN'].keys() \
                                  if boloprops[bolo].band / core.G3Units.GHz == band])
    except KeyError:
        return "CalibratorResponseSN does not exist for this observation."

    fig = plt.figure()
    for band in bands:
        plt.hist(cal_dict[band][np.isfinite(cal_dict[band])],
                 bins=np.linspace(0,400,50),
                 label='{} GHz'.format(band),
                 histtype='step')
    plt.legend()

    plt.xlabel('calibrator S/N')
    plt.title('Calibrator S/N for observation {}\n'
              'chop freq. = {:.1f} Hz'.format(request['observation'],
                                              data[0]["CalibratorResponseFrequency"] / core.G3Units.Hz))
    plt.tight_layout()
    plt.savefig(str(obsid)+'_'+'CalSNHistogram_'+obstype+'.png', bbox_inches='tight')
    plt.close()

    
    
# setting up the list of available files and observations
Bolo_input_files=glob('*Properties.g3')
all_files=glob('*.g3')
input_files=[item for item in all_files if item not in Bolo_input_files]
input_files= [item for item in input_files if item[:-3]+'BolometerProperties.g3' in Bolo_input_files]
input_files.sort()
t=[str(x)+'-'+y for x in [2018] for y in ['01','02','03','04','05','06','07','08','09','10','11','12']]
#TEST code

if args.daterange is not None:
    for daterange in args.daterange:

        obslist=[]
        for fil in input_files:
            obsid=int(fil.split('.')[0])
            x=obsid_to_g3time(obsid)
            if str(daterange) in Time(x.mjd,format='mjd').iso:
                obslist.append(obsid)

        obslist.sort()
        filelst=[str(obsid)+'.g3' for obsid in obslist]

        for p in obsp:
            if len(filelst)!=0:
                r={}
                r['boloprop']=filelst[0][:-3]+'BolometerProperties.g3'
                r['observation']=obslist
                r['input']=filelst
                r['daterange']=daterange
                r['parameter']=p
                HistogramOverTime(r)
                SingleWaferArray(r)


if args.obsid is not None:
    for obsid in args.obsid:
        r={}
        r['observation']=obsid
        r['input']=str(obsid)+'.g3'
        r['boloprop']=str(obsid)+'BolometerProperties.g3'
        
        for p in obsp:
            
            r['parameter']=p
            StackedHistogram(r)
            FocalPlane(r)
            
        # stuff that is not parameter dependent
        CalSNCDF(r)
        CalSNHistogram(r)
        CalHistogram(r)

#Now tar the outputfiles
outputfilelst=glob('*.png')
'''
tar = tarfile.open(obstype +".tar.gz", "w:gz")
for fil in outputfilelst:
    tar.add(fil)
tar.close()
'''
tar = tarfile.open(obstype+".tar.gz", "w:gz")
for name in outputfilelst:
    if name[:7] in t:
        tar.add(name,arcname='/'+name[:7]+'/'+name)
    else:
        tar.add(name)

tar.close()
