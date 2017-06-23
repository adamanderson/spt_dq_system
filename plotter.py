from spt3g import core, dfmux, calibration
import matplotlib
matplotlib.use('Agg') # don't use X backend on grid 
import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import os
import tarfile


def d1_plots(obsid, data_frame, varname, bolo_props, plot_dir):
    bands = np.unique([bolo_props[bname].band for bname in bolo_props.keys()])
    wafers = np.unique([bolo_props[bname].wafer_id for bname in bolo_props.keys()])
    data = data_frame[varname]
    
    # plot data by band
    f_band = plt.figure()
    d = np.array([data[bname] for bname in data.keys()])
    d = d[np.isfinite(d)]
    bdata = np.array([bolo_props[bname].band for bname in data.keys()])
    bdata = bdata[np.isfinite(d)]
    for band in bands:
        if varname in plot_ranges.keys():
            histbins = np.linspace(plot_ranges[varname][0], plot_ranges[varname][1], 101)
        else:
            histbins = np.linspace(np.min(d[(bdata==band)]), np.max(d[(bdata==band)]))
        plt.hist(d[(bdata==band)],
                 bins=histbins, alpha=0.7,
                 histtype='stepfilled', label=str(band/core.G3Units.GHz))
    if varname in plot_ranges.keys():
        plt.xlim(plot_ranges[varname])
    plt.legend()
    plt.xlabel(varname)
    plt.savefig('{}/{}_{}_band.png'.format(plot_dir, obsid, varname))
    plt.close(f_band)

    # plot data by wafer
    f_wafer = plt.figure()
    wdata = np.array([bolo_props[bname].wafer_id for bname in data.keys()])
    wdata = wdata[np.isfinite(d)]
    for wafer in wafers:
        wafer_data = np.array([data[bname] for bname in data.keys()
                              if bolo_props[bname].wafer_id == wafer])
        if varname in plot_ranges.keys():
            histbins = np.linspace(plot_ranges[varname][0], plot_ranges[varname][1], 101)
        else:
            histbins = np.linspace(np.min(d[(bdata==band)]), np.max(d[(bdata==band)]))
        plt.hist(wafer_data[np.isfinite(wafer_data)], bins=histbins, 
                 histtype='step', label=wafer)
    if varname in plot_ranges.keys():
        plt.xlim(plot_ranges[varname])
    plt.legend()
    plt.xlabel(varname)
    plt.savefig('{}/{}_{}_wafer.png'.format(plot_dir, obsid, varname))
    plt.close(f_wafer)
    
def d2_plots(xdata, ydata):
    return


def process_RCW38(data_fname, props_fname, plot_dir):
    data_file = core.G3File(data_fname)
    data_frame = data_file.next()
    props_file = core.G3File(props_fname)
    bolo_props = props_file.next()['NominalBolometerProperties']
    
    for var in data_frame.keys():
        d1_plots(os.path.basename(data_fname).rstrip('.g3'), data_frame, var, bolo_props, plot_dir)

def process_calibrator(data_fname, props_fname, plot_dir):
    return


def process_elnod(data_fname, props_fname, plot_dir):
    return


def process_CenA(data_fname, props_fname, plot_dir):
    return


def process_saturn(data_fname, props_fname, plot_dir):
    return


def process_maps(data_fname, props_fname, plot_dir):
    return


def check_plots():
    # this function should check whether a plots directory exists and is 
    # "complete"
    return


# mapping of observation type to processing function
processing_funcs = {'RCW38': process_RCW38,
                    'RCW38-pixelraster': process_RCW38,
                    'calibrator': process_calibrator,
                    'elnod': process_elnod,
                    'CenA': process_CenA,
                    'saturn': process_saturn
                    }

plot_ranges = {'RCW38FluxCalibration':(-2e2, 2e2),
               'PointingOffsetX':(-3e-2, 3e-2),
               'PointingOffsetY':(-3e-2, 3e-2)}


if __name__ == '__main__':
    P0 = ap.ArgumentParser(description='Make data-quality diagnostic plots for SPT.',
                           formatter_class=ap.RawTextHelpFormatter)
    P0.add_argument('processinglist', default=None, action='store', help='Text file '
                    'containing two columns: first column is the type of '
                    'observation to process (e.g. elnod, RCW38, etc.), second '
                    'column is the path of the file.')
    P0.add_argument('plotdir', default=None, action='store', help='Path where plots '
                    'should be stored. This directory will be created if it does '
                    'not already exist.')
    args = P0.parse_args()

    # make output directory for plots if necessary
    if not os.path.exists(args.plotdir):
        os.makedirs(args.plotdir)

    pathlist = np.loadtxt(args.processinglist, dtype=str, ndmin=2)
    for obstype, data_fname, props_fname in pathlist:
        print('Processing {}'.format(data_fname))
        if obstype in processing_funcs.keys():
            processing_funcs[obstype]('{}/{}'.format(os.getcwd(),data_fname),
                                      '{}/{}'.format(os.getcwd(),props_fname),
                                      args.plotdir)

    with tarfile.open('{}.tar.gz'.format(args.plotdir), "w:gz") as tar:
        tar.add(args.plotdir, arcname=os.path.basename(args.plotdir))
    print(os.listdir(os.getcwd()))
