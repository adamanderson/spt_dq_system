from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
from spt3g.std_processing import obsid_to_g3time
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pickle as pickle
import argparse as ap
import glob
import datetime
import os.path
import shutil
from functools import reduce
import operator

P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                          help='Function to perform. For help, call: '
                          '%(prog)s %(metavar)s -h')

timenow = datetime.datetime.now()
default_mintime = datetime.datetime(timenow.year, timenow.month, timenow.day - (timenow.weekday()+1))

S0 = S.add_parser('rebuild', help='Rebuild the pickle files from scratch.',
                  formatter_class=ap.ArgumentDefaultsHelpFormatter)
S0.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
S0.add_argument('bolodatapath', action='store', default=None,
                help='Path to bolometer data (for bolometer properties.')
S0.add_argument('outdir', action='store', default=None,
                help='Path in which to store output data.')
S0.add_argument('--min-time', action='store', default=default_mintime.strftime('%Y%m%d'),
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S0.add_argument('--max-time', action='store',
                default=timenow.strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')

S1 = S.add_parser('update', help='Updates the pickle file data skim when the '
                  'autoprocessing file timestamp is newer than the stored '
                  'timestamp or when there is no stored data for an '
                  'autoprocessing data file.',
                  formatter_class=ap.ArgumentDefaultsHelpFormatter)
S1.add_argument('caldatapath', action='store', default=None,
                help='Path to calibration data to skim.')
S1.add_argument('bolodatapath', action='store', default=None,
                help='Path to bolometer data (for bolometer properties.')
S1.add_argument('outdir', action='store', default=None,
                help='Path in which to store output data.')
S1.add_argument('--min-time', action='store', default=default_mintime.strftime('%Y%m%d'),
                help='Minimum time of observations to skim. Format: YYYYMMDD')
S1.add_argument('--max-time', action='store',
                default=timenow.strftime('%Y%m%d'),
                help='Maximum time of observations to skim. Format: YYYYMMDD')
args = P0.parse_args()


# functions that define splits
def select_band(boloprops, bolo, band):
    return boloprops[bolo].band / core.G3Units.GHz == band

def select_wafer(boloprops, bolo, wafer):
    return boloprops[bolo].wafer_id == wafer

selector_dict = {('w172', 90): (select_wafer, select_band),
                 ('w172', 150): (select_wafer, select_band),
                 ('w172', 220): (select_wafer, select_band),
                 ('w174', 90): (select_wafer, select_band),
                 ('w174', 150): (select_wafer, select_band),
                 ('w174', 220): (select_wafer, select_band),
                 ('w176', 90): (select_wafer, select_band),
                 ('w176', 150): (select_wafer, select_band),
                 ('w176', 220): (select_wafer, select_band),
                 ('w177', 90): (select_wafer, select_band),
                 ('w177', 150): (select_wafer, select_band),
                 ('w177', 220): (select_wafer, select_band),
                 ('w180', 90): (select_wafer, select_band),
                 ('w180', 150): (select_wafer, select_band),
                 ('w180', 220): (select_wafer, select_band),
                 ('w181', 90): (select_wafer, select_band),
                 ('w181', 150): (select_wafer, select_band),
                 ('w181', 220): (select_wafer, select_band),
                 ('w187', 90): (select_wafer, select_band),
                 ('w187', 150): (select_wafer, select_band),
                 ('w187', 220): (select_wafer, select_band),
                 ('w188', 90): (select_wafer, select_band),
                 ('w188', 150): (select_wafer, select_band),
                 ('w188', 220): (select_wafer, select_band),
                 ('w201', 90): (select_wafer, select_band),
                 ('w201', 150): (select_wafer, select_band),
                 ('w201', 220): (select_wafer, select_band),
                 ('w203', 90): (select_wafer, select_band),
                 ('w203', 150): (select_wafer, select_band),
                 ('w203', 220): (select_wafer, select_band)}
wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180', 'w181', 'w187', 'w188', 'w201', 'w203']


# functions that define quantities to be saved
def compute_median(frame, datakey, select_dict):
    data_list = np.array([frame[datakey][bolo] 
                          for bolo in frame[datakey].keys()])
    data_on_selection = {}

    if len(frame[datakey]) == 0:
        return np.array([])

    for select_values, f_select_list in selector_dict.items():
        selection = np.array([True for bolo in frame[datakey].keys()])
        
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # compute the data selection
        for select_val, f_select in zip(select_values, f_select_list):
            selection = np.array([f_select(boloprops, bolo, select_val)
                                  for bolo in frame[datakey].keys()]) & selection
        
        # get data that satisfies the selection and compute median
        reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
            np.median(data_list[selection][np.isfinite(data_list[selection])])

    return data_on_selection


def compute_nalive(frame, datakey, select_dict, sn_threshold):
    data_list = np.array([frame[datakey][bolo] 
                          for bolo in frame[datakey].keys()])
    nalive_on_selection = {}

    for select_values, f_select_list in selector_dict.items():
        selection = np.array([True for bolo in frame[datakey].keys()])
        
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, nalive_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], nalive_on_selection)[keylist[-1]] = {}

        # compute the data selection
        for select_val, f_select in zip(select_values, f_select_list):
            selection = np.array([f_select(boloprops, bolo, select_val)
                                  for bolo in frame[datakey].keys()]) & selection
        
        # get data that satisfies the selection and compute # alive bolos
        data_on_selection = data_list[selection][np.isfinite(data_list[selection])]
        reduce(operator.getitem, select_values[:-1], nalive_on_selection)[select_values[-1]] = \
            len(data_on_selection[data_on_selection>sn_threshold])

    return nalive_on_selection
    

def median_cal_sn_4Hz(frame, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponseSN', selector_dict)

def median_cal_response_4Hz(frame, selector_dict):
    if 'CalibratorResponse' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponse', selector_dict)

def alive_bolos_cal_4Hz(frame, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_nalive(frame, 'CalibratorResponseSN', selector_dict, 10)

def median_elnod_sn_slope(frame, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_median(frame, 'ElnodSNSlopes', selector_dict)

def median_elnod_iq_phase_angle(frame, selector_dict):
    if 'ElnodEigenvalueDominantVectorQ' not in frame.keys() and \
       'ElnodEigenvalueDominantVectorI' not in frame.keys():
        return None

    newframe = core.G3Frame()
    newframe['ElnodPhaseAngle'] = core.G3MapDouble()
    for bolo in frame['ElnodEigenvalueDominantVectorQ'].keys():
        if frame['ElnodEigenvalueDominantVectorI'][bolo] != 0 and \
                frame['ElnodEigenvalueDominantVectorQ'][bolo] != 0:
            newframe['ElnodPhaseAngle'][bolo] = 180/np.pi * \
                np.arctan(frame['ElnodEigenvalueDominantVectorQ'][bolo] / \
                          frame['ElnodEigenvalueDominantVectorI'][bolo])

    return compute_median(newframe, 'ElnodPhaseAngle', selector_dict)

def alive_bolos_elnod(frame, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_nalive(frame, 'ElnodSNSlopes', selector_dict, 20)

def median_rcw38_fluxcal(frame, selector_dict):
    if 'RCW38FluxCalibration' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38FluxCalibration', selector_dict)

def median_rcw38_intflux(frame, selector_dict):
    if 'RCW38IntegralFlux' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38IntegralFlux', selector_dict)

def rcw38_sky_transmission(frame, selector_dict):
    if 'RCW38SkyTransmission' not in frame.keys():
        return None

    data_on_selection = {}

    for select_values, f_select_list in selector_dict.items():
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # get data that satisfies the selection and compute median
        if str(select_values[-1]) in frame['RCW38SkyTransmission'].keys():
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
                frame['RCW38SkyTransmission'][str(select_values[-1])]

    return data_on_selection

    
function_dict = {'RCW38':             {'RCW38SkyTransmission': rcw38_sky_transmission},
                 'RCW38-pixelraster': {'MedianRCW38FluxCalibration': median_rcw38_fluxcal,
                                       'MedianRCW38IntegralFlux': median_rcw38_intflux},
                 'calibrator':        {'MedianCalSN_4Hz': median_cal_sn_4Hz,
                                       'MedianCalResponse_4Hz': median_cal_response_4Hz,
                                       'AliveBolosCal_4Hz': alive_bolos_cal_4Hz},
                 'elnod':             {'MedianElnodIQPhaseAngle': median_elnod_iq_phase_angle,
                                       'MedianElnodSNSlopes': median_elnod_sn_slope,
                                       'AliveBolosElnod': alive_bolos_elnod}}

                 
# loop over data by week
dt_mintime = dt.datetime(year=int(args.min_time[:4]),
                         month=int(args.min_time[4:6]),
                         day=int(args.min_time[6:8]))
dt_maxtime = dt.datetime(year=int(args.max_time[:4]),
                         month=int(args.max_time[4:6]),
                         day=int(args.max_time[6:8]))
d_to_next_week = dt.timedelta(days = 7 - dt_mintime.weekday())
date_boundaries = [dt_mintime]
next_day = dt_mintime + d_to_next_week
while next_day < dt_maxtime:
    date_boundaries.append(next_day)
    next_day = next_day + dt.timedelta(days=7)
date_boundaries.append(dt_maxtime)

# delete the full output directory tree if we are in rebuild mode
if args.mode == 'rebuild' and os.path.exists(args.outdir):
    shutil.rmtree('{}'.format(args.outdir))
    os.mkdir('{}'.format(args.outdir))

for mindate, maxdate in zip(date_boundaries[:-1], date_boundaries[1:]):
    # convert min/max times in observation IDs that we can compare with filenames
    min_obsid = time_to_obsid(core.G3Time('{}_000000'.format(mindate.strftime('%Y%m%d'))))
    max_obsid = time_to_obsid(core.G3Time('{}_000000'.format(maxdate.strftime('%Y%m%d'))))
    
    # make the subdirectory based on time
    outdir = '{}/{}_cached_dq_plots'.format(args.outdir, mindate.strftime('%Y%m%d'))
    if args.mode == 'rebuild':
        os.mkdir(outdir)
        data = {}
    if args.mode == 'update':
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if os.path.exists(os.path.join(outdir, 'data_cache.pkl')):
            data = pickle.load(open(os.path.join(outdir, 'data_cache.pkl'), 'rb'))
        else:
            data = {}

    was_data_updated = False

    for source, quantities in function_dict.items():
        calfiles = glob.glob('{}/calibration/{}/*g3'.format(args.caldatapath, source))
        files_to_parse = [fname for fname in calfiles if int(os.path.splitext(os.path.basename(fname))[0]) >= min_obsid and \
                              int(os.path.splitext(os.path.basename(fname))[0]) <= max_obsid]

        print('Analyzing source: {}'.format(source))

        if source not in data.keys():
            data[source] = {}
        for fname in files_to_parse:
            obsid = os.path.splitext(os.path.basename(fname))[0]
            print('observation: {}'.format(obsid))

            if obsid not in data[source].keys() or \
                    data[source][obsid]['timestamp'] != os.path.getctime(fname):
                was_data_updated = True
                data[source][obsid] = {'timestamp': os.path.getctime(fname)}
                d = [fr for fr in core.G3File(fname)]

                boloprops = [fr for fr in core.G3File('{}/{}/{}/nominal_online_cal.g3' \
                                                          .format(args.bolodatapath, \
                                                                      source, \
                                                                      obsid))] \
                                                                      [0]["NominalBolometerProperties"]

                for quantity_name in function_dict[source]:
                    func_result = function_dict[source][quantity_name](d[0], selector_dict)
                    if func_result:
                        data[source][obsid][quantity_name] = func_result

    with open('{}/data_cache.pkl'.format(outdir), 'wb') as f:
        pickle.dump(data, f)


    def plot_median_cal_sn_4Hz(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['calibrator']
                      if 'MedianCalSN_4Hz' in data['calibrator'][obsid]]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]:                
                median_calSN = np.array([data['calibrator'][obsid]['MedianCalSN_4Hz'][wafer][band]
                                         for obsid in data['calibrator']
                                         if 'MedianCalSN_4Hz' in data['calibrator'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_calSN)],
                         median_calSN[np.isfinite(median_calSN)],
                         'o', label='{} GHz'.format(band))

                if len(median_calSN[np.isfinite(median_calSN)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0, 250])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median calibrator S/N')
            plt.title('4.0 Hz calibrator S/N ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_cal_sn_4Hz_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_cal_response_4Hz(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['calibrator'] 
                      if 'MedianCalResponse_4Hz' in data['calibrator'][obsid]]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]:
                median_cal = np.array([data['calibrator'][obsid]['MedianCalResponse_4Hz'][wafer][band]
                                       for obsid in data['calibrator']
                                       if 'MedianCalResponse_4Hz' in data['calibrator'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_cal)],
                         1e15*median_cal[np.isfinite(median_cal)],
                         'o', label='{} GHz'.format(band))

                if len(median_cal[np.isfinite(median_cal)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0, 5])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median calibrator response [fW]')
            plt.title('4.0 Hz calibrator response ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_cal_response_4Hz_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_alive_bolos_cal_4Hz(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['calibrator'] 
                      if 'AliveBolosCal_4Hz' in data['calibrator'][obsid]]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]:
                n_alive_bolos = np.array([data['calibrator'][obsid]['AliveBolosCal_4Hz'][wafer][band]
                                          for obsid in data['calibrator']
                                          if 'AliveBolosCal_4Hz' in data['calibrator'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(n_alive_bolos)],
                         n_alive_bolos[np.isfinite(n_alive_bolos)],
                         'o', label='{} GHz'.format(band))

                if len(n_alive_bolos[np.isfinite(n_alive_bolos)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0, 600])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('number of alive bolos')
            plt.title('Number of bolos with calibrator S/N > 10 at 4.0 Hz ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/alive_bolos_cal_4Hz_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_elnod_sn(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['elnod']]
            f = plt.figure(figsize=(8,6))
            
            for band in [90, 150, 220]:
                median_elnodSN = np.array([data['elnod'][obsid]['MedianElnodSNSlopes'][wafer][band] for obsid in data['elnod']])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_elnodSN)],
                         median_elnodSN[np.isfinite(median_elnodSN)],
                         'o', label='{} GHz'.format(band))

                if len(median_elnodSN[np.isfinite(median_elnodSN)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0, 2000])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median elnod S/N')
            plt.title('Elnod S/N ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_elnod_sn_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_elnod_iq_phase(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['elnod']]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]: 
                median_elnod_iq = np.array([data['elnod'][obsid]['MedianElnodIQPhaseAngle'][wafer][band]
                                   for obsid in data['elnod']
                                   if 'MedianElnodIQPhaseAngle' in data['elnod'][obsid].keys() and \
                                       data['elnod'][obsid]['MedianElnodIQPhaseAngle'][wafer][band] != None])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianElnodIQPhaseAngle' in data['elnod'][obsid].keys() and \
                              data['elnod'][obsid]['MedianElnodIQPhaseAngle'][wafer][band] != None]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_elnod_iq)],
                         median_elnod_iq[np.isfinite(median_elnod_iq)],
                         'o', label='{} GHz'.format(band))

                if len(median_elnod_iq[np.isfinite(median_elnod_iq)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([-90, 90])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median elnod IQ phase [deg]')
            plt.title('Elnod IQ phase angle ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_elnod_iq_phase_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_alive_bolos_elnod(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['elnod']]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]:
                alive_bolos_elnod = np.array([data['elnod'][obsid]['AliveBolosElnod'][wafer][band] for obsid in data['elnod']])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(alive_bolos_elnod)],
                         alive_bolos_elnod[np.isfinite(alive_bolos_elnod)],
                         'o', label='{} GHz'.format(band))

                if len(alive_bolos_elnod[np.isfinite(alive_bolos_elnod)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0, 600])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('number of alive bolos')
            plt.title('Number of bolos with elnod S/N>20 ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/alive_bolos_elnod_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_rcw38_fluxcal(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['RCW38-pixelraster']]
            f = plt.figure(figsize=(8,6))

            for band in [90, 150, 220]:
                median_rcw38 = np.array([data['RCW38-pixelraster'][obsid]['MedianRCW38FluxCalibration'][wafer][band]
                                for obsid in data['RCW38-pixelraster']])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_rcw38)],
                         median_rcw38[np.isfinite(median_rcw38)],
                         'o', label='{} GHz'.format(band))

                if len(median_rcw38[np.isfinite(median_rcw38)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([-100, 0])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median RCW38 flux calibration')
            plt.title('RCW38 Flux Calibration ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_rcw38_fluxcal_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_rcw38_intflux(data):
        for wafer in wafer_list:   
            obsids = [obsid for obsid in data['RCW38-pixelraster']]
            f = plt.figure(figsize=(8,6))
            
            for band in [90, 150, 220]:
                median_rcw38 = np.array([data['RCW38-pixelraster'][obsid]['MedianRCW38IntegralFlux'][wafer][band]
                                for obsid in data['RCW38-pixelraster']])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_rcw38)],
                         median_rcw38[np.isfinite(median_rcw38)],
                         'o', label='{} GHz'.format(band))

                if len(median_rcw38[np.isfinite(median_rcw38)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([2e-7, 7e-7])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median RCW38 integral flux')
            plt.title('RCW38 Integral Flux ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_rcw38_intflux_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_rcw38_sky_transmission(data):
        for wafer in wafer_list:   
            obsids = [obsid for obsid in data['RCW38']]
            f = plt.figure(figsize=(8,6))
            
            for band in [90, 150, 220]:
                rcw38_skytrans = np.array([data['RCW38'][obsid]['RCW38SkyTransmission'][wafer][band]
                                  for obsid in obsids
                                  if type(data['RCW38'][obsid]['RCW38SkyTransmission'][wafer][band])==float])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                              for obsid in obsids
                              if type(data['RCW38'][obsid]['RCW38SkyTransmission'][wafer][band])==float]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)
                
                plt.plot(datenums[np.isfinite(rcw38_skytrans)],
                         rcw38_skytrans[np.isfinite(rcw38_skytrans)],
                         'o', label='{} GHz'.format(band))

                if len(rcw38_skytrans[np.isfinite(rcw38_skytrans)])>0:
                    xfmt = mdates.DateFormatter('%m-%d %H:%M')
                    plt.gca().xaxis.set_major_formatter(xfmt)
                    plt.xticks(rotation=25)
                    plt.ylim([0.85, 1.25])
                    plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('RCW38 sky transmission')
            plt.title('RCW38 Sky Transmission ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/rcw38_sky_transmission_{}.png'.format(outdir, wafer))
            plt.close()


    # only update figures if the underlying data actually changed.
    if was_data_updated:
        # create the plots
        plot_median_cal_sn_4Hz(data)
        plot_median_cal_response_4Hz(data)
        plot_alive_bolos_cal_4Hz(data)
        plot_median_elnod_sn(data)
        plot_median_elnod_iq_phase(data)
        plot_alive_bolos_elnod(data)
        plot_median_rcw38_fluxcal(data)
        plot_median_rcw38_intflux(data)
        plot_rcw38_sky_transmission(data)

        # create symlink from latest data directory to current
        symlinkname = '{}/current'.format(args.outdir)
        os.symlink(outdir, '{}/temp'.format(args.outdir))
        os.rename('{}/temp'.format(args.outdir), '{}/current'.format(args.outdir))
