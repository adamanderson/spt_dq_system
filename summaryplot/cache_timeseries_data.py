import matplotlib
matplotlib.use('Agg')

from spt3g import core
from spt3g.std_processing.utils import time_to_obsid
from spt3g.std_processing import obsid_to_g3time
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pickle as pickle
import argparse as ap
import glob
import datetime
import os.path
import shutil
from calibrator import *
from elnod import *
from mat5a import *
from noise import *
from rcw38 import *


P0 = ap.ArgumentParser(description='',
                       formatter_class=ap.ArgumentDefaultsHelpFormatter)
S = P0.add_subparsers(dest='mode', metavar='MODE', title='subcommands',
                          help='Function to perform. For help, call: '
                          '%(prog)s %(metavar)s -h')

timenow = datetime.datetime.now()
dt = datetime.timedelta(-1*(timenow.weekday()+1))
default_mintime = timenow + dt
cache_dir_stub = 'cached_dq_plots'

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
S1.add_argument('--new-plots', action='store_true', default=False,
                help='Force regeneration of new plots, even if data has not changed.')
args = P0.parse_args()


# functions that define splits
def select_band(boloprops, bolo, band):
    try:
        return boloprops[bolo].band / core.G3Units.GHz == band
    except:
        return False

def select_wafer(boloprops, bolo, wafer):
    if wafer == 'all':
        return True
    else:
        try:
            return boloprops[bolo].wafer_id == wafer
        except:
            return False

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
                 ('w203', 220): (select_wafer, select_band),
                 ('w204', 90): (select_wafer, select_band),
                 ('w204', 150): (select_wafer, select_band),
                 ('w204', 220): (select_wafer, select_band),
                 ('w206', 90): (select_wafer, select_band),
                 ('w206', 150): (select_wafer, select_band),
                 ('w206', 220): (select_wafer, select_band),
                 ('all', 90): (select_wafer, select_band),
                 ('all', 150): (select_wafer, select_band),
                 ('all', 220): (select_wafer, select_band)}
wafer_list = ['w172', 'w174', 'w176', 'w177', 'w180',
              'w181', 'w187', 'w188', 'w201', 'w203',
              'w204', 'w206', 'all']
    
function_dict = {'RCW38':             {'RCW38SkyTransmission': rcw38_sky_transmission},
                 'RCW38-pixelraster': {'MedianRCW38FluxCalibration': median_rcw38_fluxcal,
                                       'MedianRCW38IntegralFlux': median_rcw38_intflux},
                 'MAT5A':             {'MAT5ASkyTransmission': mat5a_sky_transmission},
                 'MAT5A-pixelraster': {'MedianMAT5AFluxCalibration': median_mat5a_fluxcal,
                                       'MedianMAT5AIntegralFlux': median_mat5a_intflux},
                 'calibrator':        {'MedianCalSN_4Hz': median_cal_sn_4Hz,
                                       'MedianCalResponse_4Hz': median_cal_response_4Hz,
                                       'AliveBolosCal_4Hz': alive_bolos_cal_4Hz},
                 'elnod':             {'MedianElnodIQPhaseAngle': median_elnod_iq_phase_angle,
                                       'MedianElnodSNSlopes': median_elnod_sn_slope,
                                       'AliveBolosElnod': alive_bolos_elnod},
                 'noise':             {'NEI_0.1Hz_to_0.5Hz': median_nei_01Hz_to_05Hz,
                                       'NEI_1.0Hz_to_2.0Hz': median_nei_1Hz_to_2Hz,
                                       'NEI_3.0Hz_to_5.0Hz': median_nei_3Hz_to_5Hz,
                                       'NEI_10.0Hz_to_15.0Hz': median_nei_10Hz_to_15Hz,
                                       'NEP_0.1Hz_to_0.5Hz': median_nep_01Hz_to_05Hz,
                                       'NEP_1.0Hz_to_2.0Hz': median_nep_1Hz_to_2Hz,
                                       'NEP_3.0Hz_to_5.0Hz': median_nep_3Hz_to_5Hz,
                                       'NEP_10.0Hz_to_15.0Hz': median_nep_10Hz_to_15Hz,
                                       'NET_0.1Hz_to_0.5Hz': median_net_01Hz_to_05Hz,
                                       'NET_1.0Hz_to_2.0Hz': median_net_1Hz_to_2Hz,
                                       'NET_3.0Hz_to_5.0Hz': median_net_3Hz_to_5Hz,
                                       'NET_10.0Hz_to_15.0Hz': median_net_10Hz_to_15Hz}}

                 
# loop over data by week
dt_mintime = datetime.datetime(year=int(args.min_time[:4]),
                         month=int(args.min_time[4:6]),
                         day=int(args.min_time[6:8]))
dt_maxtime = datetime.datetime(year=int(args.max_time[:4]),
                         month=int(args.max_time[4:6]),
                         day=int(args.max_time[6:8]))
d_to_next_week = datetime.timedelta(days = 7 - dt_mintime.weekday())
date_boundaries = [dt_mintime]
next_day = dt_mintime + d_to_next_week
while next_day < dt_maxtime:
    date_boundaries.append(next_day)
    next_day = next_day + datetime.timedelta(days=7)
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
    outdir = '{}/{}_{}'.format(args.outdir, mindate.strftime('%Y%m%d'), cache_dir_stub)
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
                    func_result = function_dict[source][quantity_name](d[0], boloprops, selector_dict)
                    if func_result:
                        data[source][obsid][quantity_name] = func_result

    with open('{}/data_cache.pkl'.format(outdir), 'wb') as f:
        pickle.dump(data, f)


    def plot_median_cal_sn_4Hz(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['calibrator']
                      if 'MedianCalSN_4Hz' in data['calibrator'][obsid]]
            f = plt.figure(figsize=(8,6))

            is_empty = True
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
                    is_empty = False
                    
            if is_empty == False:
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

            is_empty = True
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
                    is_empty = False

            if is_empty == False:
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

            is_empty = True
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
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                if wafer == 'all':
                    plt.ylim([0, 6000])
                else:
                    plt.ylim([0, 600])
                plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('number of alive bolos')
            plt.title('Number of bolos with calibrator S/N > 20 at 4.0 Hz ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/alive_bolos_cal_4Hz_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_elnod_sn(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['elnod']]
            f = plt.figure(figsize=(8,6))
            
            is_empty = True
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
                    is_empty = False

            if is_empty == False:
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

            is_empty = True
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
                    is_empty = False

            if is_empty == False:
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

            is_empty = True
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
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                if wafer == 'all':
                    plt.ylim([0, 6000])
                else:
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

            is_empty = True
            for band in [90, 150, 220]:
                median_rcw38 = np.array([data['RCW38-pixelraster'][obsid]['MedianRCW38FluxCalibration'][wafer][band]
                                         for obsid in obsids
                                         if 'MedianRCW38FluxCalibration' in data['RCW38-pixelraster'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianRCW38FluxCalibration' in data['RCW38-pixelraster'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_rcw38)],
                         median_rcw38[np.isfinite(median_rcw38)],
                         'o', label='{} GHz'.format(band))

                if len(median_rcw38[np.isfinite(median_rcw38)])>0:
                    is_empty = False

            if is_empty == False:
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
            
            is_empty = True
            for band in [90, 150, 220]:
                median_rcw38 = np.array([data['RCW38-pixelraster'][obsid]['MedianRCW38IntegralFlux'][wafer][band]
                                         for obsid in obsids
                                         if 'MedianRCW38IntegralFlux' in data['RCW38-pixelraster'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianRCW38IntegralFlux' in data['RCW38-pixelraster'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_rcw38)],
                         median_rcw38[np.isfinite(median_rcw38)],
                         'o', label='{} GHz'.format(band))

                if len(median_rcw38[np.isfinite(median_rcw38)])>0:
                    is_empty = False

            if is_empty == False:
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
            
            is_empty = True
            for band in [90, 150, 220]:
                rcw38_skytrans = np.array([data['RCW38'][obsid]['RCW38SkyTransmission'][wafer][band]
                                           for obsid in obsids
                                           if 'RCW38SkyTransmission' in data['RCW38'][obsid]])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                              for obsid in obsids
                              if 'RCW38SkyTransmission' in data['RCW38'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(rcw38_skytrans)],
                         rcw38_skytrans[np.isfinite(rcw38_skytrans)],
                         'o', label='{} GHz'.format(band))

                if len(rcw38_skytrans[np.isfinite(rcw38_skytrans)])>0:
                    is_empty = False

            if is_empty == False:
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

    def plot_median_mat5a_fluxcal(data):
        for wafer in wafer_list:
            obsids = [obsid for obsid in data['MAT5A-pixelraster']]
            f = plt.figure(figsize=(8,6))

            is_empty = True
            for band in [90, 150, 220]:
                median_mat5a = np.array([data['MAT5A-pixelraster'][obsid]['MedianMAT5AFluxCalibration'][wafer][band]
                                         for obsid in obsids
                                         if 'MedianMAT5AFluxCalibration' in data['MAT5A-pixelraster'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianMAT5AFluxCalibration' in data['MAT5A-pixelraster'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_mat5a)],
                         median_mat5a[np.isfinite(median_mat5a)],
                         'o', label='{} GHz'.format(band))

                if len(median_mat5a[np.isfinite(median_mat5a)])>0:
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                plt.ylim([-100, 0])
                plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median MAT5A flux calibration')
            plt.title('MAT5A Flux Calibration ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_mat5a_fluxcal_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_mat5a_intflux(data):
        for wafer in wafer_list:   
            obsids = [obsid for obsid in data['MAT5A-pixelraster']]
            f = plt.figure(figsize=(8,6))
            
            is_empty = True
            for band in [90, 150, 220]:
                median_mat5a = np.array([data['MAT5A-pixelraster'][obsid]['MedianMAT5AIntegralFlux'][wafer][band]
                                         for obsid in obsids 
                                         if 'MedianMAT5AIntegralFlux' in data['MAT5A-pixelraster'][obsid]])

                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianMAT5AIntegralFlux' in data['MAT5A-pixelraster'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(median_mat5a)],
                         median_mat5a[np.isfinite(median_mat5a)],
                         'o', label='{} GHz'.format(band))

                if len(median_mat5a[np.isfinite(median_mat5a)])>0:
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                plt.ylim([2e-7, 7e-7])
                plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('median MAT5A integral flux')
            plt.title('MAT5A Integral Flux ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/median_mat5a_intflux_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_mat5a_sky_transmission(data):
        for wafer in wafer_list:   
            obsids = [obsid for obsid in data['MAT5A']]
            f = plt.figure(figsize=(8,6))
            
            is_empty = True
            for band in [90, 150, 220]:
                mat5a_skytrans = np.array([data['MAT5A'][obsid]['MAT5ASkyTransmission'][wafer][band]
                                           for obsid in obsids
                                           if 'MAT5ASkyTransmission' in data['MAT5A'][obsid]])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                              for obsid in obsids
                              if 'MAT5ASkyTransmission' in data['MAT5A'][obsid]]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)

                plt.plot(datenums[np.isfinite(mat5a_skytrans)],
                         mat5a_skytrans[np.isfinite(mat5a_skytrans)],
                         'o', label='{} GHz'.format(band))

                if len(mat5a_skytrans[np.isfinite(mat5a_skytrans)])>0:
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                plt.ylim([0.85, 1.25])
                plt.legend()
            plt.xlabel('observation time')
            plt.ylabel('MAT5A sky transmission')
            plt.title('MAT5A Sky Transmission ({})'.format(wafer))
            plt.tight_layout()
            plt.savefig('{}/mat5a_sky_transmission_{}.png'.format(outdir, wafer))
            plt.close()

    def plot_median_noise(data, noise_type):
        nex_name = noise_type.split('_')[0]
        labels = {'NEI': 'NEI [pA / sqrt(Hz)]',
                  'NET': 'NET [uK rtsec]',
                  'NEP': 'NEP [aW / sqrt(Hz)]'}
        limits = {'NEI': [0, 100],
                  'NET': [0, 5000],
                  'NEP': [0, 200]}
        units  = {'NEI': core.G3Units.amp*1e-12 / np.sqrt(core.G3Units.Hz),
                  'NET': core.G3Units.microkelvin * np.sqrt(core.G3Units.sec),
                  'NEP': core.G3Units.attowatt / np.sqrt(core.G3Units.Hz)}

        for wafer in wafer_list: 
            obsids = [obsid for obsid in data['noise']]
            f = plt.figure(figsize=(8,6))
            
            is_empty = True
            for band in [90, 150, 220]:
                noise = np.array([data['noise'][obsid][noise_type][wafer][band] / units[nex_name]
                                  for obsid in obsids
                                  if noise_type in data['noise'][obsid].keys()])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds
                              for obsid in obsids
                              if noise_type in data['noise'][obsid].keys()]
                dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
                datenums = mdates.date2num(dts)
                
                plt.plot(datenums[np.isfinite(noise)],
                         noise[np.isfinite(noise)],
                         'o', label='{} GHz'.format(band))

                if len(noise[np.isfinite(noise)])>0:
                    is_empty = False

            if is_empty == False:
                xfmt = mdates.DateFormatter('%m-%d %H:%M')
                plt.gca().xaxis.set_major_formatter(xfmt)
                plt.xticks(rotation=25)
                plt.ylim(limits[nex_name])
                plt.legend()
            plt.xlabel('observation time')
            plt.ylabel(labels[nex_name])
            plt.title('{} ({})'.format(noise_type.replace('_', ' '), wafer))
            plt.tight_layout()
            plt.savefig('{}/median_{}_{}.png'.format(outdir, noise_type, wafer))
            plt.close()


    # only update figures if the underlying data actually changed.
    if was_data_updated or (args.mode == 'update' and args.new_plots):
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
        plot_median_mat5a_fluxcal(data)
        plot_median_mat5a_intflux(data)
        plot_mat5a_sky_transmission(data)
        plot_median_noise(data, 'NEI_0.1Hz_to_0.5Hz')
        plot_median_noise(data, 'NEI_1.0Hz_to_2.0Hz')
        plot_median_noise(data, 'NEI_3.0Hz_to_5.0Hz')
        plot_median_noise(data, 'NEI_10.0Hz_to_15.0Hz')
        plot_median_noise(data, 'NEP_0.1Hz_to_0.5Hz')
        plot_median_noise(data, 'NEP_1.0Hz_to_2.0Hz')
        plot_median_noise(data, 'NEP_3.0Hz_to_5.0Hz')
        plot_median_noise(data, 'NEP_10.0Hz_to_15.0Hz')
        plot_median_noise(data, 'NET_0.1Hz_to_0.5Hz')
        plot_median_noise(data, 'NET_1.0Hz_to_2.0Hz')
        plot_median_noise(data, 'NET_3.0Hz_to_5.0Hz')
        plot_median_noise(data, 'NET_10.0Hz_to_15.0Hz')

        # Create symlink from latest data directory to current
        # The latest data directory is not necessarily reliably identified by
        # the modification timestamp; instead, it is defined by having the
        # latest value of timestamp in its filename. Based on the naming scheme
        # for directories, an alphanumeric sort should also produce
        # chronological ordering, so we'll rely on this assumption.
        symlinkname = '{}/current'.format(args.outdir)
        dirnames = np.sort(glob.glob('{}/*{}'.format(args.outdir, cache_dir_stub)))
        latest_dirname = dirnames[-1]
        os.symlink(latest_dirname, '{}/temp'.format(args.outdir))
        os.rename('{}/temp'.format(args.outdir), '{}/current'.format(args.outdir))
