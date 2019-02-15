from statistics import compute_median, compute_nalive
import numpy as np
import matplotlib.pyplot as plt
from spt3g.std_processing import obsid_to_g3time
from spt3g import core
import datetime
import matplotlib.dates as mdates
import numpy as np


def median_cal_sn_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponseSN', boloprops, selector_dict)

def median_cal_response_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponse' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponse', boloprops, selector_dict)

def alive_bolos_cal_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_nalive(frame, 'CalibratorResponseSN', boloprops, selector_dict, 20)

def mean_cal_elevation(rawpath, boloprops):
    class ElExtractor(object):
        def __init__(self):
            self.mean_el = np.nan
        def __call__(self, frame):
            if frame.type == core.G3FrameType.Scan and \
               "OnlineBoresightEl" in frame.keys():
                self.mean_el = np.mean(frame["OnlineBoresightEl"]) / core.G3Units.degree

    el_extractor = ElExtractor()
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=rawpath)
    pipe.Add(el_extractor)
    pipe.Run()

    return el_extractor.mean_el


def plot_median_cal_sn_4Hz(data, wafers, outdir, el):
    # min/max for plotting purposes
    ymin = 0
    ymax = 400
    lines = {}
    
    for wafer in wafers:
        l_nan = None

        obsids = [obsid for obsid in data['calibrator']
                  if 'MedianCalSN_4Hz' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                median_calSN = np.array([data['calibrator'][obsid]['MedianCalSN_4Hz'][wafer][band]
                                         for obsid in data['calibrator']
                                         if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                                         data['calibrator'][obsid]['elevation']>40 and \
                                         data['calibrator'][obsid]['elevation']<56])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>40 and \
                              data['calibrator'][obsid]['elevation']<56]
            elif el == 'high':
                median_calSN = np.array([data['calibrator'][obsid]['MedianCalSN_4Hz'][wafer][band]
                                         for obsid in data['calibrator']
                                         if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                                         data['calibrator'][obsid]['elevation']>56 and \
                                         data['calibrator'][obsid]['elevation']<72])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>56 and \
                              data['calibrator'][obsid]['elevation']<72]

            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
            datenums = mdates.date2num(dts)

            lines[band], = plt.plot(datenums[np.isfinite(median_calSN)],
                                    median_calSN[np.isfinite(median_calSN)],
                                    'o', label='{} GHz'.format(band))

            # plot light dashed lines for NaNs
            nan_dates = datenums[~np.isfinite(median_calSN)]
            for date in nan_dates:
                l_nan, = plt.plot([date, date], [ymin, ymax], 'r--', linewidth=0.5)

            if len(median_calSN)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
            plt.ylim([ymin, ymax])
            if l_nan != None:
                plt.legend((lines[90], lines[150], lines[220], l_nan),
                           ('90 GHz', '150 GHz', '220 GHz', 'NaNs'))
            else:
                plt.legend((lines[90], lines[150], lines[220]),
                           ('90 GHz', '150 GHz', '220 GHz'))

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median calibrator S/N')
        if el == 'low':
            el_title = '40 < el < 56'
        elif el == 'high':
            el_title = '56 < el < 72'
        plt.title('4.0 Hz calibrator S/N ({})\n{}'
                  .format(wafer, el_title))
        plt.tight_layout()
        plt.savefig('{}/median_cal_sn_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()


def plot_median_cal_response_4Hz(data, wafers, outdir, el):
    # min/max for plotting purposes
    ymin = 0
    ymax = 5
    lines = {}
    
    for wafer in wafers:
        l_nan = None

        obsids = [obsid for obsid in data['calibrator']
                  if 'MedianCalResponse_4Hz' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                median_cal = np.array([data['calibrator'][obsid]['MedianCalResponse_4Hz'][wafer][band]
                                       for obsid in data['calibrator']
                                       if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                                       data['calibrator'][obsid]['elevation']>40 and \
                                       data['calibrator'][obsid]['elevation']<56])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>40 and \
                              data['calibrator'][obsid]['elevation']<56]
            elif el == 'high':
                median_cal = np.array([data['calibrator'][obsid]['MedianCalResponse_4Hz'][wafer][band]
                                       for obsid in data['calibrator']
                                       if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                                       data['calibrator'][obsid]['elevation']>56 and \
                                       data['calibrator'][obsid]['elevation']<72])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>56 and \
                              data['calibrator'][obsid]['elevation']<72]

            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
            datenums = mdates.date2num(dts)

            lines[band], = plt.plot(datenums[np.isfinite(median_cal)],
                                    1e15*median_cal[np.isfinite(median_cal)],
                                    'o', label='{} GHz'.format(band))

            # plot light dashed lines for NaNs
            nan_dates = datenums[~np.isfinite(median_cal)]
            if len(nan_dates) > 0:
                print(wafer)
                print(nan_dates)
            for date in nan_dates:
                l_nan, = plt.plot([date, date], [ymin, ymax], 'r--', linewidth=0.5)

            if len(median_cal)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
            plt.ylim([ymin, ymax])
            if l_nan != None:
                plt.legend((lines[90], lines[150], lines[220], l_nan),
                           ('90 GHz', '150 GHz', '220 GHz', 'NaNs'))
            else:
                plt.legend((lines[90], lines[150], lines[220]),
                           ('90 GHz', '150 GHz', '220 GHz'))
        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median calibrator response [fW]')
        if el == 'low':
            el_title = '40 < el < 56'
        elif el == 'high':
            el_title = '56 < el < 72'
        plt.title('4.0 Hz calibrator response ({})\n{}'.
                  format(wafer, el_title))
        plt.tight_layout()
        plt.savefig('{}/median_cal_response_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()


def plot_alive_bolos_cal_4Hz(data, wafers, outdir, el):
    # min/max for plotting purposes
    ymin = 0
    ymax = 600
    ymax_all = 6000
    lines = {}
    
    for wafer in wafers:
        l_nan = None
        
        obsids = [obsid for obsid in data['calibrator']
                  if 'AliveBolosCal_4Hz' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                n_alive_bolos = np.array([data['calibrator'][obsid]['AliveBolosCal_4Hz'][wafer][band]
                                          for obsid in data['calibrator']
                                          if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                                          data['calibrator'][obsid]['elevation']>40 and \
                                          data['calibrator'][obsid]['elevation']<56])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>40 and \
                              data['calibrator'][obsid]['elevation']<56]
            elif el == 'high':
                n_alive_bolos = np.array([data['calibrator'][obsid]['AliveBolosCal_4Hz'][wafer][band]
                                          for obsid in data['calibrator']
                                          if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                                          data['calibrator'][obsid]['elevation']>56 and \
                                          data['calibrator'][obsid]['elevation']<72])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>56 and \
                              data['calibrator'][obsid]['elevation']<72]
            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
            datenums = mdates.date2num(dts)

            lines[band], = plt.plot(datenums[np.isfinite(n_alive_bolos)],
                                    n_alive_bolos[np.isfinite(n_alive_bolos)],
                                    'o', label='{} GHz'.format(band))

            # plot light dashed lines for NaNs
            nan_dates = datenums[~np.isfinite(n_alive_bolos)]
            if wafer == 'all':
                for date in nan_dates:
                    l_nan, = plt.plot([date, date], [ymin, ymax], 'r--', linewidth=0.5)
            else:
                for date in nan_dates:
                    l_nan, = plt.plot([date, date], [ymin, ymax_all], 'r--', linewidth=0.5)

            if len(n_alive_bolos)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
            if wafer == 'all':
                plt.ylim([0, ymax_all])
            else:
                plt.ylim([0, ymax])
            if l_nan != None:
                plt.legend((lines[90], lines[150], lines[220], l_nan),
                           ('90 GHz', '150 GHz', '220 GHz', 'NaNs'))
            else:
                plt.legend((lines[90], lines[150], lines[220]),
                           ('90 GHz', '150 GHz', '220 GHz'))
        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('number of alive bolos')
        if el == 'low':
            el_title = '40 < el < 56'
        elif el == 'high':
            el_title = '56 < el < 72'
        plt.title('Number of bolos with calibrator S/N > 20 at 4.0 Hz ({})\n{}'
                  .format(wafer, el_title))
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_cal_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()
