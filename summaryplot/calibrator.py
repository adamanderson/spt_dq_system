from summaryplot.statistics import compute_median, compute_nalive
import numpy as np
import matplotlib.pyplot as plt
from spt3g.std_processing import obsid_to_g3time
from spt3g import core
import datetime
import matplotlib.dates as mdates
import numpy as np
from summaryplot.plot_util import plot_timeseries
import operator

low_el_range  = [20, 56]
high_el_range = [56, 72]

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
    return compute_nalive(frame, 'CalibratorResponseSN', boloprops, selector_dict, 20, operator.gt)

def mean_cal_elevation(frame, boloprops, selector_dict):
    if 'CalibratorResponseEl' not in frame.keys():
        return None
    return frame['CalibratorResponseEl'] / core.G3Units.degree


def plot_median_cal_sn_4Hz(data, wafers, outdir, el, xlims=None, ylims=[0, 400]):
    lines = {}
    colors = {90: 'C0', 150: 'C1', 220: 'C2'}
    
    for wafer in wafers:
        obsids = [obsid for obsid in data['calibrator']
                  if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                  'elevation' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                median_calSN = np.array([data['calibrator'][obsid]['MedianCalSN_4Hz'][wafer][band]
                                         for obsid in obsids
                                         if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                                         data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                                         data['calibrator'][obsid]['elevation']<low_el_range[1]])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<low_el_range[1]]
            elif el == 'high':
                median_calSN = np.array([data['calibrator'][obsid]['MedianCalSN_4Hz'][wafer][band]
                                         for obsid in obsids
                                         if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                                         data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                                         data['calibrator'][obsid]['elevation']<high_el_range[1]])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalSN_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<high_el_range[1]]

            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_calSN, band, xlims=xlims, ylims=ylims)

            if len(median_calSN)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median calibrator S/N')
        if el == 'low':
            el_title = '{} < el < {}'.format(low_el_range[0], low_el_range[1])
        elif el == 'high':
            el_title = '{} < el < {}'.format(high_el_range[0], high_el_range[1])
        plt.title('4.0 Hz calibrator S/N ({})\n{}'
                  .format(wafer, el_title))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_cal_sn_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()


def plot_median_cal_response_4Hz(data, wafers, outdir, el, xlims=None, ylims=[0, 5]):
    lines = {}
    
    colors = {90: 'C0', 150: 'C1', 220: 'C2'}

    for wafer in wafers:
        obsids = [obsid for obsid in data['calibrator']
                  if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                  'elevation' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                median_cal = np.array([data['calibrator'][obsid]['MedianCalResponse_4Hz'][wafer][band] / \
                                       (core.G3Units.watt*1e-15)  
                                       for obsid in obsids
                                       if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                                       data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                                       data['calibrator'][obsid]['elevation']<low_el_range[1]])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<low_el_range[1]]
            elif el == 'high':
                median_cal = np.array([data['calibrator'][obsid]['MedianCalResponse_4Hz'][wafer][band] / \
                                       (core.G3Units.watt*1e-15)
                                       for obsid in obsids
                                       if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                                       data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                                       data['calibrator'][obsid]['elevation']<high_el_range[1]])                
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'MedianCalResponse_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<high_el_range[1]]

            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_cal, band, xlims=xlims, ylims=ylims)

            if len(median_cal)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median calibrator response [fW]')
        if el == 'low':
            el_title = '{} < el < {}'.format(low_el_range[0], low_el_range[1])
        elif el == 'high':
            el_title = '{} < el < {}'.format(high_el_range[0], high_el_range[1])
        plt.title('4.0 Hz calibrator response ({})\n{}'.
                  format(wafer, el_title))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_cal_response_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()


def plot_alive_bolos_cal_4Hz(data, wafers, outdir, el, xlims=None, ylims=[0, 600], ylims_all=[0, 6000]):
    lines = {}
    
    for wafer in wafers:
        obsids = [obsid for obsid in data['calibrator']
                  if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                  'elevation' in data['calibrator'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            if el == 'low':
                n_alive_bolos = np.array([data['calibrator'][obsid]['AliveBolosCal_4Hz'][wafer][band]
                                          for obsid in obsids
                                          if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                                          data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                                          data['calibrator'][obsid]['elevation']<low_el_range[1]])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>low_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<low_el_range[1]]
            elif el == 'high':
                n_alive_bolos = np.array([data['calibrator'][obsid]['AliveBolosCal_4Hz'][wafer][band]
                                          for obsid in obsids
                                          if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                                          data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                                          data['calibrator'][obsid]['elevation']<high_el_range[1]])
                timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                              for obsid in obsids
                              if 'AliveBolosCal_4Hz' in data['calibrator'][obsid] and \
                              data['calibrator'][obsid]['elevation']>high_el_range[0] and \
                              data['calibrator'][obsid]['elevation']<high_el_range[1]]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            if wafer == 'all':
                plot_timeseries(dts, n_alive_bolos, band, xlims=xlims, ylims=ylims_all)
            else:
                plot_timeseries(dts, n_alive_bolos, band, xlims=xlims, ylims=ylims)

            if len(n_alive_bolos)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
            
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('number of alive bolos')
        if el == 'low':
            el_title = '{} < el < {}'.format(low_el_range[0], low_el_range[1])
        elif el == 'high':
            el_title = '{} < el < {}'.format(high_el_range[0], high_el_range[1])
        plt.title('Number of bolos with calibrator S/N > 20 at 4.0 Hz ({})\n{}'
                  .format(wafer, el_title))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_cal_4Hz_{}el_{}.png'.format(outdir, el, wafer))
        plt.close()
