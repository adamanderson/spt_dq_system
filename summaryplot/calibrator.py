from statistics import compute_median, compute_nalive
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core


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
