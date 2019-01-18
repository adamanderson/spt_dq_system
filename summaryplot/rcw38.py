from statistics import compute_median
from functools import reduce
import operator
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core

def median_rcw38_fluxcal(frame, boloprops, selector_dict):
    if 'RCW38FluxCalibration' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38FluxCalibration', boloprops, selector_dict)


def median_rcw38_intflux(frame, boloprops, selector_dict):
    if 'RCW38IntegralFlux' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38IntegralFlux', boloprops, selector_dict)


def rcw38_sky_transmission(frame, boloprops, selector_dict):
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
        else:
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = np.nan

    return data_on_selection


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

