from statistics import compute_median, compute_nalive
from functools import reduce
import operator
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core
from spt3g.std_processing import obsid_to_g3time
import datetime
import matplotlib.dates as mdates
from plot_util import plot_timeseries
import operator


def alive_bolos_rcw38_fluxcal(frame, boloprops, selector_dict):
    if 'RCW38FluxCalibration' not in frame.keys():
        return None
    return compute_nalive(frame, 'RCW38FluxCalibration', boloprops, selector_dict, 0, operator.lt)


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


def plot_median_rcw38_fluxcal(data, wafers, outdir, xlims=None, ylims=[-100, 0]):
    lines = {}
    
    for wafer in wafers:
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
            plot_timeseries(dts, median_rcw38, band, xlims=xlims, ylims=ylims)

            if len(median_rcw38)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median RCW38 flux calibration')
        plt.title('RCW38 Flux Calibration ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_rcw38_fluxcal_{}.png'.format(outdir, wafer))
        plt.close()


def plot_median_rcw38_intflux(data, wafers, outdir, xlims=None, ylims=[2e-7, 7e-7]):
    lines = {}
    
    for wafer in wafers:   
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
            plot_timeseries(dts, median_rcw38, band, xlims=xlims, ylims=ylims)

            if len(median_rcw38)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median RCW38 integral flux')
        plt.title('RCW38 Integral Flux ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_rcw38_intflux_{}.png'.format(outdir, wafer))
        plt.close()


def plot_rcw38_sky_transmission(data, wafers, outdir, xlims=None, ylims=[0.80, 1.10]):
    lines = {}
    
    for wafer in wafers:   
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
            plot_timeseries(dts, rcw38_skytrans, band, xlims=xlims, ylims=ylims)

            if len(rcw38_skytrans)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('RCW38 sky transmission')
        plt.title('RCW38 Sky Transmission ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/rcw38_sky_transmission_{}.png'.format(outdir, wafer))
        plt.close()


def plot_alive_bolos_rcw38(data, wafers, outdir, xlims=None, ylims=[0, 600], ylims_all=[0, 6000]):
    lines = {}

    for wafer in wafers:
        obsids = [obsid for obsid in data['RCW38-pixelraster']
                  if 'AliveBolosRCW38' in data['RCW38-pixelraster'][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            n_alive_bolos = np.array([data['RCW38-pixelraster'][obsid]['AliveBolosRCW38'][wafer][band]
                                      for obsid in obsids
                                      if 'AliveBolosRCW38' in data['RCW38-pixelraster'][obsid]])
            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                          for obsid in obsids
                          if 'AliveBolosRCW38' in data['RCW38-pixelraster'][obsid]]

            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
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
        plt.xlabel('observation time')
        plt.ylabel('number of alive bolos')
        plt.title('Number of bolos with RCW38 flux calibration < 0 ({})'
                  .format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_rcw38_{}.png'.format(outdir, wafer))
        plt.close()
