from summaryplot.statistics import compute_median, compute_nalive
from functools import reduce
import operator
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core
from spt3g.std_processing import obsid_to_g3time
import datetime
import matplotlib.dates as mdates
from summaryplot.plot_util import plot_timeseries


def alive_bolos_htwo_fluxcal(frame, boloprops, selector_dict):
    for k in frame.keys():
        if k.endswith('FluxCalibration'):
            cal_key = k
            break
    return compute_nalive(frame, cal_key, boloprops, selector_dict, 0, operator.lt)


def median_htwo_fluxcal(frame, boloprops, selector_dict):
    for k in frame.keys():
        if k.endswith('FluxCalibration'):
            cal_key = k
            break
    return compute_median(frame, cal_key, boloprops, selector_dict)


def median_htwo_intflux(frame, boloprops, selector_dict):
    for k in frame.keys():
        if k.endswith('IntegralFlux'):
            cal_key = k
            break
    return compute_median(frame, cal_key, boloprops, selector_dict)


def htwo_sky_transmission(frame, boloprops, selector_dict):
    for k in frame.keys():
        if k.endswith('SkyTransmission'):
            cal_key = k
            break

    data_on_selection = {}

    for select_values, f_select_list in selector_dict.items():
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # get data that satisfies the selection and compute median
        if str(select_values[-1]) in frame[cal_key].keys():
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = frame[cal_key][str(select_values[-1])]
        else:
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = np.nan

    return data_on_selection


def plot_median_htwo_fluxcal(src, data, wafers, outdir, xlims=None, ylims=[-100, 0]):
    lines = {}
    cal_k = 'Median' + src + 'FluxCalibration'
    
    for wafer in wafers:
        src += '-pixelraster'
        obsids = [obsid for obsid in data[src]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            medians = np.array([data[src][obsid][cal_k][wafer][band]
                                for obsid in obsids
                                if cal_k in data[src][obsid]])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.sec
                          for obsid in obsids
                          if cal_k in data[src][obsid]]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts)
                            for ts in timestamps])
            plot_timeseries(dts, medians, band, xlims=xlims, ylims=ylims)

            if len(medians)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        src = src.replace('-pixelraster', '')
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median {} flux calibration'.format(src))
        plt.title('{} Flux Calibration ({})'.format(src, wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_{}_fluxcal_{}.png'.format(
                    outdir, src.lower(), wafer))
        plt.close()

def plot_median_htwo_intflux(src, data, wafers, outdir, xlims=None, ylims=[1, 8]):
    lines = {}
    cal_k = 'Median' + src + 'IntegralFlux'
    
    for wafer in wafers:   
        src += '-pixelraster'
        obsids = [obsid for obsid in data[src]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            medians = np.array([data[src][obsid][cal_k][wafer][band]
                                for obsid in obsids
                                if  cal_k in data[src][obsid]])
            medians /= core.G3Units.arcmin**2

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.sec
                          for obsid in obsids
                          if  cal_k in data[src][obsid]]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts)
                            for ts in timestamps])
            if 'RCW38' in src:
                ylims = [3.25, 6.75]
            elif 'MAT5A' in src:
                ylims = [1.75, 5.25]
            plot_timeseries(dts, medians, band, xlims=xlims, ylims=ylims)

            if len(medians)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        src = src.replace('-pixelraster', '')
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median {} integral flux [arcmin^2]'.format(src))
        plt.title('{} Integral Flux ({})'.format(src, wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_{}_intflux_{}.png'.format(
                    outdir, src.lower(), wafer))
        plt.close()


def plot_htwo_skytrans(src, data, wafers, outdir, xlims=None, ylims=[0.80, 1.20]):
    lines = {}
    cal_k = src + 'SkyTransmission'
    
    for wafer in wafers:   
        obsids = [obsid for obsid in data[src]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            transs = np.array([data[src][obsid][cal_k][wafer][band]
                               for obsid in obsids
                               if  cal_k in data[src][obsid]])
            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.sec
                          for obsid in obsids
                          if  cal_k in data[src][obsid]]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts)
                  for ts in timestamps])
            plot_timeseries(dts, transs, band, xlims=xlims, ylims=ylims)

            if len(transs)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('{} sky transmission'.format(src))
        plt.title('{} Sky Transmission ({})'.format(src, wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/{}_sky_transmission_{}.png'.format(
                     outdir, src.lower(), wafer))
        plt.close()


def plot_alive_bolos_htwo(src, data, wafers, outdir, xlims=None, ylims=[0, 600], ylims_all=[0, 6000]):
    lines = {}
    cal_k = 'AliveBolos' + src

    for wafer in wafers:
        src += '-pixelraster'
        obsids = [obsid for obsid in data[src]
                  if cal_k in data[src][obsid]]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            n_alive_bolos = np.array([data[src][obsid][cal_k][wafer][band]
                                      for obsid in obsids
                                      if  cal_k in data[src][obsid]])
            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.sec
                          for obsid in obsids
                          if  cal_k in data[src][obsid]]

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

        src = src.replace('-pixelraster', '')
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('number of alive bolos')
        plt.title('Number of bolos with {} flux calibration < 0 ({})'
                  .format(src, wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_{}_{}.png'.format(
                    outdir, src.lower(), wafer))
        plt.close()
