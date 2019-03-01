from statistics import compute_median, compute_nalive
from spt3g import core
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core
from spt3g.std_processing import obsid_to_g3time
import datetime
import matplotlib.dates as mdates
from plot_util import plot_timeseries
import operator

def median_elnod_sn_slope(frame, boloprops, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_median(frame, 'ElnodSNSlopes', boloprops, selector_dict)


def median_elnod_iq_phase_angle(frame, boloprops, selector_dict):
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

    return compute_median(newframe, 'ElnodPhaseAngle', boloprops, selector_dict)



def alive_bolos_elnod(frame, boloprops, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_nalive(frame, 'ElnodSNSlopes', boloprops, selector_dict, 20, operator.gt)


def plot_median_elnod_sn(data, wafers, outdir):
    # min/max for plotting purposes
    ymin = 0
    ymax = 2000
    lines = {}

    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            median_elnodSN = np.array([data['elnod'][obsid]['MedianElnodSNSlopes'][wafer][band] for obsid in data['elnod']])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                          for obsid in obsids]
            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_elnodSN, [ymin, ymax], band)

            if len(median_elnodSN)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median elnod S/N')
        plt.title('Elnod S/N ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_sn_{}.png'.format(outdir, wafer))
        plt.close()


def plot_median_elnod_iq_phase(data, wafers, outdir):
    # min/max for plotting purposes
    ymin = -20
    ymax = 20
    lines = {}
    
    for wafer in wafers:
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
            plot_timeseries(dts, median_elnod_iq, [ymin, ymax], band)

            if len(median_elnod_iq)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('median elnod IQ phase [deg]')
        plt.title('Elnod IQ phase angle ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_iq_phase_{}.png'.format(outdir, wafer))
        plt.close()


def plot_alive_bolos_elnod(data, wafers, outdir):
    # min/max for plotting purposes
    ymin = 0
    ymax = 600
    ymax_all = 6000
    lines = {}

    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            alive_bolos_elnod = np.array([data['elnod'][obsid]['AliveBolosElnod'][wafer][band] for obsid in data['elnod']])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                          for obsid in obsids]
            dts = np.array([datetime.datetime.fromtimestamp(ts) for ts in timestamps])
            if wafer == 'all':
                plot_timeseries(dts, alive_bolos_elnod, [ymin, ymax_all], band)
            else:
                plot_timeseries(dts, alive_bolos_elnod, [ymin, ymax], band)

            if len(alive_bolos_elnod)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel('number of alive bolos')
        plt.title('Number of bolos with elnod S/N>20 ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_elnod_{}.png'.format(outdir, wafer))
        plt.close()
