from summaryplot.statistics import compute_median, compute_nalive
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, dfmux
from spt3g.std_processing import obsid_to_g3time
import datetime
import matplotlib.dates as mdates
from summaryplot.plot_util import plot_timeseries
import operator
import pickle

opt_eff_dict = None


def median_pelec_per_airmass(frame, boloprops, selector_dict):
    if 'ElnodSlopesPower' not in frame:
        return None
    return compute_median(frame, 'ElnodSlopesPower', boloprops, selector_dict)


def median_elnod_opacity(frame, boloprops, selector_dict):
    if 'ElnodSlopesPower' not in frame:
        return None

    pelec_per_airmass_dict = frame['ElnodSlopesPower']

    global opt_eff_dict
    if opt_eff_dict is None:
        with open('summaryplot/optical_efficiency.pickle', 'rb') as fobj:
            opt_eff_dict = pickle.load(fobj)

    opacities = core.G3MapDouble()
    obsid = frame['ObservationID']
    month = str(obsid_to_g3time(obsid)).split('-')[1].lower()
    tsky_dict = {'jan':247, 'feb':235, 'mar':221, 'apr':217,
                 'may':217, 'jun':217, 'jul':215, 'aug':215,
                 'sep':215, 'oct':223, 'nov':236, 'dec':247}
    deltanu_dict = {90.0: 28.2e9, 150.0: 42.6e9, 220.0: 52.1e9}
    kb = 1.38064852e-23
    tsky = tsky_dict[month]
    for bolo, ppa in pelec_per_airmass_dict.items():
        band = boloprops[bolo].band / core.G3Units.GHz
        if band in [90.0, 150.0, 220.0]:
            deltanu = deltanu_dict[band]
            try:
                eta = opt_eff_dict[bolo]
                tau = (-1.0*ppa/core.G3Units.W) / (eta*deltanu*kb*tsky)
                opacities[bolo] = tau
            except KeyError:
                pass

    frame['ElnodOpacity'] = opacities
    return compute_median(frame, 'ElnodOpacity', boloprops, selector_dict)


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


def plot_median_elnod_response(data, wafers, outdir, xlims=None, ylims=[-6.0, 0.0]):
    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))
        
        is_empty = True
        for band in [90, 150, 220]:
            median_ppa = np.array([data['elnod'][obsid]['MedianPelecPerAirmass'][wafer][band] for obsid in data['elnod']]) / core.G3Units.pW
            
            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds for obsid in obsids]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_ppa, band, xlims=xlims, ylims=ylims)
            
            if len(median_ppa)>0:
                is_empty = False
        
        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
        
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median elnod slope [pW/airmass]')
        plt.title('Elnod slope ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_response_{}.png'.format(outdir, wafer))
        plt.close()


def plot_median_elnod_sn(data, wafers, outdir, xlims=None, ylims=[0, 4000]):
    lines = {}

    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            median_elnodSN = np.array([data['elnod'][obsid]['MedianElnodSNSlopes'][wafer][band] for obsid in data['elnod']])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                          for obsid in obsids]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_elnodSN, band, xlims=xlims, ylims=ylims)

            if len(median_elnodSN)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median elnod slope S/N')
        plt.title('Elnod slope S/N ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_sn_{}.png'.format(outdir, wafer))
        plt.close()


def plot_median_elnod_opacity(data, wafers, outdir, xlims=None, ylims=[0.0, 0.20]):
    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))
        
        is_empty = True
        for band in [90, 150, 220]:
            median_ppa = np.array([data['elnod'][obsid]['MedianElnodOpacity'][wafer][band] for obsid in data['elnod']])
            
            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds for obsid in obsids]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_ppa, band, xlims=xlims, ylims=ylims)
            
            if len(median_ppa)>0:
                is_empty = False
        
        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
        
        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median opacity [airmass^-1]')
        plt.title('Sky opacity (optical depth) derived from elnod slope ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_opacity_{}.png'.format(outdir, wafer))
        plt.close()


def plot_median_elnod_iq_phase(data, wafers, outdir, xlims=None, ylims=[-20, 20]):
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
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, median_elnod_iq, band, xlims=xlims, ylims=ylims)

            if len(median_elnod_iq)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('median elnod IQ phase [deg]')
        plt.title('Elnod IQ phase angle ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_elnod_iq_phase_{}.png'.format(outdir, wafer))
        plt.close()


def plot_alive_bolos_elnod(data, wafers, outdir, xlims=None, ylims=[0, 600], ylims_all=[0, 6000]):
    lines = {}

    for wafer in wafers:
        obsids = [obsid for obsid in data['elnod']]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            alive_bolos_elnod = np.array([data['elnod'][obsid]['AliveBolosElnod'][wafer][band] for obsid in data['elnod']])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.seconds \
                          for obsid in obsids]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            if wafer == 'all':
                plot_timeseries(dts, alive_bolos_elnod, band, xlims=xlims, ylims=ylims_all)
            else:
                plot_timeseries(dts, alive_bolos_elnod, band, xlims=xlims, ylims=ylims)

            if len(alive_bolos_elnod)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel('number of alive bolos')
        plt.title('Number of bolos with elnod S/N>20 ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/alive_bolos_elnod_{}.png'.format(outdir, wafer))
        plt.close()
