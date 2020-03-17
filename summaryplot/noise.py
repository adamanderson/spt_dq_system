from summaryplot.statistics import compute_median
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core, calibration, dfmux, todfilter
from spt3g.std_processing import obsid_to_g3time
import gc
import matplotlib.dates as mdates
import datetime
from summaryplot.plot_util import plot_timeseries


def median_nei_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NEI_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_0.1Hz_to_0.5Hz', boloprops, selector_dict)    


def median_nei_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NEI_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_1.0Hz_to_2.0Hz', boloprops, selector_dict)    


def median_nei_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NEI_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_3.0Hz_to_5.0Hz', boloprops, selector_dict)    


def median_nei_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NEI_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_10.0Hz_to_15.0Hz', boloprops, selector_dict)    


def median_nep_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NEP_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_0.1Hz_to_0.5Hz', boloprops, selector_dict)    


def median_nep_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NEP_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_1.0Hz_to_2.0Hz', boloprops, selector_dict)    


def median_nep_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NEP_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_3.0Hz_to_5.0Hz', boloprops, selector_dict)    


def median_nep_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NEP_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_10.0Hz_to_15.0Hz', boloprops, selector_dict)    


def median_net_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NET_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_0.1Hz_to_0.5Hz', boloprops, selector_dict)    


def median_net_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NET_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_1.0Hz_to_2.0Hz', boloprops, selector_dict)



def number_of_lines_in_median_psds(rawpath, cal_fname, fname, boloprops, selector_dict):

    pipeline = core.G3Pipeline()
    pipeline.Add(core.G3Reader,
                 filename=[fname, cal_fname, rawpath])
    pipeline.Add(core.Delete,
                 keys=['RawTimestreams_Q'])
    pipeline.Add(dfmux.ConvertTimestreamUnits,
                 Input='RawTimestreams_I',
                 Output='Timestreams_I_inPower',
                 Units=core.G3TimestreamUnits.Power)
    pipeline.Add(core.Delete,
                 keys=['RawTimestreams_I'])
    class ExtractVariousMaps(object):
        def __init__(self):
            self.tm  = None
            self.bpm = None
            self.wm  = None
        def __call__(self, frame):
            if 'Timestreams_I_inPower' in frame.keys():
                tm = frame['Timestreams_I_inPower']
                self.tm = core.G3TimestreamMap()
                for k, ts in tm.items():
                    self.tm[k] = ts - np.mean(ts)
                return
            if 'NominalBolometerProperties' in frame.keys():
                self.bpm = frame['NominalBolometerProperties']
                return
            if 'WiringMap' in frame.keys():
                self.wm = frame['WiringMap']
                return
    mapsxtractor = ExtractVariousMaps()
    pipeline.Add(mapsxtractor)
    pipeline.Run()

    lines = todfilter.notchfilter.find_lines_from_PSDs(
                mapsxtractor.tm, [mapsxtractor.tm.n_samples], mapsxtractor.bpm, mapsxtractor.wm,
                from_by_group_median_PSD_only=True, by_band=True, by_wafer=True,
                factor_for_thresholds=2.00,
                factor_for_baseline=1.20,
                power_law_fit_starting_frequency=0.06,
                power_law_fit_ending_frequency=0.60,
                white_noise_estimate_starting_frequency=4.5,
                white_noise_estimate_ending_frequency=7.5,
                low_frequency_cut=0.99,
                maximum_number_of_lines=50)
    del mapsxtractor
    gc.collect()

    returndict = {}
    for groupdest in selector_dict.keys():
        for groupsource in lines['LineLocations'].keys():
            band  = float(groupsource.split('_')[0])
            wafer = groupsource.split('_')[1]
            if (groupdest[0] == wafer) and (float(groupdest[1]) == band):
                if wafer not in returndict.keys():
                    returndict[wafer] = {}
                returndict[wafer][int(band)] = len(lines['LineLocations'][groupsource])
    for groupdest in selector_dict.keys():
        wafer = groupdest[0]
        if wafer not in returndict.keys():
            returndict[wafer] = {}
            returndict[wafer] = {90: 0, 150: 0, 220: 0}
    totlines = {}
    for oneband in [90, 150, 220]:
        totlines[oneband] = np.sum([lines3bands[oneband] for wafer, lines3bands in returndict.items()])
    returndict['all'] = totlines
    return [returndict]



def median_net_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NET_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_3.0Hz_to_5.0Hz', boloprops, selector_dict)    


def median_net_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NET_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_10.0Hz_to_15.0Hz', boloprops, selector_dict)    


def plot_median_noise(data, noise_type, wafers, outdir, xlims=None,
                      ylims={'NEI': [0, 100], 'NET': [0, 5000], 'NEP': [0, 200]}):
    # min/max for plotting purposes
    lines = {}
    
    nex_name = noise_type.split('_')[0]
    labels = {'NEI': 'NEI [pA / sqrt(Hz)]',
              'NET': 'NET [uK rtsec]',
              'NEP': 'NEP [aW / sqrt(Hz)]'}
    units  = {'NEI': core.G3Units.amp*1e-12 / np.sqrt(core.G3Units.Hz),
              'NET': core.G3Units.microkelvin * np.sqrt(core.G3Units.sec),
              'NEP': core.G3Units.attowatt / np.sqrt(core.G3Units.Hz)}

    for wafer in wafers: 
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
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            plot_timeseries(dts, noise, band, xlims=xlims, ylims=ylims[nex_name])

            if len(noise)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        plt.ylabel(labels[nex_name])
        plt.title('{} ({})'.format(noise_type.replace('_', ' '), wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/median_{}_{}.png'.format(outdir, noise_type, wafer))
        plt.close()



def plot_number_of_lines(data, wafers, outdir,
                         xlims=None, ylims=[-0.3, 10], ylims_all=[-0.9, 30]):
    for wafer in wafers:
        obsids = [obsid for obsid in data['noise']]
        f = plt.figure(figsize=(8,6))

        is_empty = True
        for band in [90, 150, 220]:
            number_lines = np.array([data['noise'][obsid]['NumberOfLinesInMedianPSDs'][wafer][band]
                                     for obsid in data['noise']])

            timestamps = [obsid_to_g3time(int(obsid)).time / core.G3Units.sec
                          for obsid in obsids]
            dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])
            if wafer == 'all':
                plot_timeseries(dts, number_lines, band, xlims=xlims, ylims=ylims_all)
            else:
                plot_timeseries(dts, number_lines, band, xlims=xlims, ylims=ylims)

            if len(number_lines)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)

        plt.grid()
        plt.xlabel('observation time (UTC)')
        if wafer == 'all':
            plt.ylabel('number of lines found (sum from all wafers)')
        else:
            plt.ylabel('number of lines found')
        plt.title('Number of lines found in median PSD ({})'.format(wafer))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/number_of_lines_found_{}.png'.format(outdir, wafer))
        plt.close()
