from statistics import compute_median
import numpy as np
import matplotlib.pyplot as plt
from spt3g import core
from spt3g.std_processing import obsid_to_g3time
import matplotlib.dates as mdates
import datetime

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


def median_net_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NET_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_3.0Hz_to_5.0Hz', boloprops, selector_dict)    


def median_net_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NET_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_10.0Hz_to_15.0Hz', boloprops, selector_dict)    


def plot_median_noise(data, noise_type, wafers, outdir):
    # min/max for plotting purposes
    lines = {}
    
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

    for wafer in wafers: 
        l_nan = None

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

            lines[band], = plt.plot(datenums[np.isfinite(noise)],
                                    noise[np.isfinite(noise)],
                                    'o', label='{} GHz'.format(band))
            
            # plot light dashed lines for NaNs
            nan_dates = datenums[~np.isfinite(noise)]
            for date in nan_dates:
                l_nan, = plt.plot([date, date],
                                  limits[nex_name],
                                  'r--', linewidth=0.5)

            if len(noise)>0:
                is_empty = False

        if is_empty == False:
            xfmt = mdates.DateFormatter('%m-%d %H:%M')
            plt.gca().xaxis.set_major_formatter(xfmt)
            plt.xticks(rotation=25)
            plt.ylim(limits[nex_name])
            if l_nan != None:
                plt.legend((lines[90], lines[150], lines[220], l_nan),
                           ('90 GHz', '150 GHz', '220 GHz', 'NaNs'))
            else:
                plt.legend((lines[90], lines[150], lines[220]),
                           ('90 GHz', '150 GHz', '220 GHz'))

        plt.grid()
        plt.xlabel('observation time')
        plt.ylabel(labels[nex_name])
        plt.title('{} ({})'.format(noise_type.replace('_', ' '), wafer))
        plt.tight_layout()
        plt.savefig('{}/median_{}_{}.png'.format(outdir, noise_type, wafer))
        plt.close()
