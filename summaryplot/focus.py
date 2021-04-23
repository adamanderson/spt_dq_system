from spt3g import core, mapmaker, util, pointing, std_processing
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime
import os
from summaryplot.plot_util import plot_timeseries
from spt3g.std_processing import obsid_to_g3time
from spt3g.std_processing import time_to_obsid


focus_sources = [
    'PMNJ0210-5101',
    'PMNJ0522-3628',
    'ra0hdec-52.25',
    'ra1h40dec-33.25',
    'ra5hdec-33.25',
    'ra12h30dec-33.25',
]


def extract_benchpos(frames, boloprops=None, selector_dict=None):

    for frame in frames:
        if frame.type == core.G3FrameType.Observation:
            bp = frame['BenchCommandedPosition']
            bp_zeros = frame['BenchZeros']
            for k in bp.keys():
                bp[k] -= bp_zeros[k]
            optical_axis_pos = pointing.focus.bench2optical(bp['x4'], bp['y1'], bp['z6']).A1[2]
            break

    return optical_axis_pos



def benchpos_min_fwhm_ellip(frames, boloprops=None, selector_dict=None):

    out = {}
    out['FWHM'] = {}
    out['ellipticity'] = {}
    out['bench position'] = extract_benchpos(frames)

    for frame in frames:
        if frame.type == core.G3FrameType.Map:
            if not ('Id' in frame and 'GHz' in frame['Id']):
                continue

            band = frame['Id'][-6:].replace('-','')
            tmap = frame['T'] / frame['Wunpol'].TT
            res = tmap.res
            
            old_shape = tmap.shape
            yc = old_shape[0] // 2 
            xc = old_shape[1] // 2
            width = int(np.ceil(45.*core.G3Units.arcmin / res / 2.))
            tmap  = np.asarray(tmap)[yc-width:yc+width, xc-width:xc+width]

            new_shape = tmap.shape
            width = int(np.ceil(10.*core.G3Units.arcmin / res / 2.))
            cy, cx = np.unravel_index(np.nanargmax(tmap), new_shape)

            tmap = np.asarray(tmap)[cy-width : cy+width, cx-width:cx+width]
            try:
                params = util.fitting.fit_gaussian2d(tmap)
                (amp, x0, y0, sigma_x, sigma_y, theta, height) = params
                sigma_to_fwhm = 2.*np.sqrt(2.*np.log(2))
                fwhm_px = sigma_to_fwhm*np.sqrt(0.5*sigma_x**2 + 0.5*sigma_y**2)
                fwhm_am = fwhm_px * res / core.G3Units.arcmin
                ellipticity = np.sqrt(np.abs(sigma_x**2-sigma_y**2) / np.max((sigma_x**2, sigma_y**2)))
            except:
                fwhm_am = np.nan
                ellipticity = np.nan
            out['FWHM'][band] = fwhm_am
            out['ellipticity'][band] = ellipticity

    return out



def group_by_time(arr, obs, delta_t = 86400):
    
    if len(arr) < 1:
        return [[obs]]
    for subarr in arr:
        if any([np.isclose(int(a),int(obs),atol=delta_t) for a in subarr]):
            subarr.append(obs)
            return arr
    arr.append([obs])
    
    return arr



def plot_focus_quasar_fitting_results(data, outdir, xlims=None, interpolate_minima=True):

    missingdata = False
    for source in focus_sources:
        if source not in data.keys():
            missingdata = True
    if missingdata:
        return

    pdata = {**data['PMNJ0210-5101'], **data['PMNJ0522-3628']}
    obsids = np.unique(sorted([int(obsid) for obsid in pdata]))

    grouped_obsids = reduce(group_by_time, obsids, [])
    grouped_obsids = [group for group in grouped_obsids if len(group) == 5]

    fignames = ["BenchF", "BenchE", "FWHM", "Ellipticity"]
    figbenchf = plt.figure("BenchF", figsize=(8,6))
    figbenche = plt.figure("BenchE", figsize=(8,6))
    figfwhm = plt.figure("FWHM", figsize=(8,6))
    figell = plt.figure("Ellipticity", figsize=(8,6))

    ylims = {"BenchF": {90: (20., 33.),
                        150: (20., 33.),
                        220: (20., 33.),
                        'all': (22., 33.)},
             "BenchE": {90: (20., 33.),
                        150: (20., 33.),
                        220: (20., 33.),
                        'all': (22., 33.)},
             "FWHM": {90: (1.3, 1.9),
                      150: (0.9, 1.5),
                      220: (0.7, 1.3),
                      'all': (0.5, 2.1)},
             "Ellipticity": {90: (-0.05, 1.05),
                             150: (-0.05, 1.05),
                             220: (-0.05, 1.05),
                             'all': (-0.05, 1.05)}}
    ylabels = {"BenchF": "Position [mm]",
               "BenchE": "Position [mm]",
               "FWHM": "FWHM [arcmin]", 
               "Ellipticity": "Ellipticity"}
    titles = {"BenchF": "Bench position (along the optical axis) resulting in the minimum FWHM (all)",
              "BenchE": "Bench position (along the optical axis) resulting in the minimum ellipticity (all)",
              "FWHM": "Minimum beam FWHM (all)",
              "Ellipticity": "Minimum beam ellipticity (all)"}
    colors  = {90:'C0', 150:'C1', 220:'C2'}
    markers = {0: 'o', 1: 's', 2: '+', 3: 'x', 4: '<', 5: '>'} 

    minobsid = time_to_obsid(core.G3Time('{}_000000'.format(xlims[0].strftime('%Y%m%d'))))
    maxobsid = time_to_obsid(core.G3Time('{}_000000'.format(xlims[1].strftime('%Y%m%d'))))
    
    nominalbp = []
    nominalts = []
    ndata = {**data['ra0hdec-52.25'],
             **data['ra1h40dec-33.25'], **data['ra5hdec-33.25'], **data['ra12h30dec-33.25']}
    for obsid in ndata.keys():
        nominalts.append(obsid_to_g3time(obsid).time / core.G3Units.s)
        nominalbp.append(ndata[obsid]['NominalBenchPosition'])
    idxordered = np.argsort(nominalts)
    nominalts = np.asarray(nominalts)[idxordered]
    nominalbp = np.asarray(nominalbp)[idxordered]
    nominaldates = np.array([datetime.datetime.utcfromtimestamp(nts) for nts in nominalts])

    for band in [90, 150, 220]:
        timestamps = []
        benchf_min = []
        benche_min = []
        fwhm_min = []
        e_min = []
        
        fig_fwhm_vs_pos = plt.figure("PosVsWidth", figsize=(8,6))
        fig_elli_vs_pos = plt.figure("PosVsEllip", figsize=(8,6))
        
        for counter_g, group in enumerate(grouped_obsids):
            timestamps.append(obsid_to_g3time(int(group[0])).time / core.G3Units.seconds)

            bench = np.array([pdata[obs]['BenchPosAndFittingResults']['bench position'] for obs in group])
            fwhm = np.array([pdata[obs]['BenchPosAndFittingResults']['FWHM']['%dGHz'%band] for obs in group])
            ellipticity = np.array([pdata[obs]['BenchPosAndFittingResults']['ellipticity']['%dGHz'%band] for obs in group])

            try:
                fcoeffs = np.polyfit(bench, fwhm, 2)
                ecoeffs = np.polyfit(bench, ellipticity, 2)
            except ValueError:
                # possibly due to some NaNs
                fcoeffs = [np.nan, np.nan]
                ecoeffs = [np.nan, np.nan]

            bminf = - 0.5*fcoeffs[1]/fcoeffs[0]
            benchf_min.append(bminf)
            bmine = - 0.5*ecoeffs[1]/ecoeffs[0]
            benche_min.append(bmine)
            
            for ydata, coeffs, minloc, xlabel, title \
            in zip([fwhm, ellipticity], [fcoeffs, ecoeffs], [bminf, bmine],
                   ["FWHM [arcmin]", "Ellipticity"], ["FWHM", "Ellipticity"]):
                f = plt.figure(num="parabola", figsize=(8, 6))
                p = f.add_subplot(111)
                p.plot(bench, ydata, linestyle="none", marker=".", markersize=10,
                       color={90: "C0", 150: "C1", 220: "C2"}[band])
                xp = np.linspace(23, 32, num=200)
                yp = coeffs[0]*xp**2 + coeffs[1]*xp + coeffs[2]
                p.plot(xp, yp, color="black")
                p.axvline(minloc, linestyle="dashed", color="black")
                p.set_xlim(left=20, right=33)
                p.set_xticks(np.arange(21, 32), minor=False)
                p.set_xticks(np.linspace(21.5, 31.5, 21), minor=True)
                p.grid(which="major", linestyle="dotted", linewidth=0.75)
                p.grid(which="minor", linestyle="dotted", linewidth=0.45)
                p.tick_params(labelsize=14)
                p.set_xlabel("\nBench position [mm]", fontsize=14)
                p.set_ylabel("{}\n".format(xlabel), fontsize=14)
                p.set_title("{} GHz, Beam {} vs Bench position and the fitting results\nSchedule start time: {}\n".format(
                                band, title, std_processing.obsid_to_g3time(int(group[0]))), fontsize=14)
                f.savefig(os.path.join(outdir, "{}_{:09d}_{:03d}GHz.png".format(title, int(group[0]), band)),
                          bbox_inches="tight")
                plt.close(f)
            
            if interpolate_minima:
                fwhm_min.append(np.poly1d(fcoeffs)(bminf))
                e_min.append(np.poly1d(ecoeffs)(bmine))
            else:
                fwhm_min.append(fwhm.min())
                e_min.append(ellipticity.min())

            plt.figure("PosVsWidth")
            plt.plot(bench, fwhm, label="Schedule "+str(group[0]),
                     marker=markers[counter_g%len(markers)], color=colors[band], alpha=0.65)
            plt.figure("PosVsEllip")
            plt.plot(bench, ellipticity, label="Schedule "+str(group[0]),
                     marker=markers[counter_g%len(markers)], color=colors[band], alpha=0.65)

        plt.figure("PosVsWidth")
        plt.xlim(left=ylims["BenchF"][band][0], right=ylims["BenchF"][band][1])
        plt.ylim(bottom=ylims["FWHM"][band][0], top=ylims["FWHM"][band][1])
        plt.grid()
        plt.legend(loc="upper right")
        plt.xlabel("Position [mm]")
        plt.ylabel("FWHM [arcmin]")
        plt.title("Beam FWHM versus bench position along optical axis, "+str(band)+" GHz (all)")
        plt.tight_layout()
        plt.savefig('{}/focus_fwhm_vs_bench_{}.png'.format(outdir, band)) 
        plt.close()
        plt.figure("PosVsEllip")
        plt.xlim(left=ylims["BenchF"][band][0], right=ylims["BenchF"][band][1])
        plt.ylim(bottom=ylims["Ellipticity"][band][0], top=ylims["Ellipticity"][band][1])
        plt.grid()
        plt.legend(loc="upper right")
        plt.xlabel("Position [mm]")
        plt.ylabel("Ellipticity")
        plt.title("Beam ellipticity versus bench position along optical axis, "+str(band)+" GHz (all)")
        plt.tight_layout()
        plt.savefig('{}/focus_ellip_vs_bench_{}.png'.format(outdir, band))
        plt.close()

        dts = np.array([datetime.datetime.utcfromtimestamp(ts) for ts in timestamps])

        for figname, figarr in zip(fignames, [benchf_min, benche_min, fwhm_min, e_min]):
            plt.figure(figname)
            plot_timeseries(dts, np.array(figarr), band, xlims=xlims, ylims=ylims[figname]['all'], alpha=0.65)
            """
            # useful if we want the bench positions to be printed
            if "Bench" in figname:
                print("=====", figname, "=====\n")
                print("-", band, "GHz")
                print([str(d) for d in dts])
                print(figarr)
                print("\n")
            """

    for figname, figarr in zip(fignames[:2], [benchf_min, benche_min]):
        plt.figure(figname)
        nominaldatenums = mdates.date2num(nominaldates)
        plt.plot(nominaldatenums, nominalbp, color='black', linestyle='dashed', alpha=0.50, label='Nominal pos.')
        plt.legend()

    for figname in fignames:
        plt.figure(figname)
        plt.grid()
        plt.xticks(rotation=25)
        plt.xlabel("observation time (UTC)")
        plt.ylabel(ylabels[figname])
        plt.title(titles[figname])
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}/focus_{}.png'.format(outdir, figname))
        plt.close()
