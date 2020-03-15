import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from spt3g import core, mapmaker, util, pointing



def BeamFWHMandEllipticity(request):
    try:
        data = [fr for fr in core.G3File('{}/calibration/{}/{}.g3' \
                                          .format(request['caldatapath'],
                                                  request['source'],
                                                  request['observation']))]
    except RuntimeError:
        return "Could not find data file."

    colormap = 'viridis'
    textcolor = 'white' 
    figwidth_arcmin = 10.
    orientation = 'horizontal'

    def get_benchpos(frame):
        bp = frame['BenchCommandedPosition']
        bp_zeros = frame['BenchZeros']
        for k in bp.keys():
            bp[k] -= bp_zeros[k]
        return bp


    maps = {}
    benchpos = None
    mapfile = data
    for frame in data:
        if frame.type == core.G3FrameType.Observation:
            benchpos = get_benchpos(frame)
        if 'Id' in frame and 'GHz' in frame['Id']:
            band = re.search(r"[0-9]{2,3}GHz", frame['Id']).group(0)
            maps[band] = frame['T']/frame['Wunpol'].TT

    if benchpos is None:
        core.log_warn("Couldn\'t get bench position from input file. Trying to retrieve original data file")
        try:
            rawdatapath = '{}/{}/{}/0000.g3'.format(
                              request['bolodatapath'],
                              request['source'],
                              request['observation'])
            obsfr = core.G3File(rawdatapath).next()
            benchpos = get_benchpos(obsfr)
        except (AttributeError, RuntimeError):
            core.log_warn("Couldn\'t get bench position from original data file. Proceeding without.")

    if orientation == 'vertical':
        nrows = 3
        ncols = 1
        figsize = (8, 6)
    else:
        nrows = 1
        ncols = 3
        figsize = (8, 6)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, gridspec_kw={'hspace':0.05, 'wspace':0.05})
    
    for ax, band in zip(axes, ['90GHz', '150GHz', '220GHz']):

        res = maps[band].res
        width = int(np.ceil(figwidth_arcmin*core.G3Units.arcmin / res / 2.))
        cy, cx = np.unravel_index(np.nanargmax(maps[band]), maps[band].shape)

        tmap = np.asarray(maps[band])[cy-width : cy+width, cx-width:cx+width]
        tmap /= core.G3Units.mK
        try:
            params = util.fitting.fit_gaussian2d(tmap)
            (amp, x0, y0, sigma_x, sigma_y, theta, height) = params
            px, py = (int(x0), int(y0))

            fit = util.maths.gaussian2d(*params)(*np.indices(tmap.shape))
            magic = 2.*np.sqrt(2.*np.log(2)) # sigma to FWHM conversion
            fwhm_px = magic*np.sqrt(0.5*sigma_x**2 + 0.5*sigma_y**2)
            fwhm_am = fwhm_px * res / core.G3Units.arcmin
            ell = np.sqrt(np.abs(sigma_x**2-sigma_y**2) / np.max((sigma_x**2, sigma_y**2)))

            mask = np.ones(fit.shape, dtype=bool)
            mask[fit/fit.max() > 1e-6] = False # mask source out to -120dB, about 6 sigma
            noise = np.nanstd(tmap[mask])
            SNR = amp / noise
            fitting_success = True
        except:
            fitting_success = False

        im = ax.imshow(tmap, cmap=colormap, aspect='equal')
        if fitting_success:
            ax.contour(fit, [amp/2.], colors=[textcolor], linestyles=['solid'], linewidths=[0.5])
            ax.text(x=px-0.9*width, y=py-0.6*width, s="FWHM: %.2f\'\ne: %.2f"%(fwhm_am, ell), color=textcolor)

            px_am = core.G3Units.arcmin / res
            ct = np.cos(theta)
            st = np.sin(theta)
            sgn = 1. if sigma_x < sigma_y else -1.
            ax.plot([x0 + px_am*ct, x0 - px_am*ct], [y0 + sgn*px_am*st, y0 - sgn*px_am*st], color=textcolor, linewidth=0.5)
            ax.plot([px-0.9*width, px-0.9*width+2*px_am], [py+0.9*width, py+0.9*width], color=textcolor, linewidth=0.5)
            ax.text(x=px-0.8*width+2*px_am, y=py+0.95*width, s='2\'', color=textcolor)
            ax.text(x=px+0.45*width, y=py+0.95*width, s=band, color=textcolor)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

        divider = make_axes_locatable(ax)
        if fitting_success:
            if orientation == 'vertical':
                cax = divider.append_axes('right', size='5%', pad=0.05)
                cb = fig.colorbar(im, cax=cax, orientation='vertical')
                cb.ax.yaxis.set_ticks_position('right')
                cb.ax.yaxis.set_label_position('right')
            else:
                cax = divider.append_axes('bottom', size='5%', pad=0.05)
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                cb.ax.xaxis.set_ticks_position('bottom')
                cb.ax.xaxis.set_label_position('bottom')
            cb.set_label("SNR")
        else:
            ax.text(x=0.5, y=0.5, s=band+'\nSomething\nwent\nwrong...', transform=ax.transAxes, ha='center', va='center')

    if benchpos is not None:
        vpos = pointing.focus.bench2optical(benchpos['x4'], benchpos['y1'], benchpos['z6'])
        fig.suptitle("Beam fitting results from\n%s Observation %s\n(Bench position along the optical axis: %.1fmm)" %(request['source'], request['observation'], vpos[2]))

    return fig
