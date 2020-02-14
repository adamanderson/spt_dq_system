# ============================================================================ #
#  This script is intended to be called by                                     #
#  spt_dq_system/summaryplot/cache_field_maps.py.                              #
#  It can also be run from command line or imported by another script.         #
#                                                                              #
#  The purpose of this script is to make figures showing the                   #
#  analysis results obtained by fields_coadding.py. As a result, the script    #
#  makes many assumptions about the structure and contents of the data it      #
#  receives, which means it is probably not very usable in other contexts.     #
#                                                                              #
#  The script mainly has four parts, each of which defines a pipeline          #
#  module. All four modules make figures, but each module makes a different    #
#  type of figures.                                                            #
#                                                                              #
#  The first module, MakeFiguresForFieldMapsAndWeightMaps, generates           #
#  colormaps from data stored in a map frame, in particular, frame["T"] and    #
#  frame["Wunpol"].TT. It also generates a figure showing a cross section of   #
#  the weight map. Currently, making figures for Q maps, U maps, and the       #
#  associated weight maps is not implemented.                                  #
#                                                                              #
#  The second module, MakeFiguresForTimeVariationsOfMapRelatedQuantities,      #
#  generates figures showing the time variations of some of the map-related    #
#  quantities that were calculated by fields_coadding.py. The quantities are:  #
#    (1) average numbers of detectors flagged by different reasons and         #
#        those not flagged                                                     #
#    (2) medians of the pW/K temperature calibration conversion factors        #
#    (3) medians of fractional changes in detectors' response to the           #
#        calibrator at different elevations                                    #
#    (4) mean weights of TT weight maps                                        #
#    (5) ra offsets and dec offsets three brightesst point sources from        #
#        each sub-field                                                        #
#    (6) averages of ratios of SPT x Planck to Planck x Planck cross spectra   #
#    (7) noise levels of individual or coadded maps                            #
#  These figures will contain data points from only those maps that were used  #
#  to form the coadded maps that were plotted by the first module (those       #
#  coadds are supposed to contain maps from certain time interval such as one  #
#  week, one month, and so on.                                                 #
#                                                                              #
#  The third module, MakeFiguresForTimeEvolutionOfMapRelatedQuantities,        #
#  generates figures showing the time evolution of the noise of the running    #
#  year-to-date coadded noise maps and those showing variances of map values   #
#  in the running coadded signal maps and coadded noise maps. These figures    #
#  are for the 'yearly' section of the webpage and thus will include all       #
#  available data points from certain year instead of just a week, a month,    #
#  and so on.                                                                  #
#                                                                              #
#  The fourth module, MakeFiguresForDistributionsOfMapRelatedQuantities,       #
#  generates figures showing distributions of some of the aforementioned       #
#  quantities by making histograms of them. These quantities are ra offsets,   #
#  dec offsets, ratios of the power spectra, fractional changes in             #
#  responsivity, noise levels of individual maps, and mean weights. These      #
#  histograms are for the 'yearly' section of the webpage and thus will        #
#  include all available data points.                                          #
#                                                                              #
#  At the end of the script, there is a function that constructs a pipeline    #
#  that utilizes these modules.                                                #
#                                                                              #
# ============================================================================ #


import  matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["mathtext.default"] = "regular"
from matplotlib import pyplot

import  os
import  gc
import  sys
import  time
import  glob
import  numpy 
import  argparse
import  logging

from    operator    import  itemgetter
from    spt3g       import  core
from    spt3g       import  mapmaker
from    spt3g       import  coordinateutils
from    spt3g       import  std_processing
from    scipy       import  ndimage
from    scipy       import  linalg
from    scipy       import  optimize

from mpl_toolkits.axes_grid1.inset_locator import inset_axes




# ==============================================================================
# Define functions and pipeline modules needed for the pipeline
# ------------------------------------------------------------------------------


# Define functions related to general plotting utilities
# ------------------------------------------------------------------------------

def get_season_based_on_fields(some_fields):
    
    winter_fields = ["ra0hdec-44.75", "ra0hdec-52.25",
                     "ra0hdec-59.75", "ra0hdec-67.25"]
    summer_fields = ["ra5hdec-24.5" , "ra5hdec-31.5",
                     "ra5hdec-38.5" , "ra5hdec-45.5",
                     "ra5hdec-52.5" , "ra5hdec-59.5"]
    if set(some_fields) <= set(winter_fields):
        return "winter"
    elif set(some_fields) <= set(summer_fields):
        return "summer"
                    

def get_figure_and_plot_objects(w=12.0, h=9.0, dpi=100):
    
    figure_obj = pyplot.figure(figsize=(w, h), dpi=dpi)
    plot_obj   = figure_obj.add_subplot(1, 1, 1)
    return figure_obj, plot_obj


def determine_various_font_sizes(title_fontsize):
    
    labels_fs  = title_fontsize * 0.96
    ticks_fs   = title_fontsize * 0.93
    legends_fs = title_fontsize * 0.90
    return labels_fs, ticks_fs, legends_fs


def determine_various_line_widths(data_linewidth):
    
    ax_h_lw = data_linewidth * 0.75
    ax_v_lw = data_linewidth * 0.75
    major_grid_lw = data_linewidth * 0.40
    minor_grid_lw = data_linewidth * 0.10
    return ax_h_lw, ax_v_lw, major_grid_lw, minor_grid_lw


def set_lims(plot_obj, l, r, b, t):
    
    if l is not None: plot_obj.set_xlim(left=l)
    if r is not None: plot_obj.set_xlim(right=r)
    if b is not None: plot_obj.set_ylim(bottom=b)
    if t is not None: plot_obj.set_ylim(top=t)


def set_ticks(plot_obj, xta, xti, xtl, yta, yti, ytl, ttl_fs,
              dat_lw, xtrot="horizontal"):
    
    lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(ttl_fs)
    axh_lw, axv_lw, gda_lw, gdi_lw = determine_various_line_widths(dat_lw)
    if xta is not None: plot_obj.set_xticks(xta, minor=False)
    if xti is not None: plot_obj.set_xticks(xti, minor=True)
    if yta is not None: plot_obj.set_yticks(yta, minor=False)
    if yti is not None: plot_obj.set_yticks(yti, minor=True)
    if xtl is not None:
        plot_obj.set_xticklabels(
            xtl, fontsize=tck_fs, rotation=xtrot, minor=False)
    if ytl is not None:
        plot_obj.set_yticklabels(ytl, fontsize=tck_fs, minor=False)
    plot_obj.tick_params(
        axis="both", which="major", direction="in", labelsize=tck_fs)
    plot_obj.tick_params(
        axis="both", which="minor", direction="in")
    plot_obj.grid(axis="x", which="major",
                  linestyle="dotted", linewidth=gda_lw)
    plot_obj.grid(axis="x", which="minor",
                  linestyle="dotted", linewidth=gdi_lw)
    plot_obj.grid(axis="y", which="major",
                  linestyle="dotted", linewidth=gda_lw)
    plot_obj.grid(axis="y", which="minor",
                  linestyle="dotted", linewidth=gdi_lw)
    plot_obj.spines["right"].set_visible(False)
    plot_obj.spines["top"].set_visible(False)


def set_ax_labels_and_title(
        plot_obj, xlabel, ylabel, title, ttl_fs, add_legend=True):
    
    lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(ttl_fs)
    if add_legend:
        plot_obj.legend(loc="upper right", fontsize=lgd_fs, framealpha=0.2)
    plot_obj.set_xlabel(xlabel, fontsize=lbl_fs)
    plot_obj.set_ylabel(ylabel, fontsize=lbl_fs)
    plot_obj.set_title(title, fontsize=ttl_fs)


def save_figure_etc(figure_obj, fig_dir, file_name,
                    rect=(0.02, 0.02, 0.97, 0.97), bbox=None):
    
    figure_obj.tight_layout(rect=rect)
    fig_path = os.path.join(fig_dir, file_name)
    figure_obj.savefig(fig_path, bbox_inches=bbox, transparent=False)
    pyplot.close(figure_obj)




# Define pipeline modules and more specific plotting functions
# ------------------------------------------------------------------------------


class MakeFiguresForFieldMapsAndWeightMaps(object):
    
    def __init__(self,
                 map_type=None, map_id=None, coadded_data=True,
                 fig_f=False, fig_w=False, fig_cr=False,
                 rebin_map_before_plotting=False, new_map_resolution=None,
                 smooth_map_with_gaussian=False, gaussian_fwhm=None,
                 custom_vmin_field_map=None,
                 custom_vmax_field_map=None,
                 figure_title_font_size=11,
                 simpler_file_names=True,
                 directory_to_save_figures="",
                 logging_function=logging.info):
        
        self.map_type = map_type
        self.map_id   = map_id
        self.coadded_data = coadded_data
        self.make_figure_for_field_map    = fig_f
        self.make_figure_for_weight_map   = fig_w
        self.make_figure_for_w_map_cr_sec = fig_cr
        self.rebin_map_before_plotting    = rebin_map_before_plotting
        self.res_to_rebin = new_map_resolution
        self.smooth_map_with_gaussian = smooth_map_with_gaussian
        self.gaussian_fwhm = gaussian_fwhm
        self.custom_vmin = custom_vmin_field_map
        self.custom_vmax = custom_vmax_field_map
        self.simpler_file_names = simpler_file_names
        self.ttl_fs = figure_title_font_size
        self.fig_dir = directory_to_save_figures
        self.log = logging_function
        self.obs_fr = None
    
    
    def rebin_map_to_diff_res(self, frame, new_res):
        
        new_res  *= core.G3Units.arcmin
        old_res   = frame["T"].res
        res_ratio = new_res / old_res
        map_parameters = {"x_len"       : int(frame["T"].shape[1]/res_ratio),
                          "y_len"       : int(frame["T"].shape[0]/res_ratio),
                          "res"         : new_res,
                          "proj"        : frame["T"].proj,
                          "alpha_center": frame["T"].alpha_center,
                          "delta_center": frame["T"].delta_center,
                          "coord_ref"   : frame["T"].coord_ref}
        rebin_factor = int(new_res/(frame["T"].res/core.G3Units.arcmin))
        print("The rebin factor is", rebin_factor)
        
        if self.map_type == "T":
            new_field_map  = coordinateutils.FlatSkyMap(**map_parameters)
            new_weight_map = coordinateutils.FlatSkyMap(**map_parameters)
            coordinateutils.reproj_map(
                frame["T"], new_field_map, rebin_factor)
            coordinateutils.reproj_map(
                frame["Wunpol"].TT, new_weight_map, rebin_factor)
            new_frame = core.G3Frame(core.G3FrameType.Map)
            new_frame["T"]         = new_field_map
            new_frame["Wunpol"]    = core.G3SkyMapWeights()
            new_frame["Wunpol"].TT = new_weight_map
        
        return new_frame
    
    
    def get_field_map(self, frame):
        
        if self.map_type == "T":
            t_map = mapmaker.mapmakerutils.remove_weight_t(
                        frame["T"], frame["Wunpol"])
            idx_zero_weights = \
                numpy.where(numpy.asarray(frame["Wunpol"].TT)==0.0)
            numpy.asarray(t_map)[idx_zero_weights] = numpy.nan
            
            return t_map, "T map"
    
    
    def get_weight_map(self, frame):
        
        if self.map_type == "T":
            return frame["Wunpol"].TT, "TT weight map"
    
    
    def get_title_and_file_name_of_figure_for_map(
            self, obs_fr, obs_id_list, mp_ty_str,
            vals_for_stats, nr, map_res, show_res=True):
        
        if self.map_id == "90GHz":
            map_id_for_fig = "95GHz"
        else:
            map_id_for_fig = self.map_id
        
        if (obs_fr is not None) and (obs_id_list is None):
            source  = obs_fr["SourceName"]
            obs_id  = str(obs_fr["ObservationID"])
            date    = str(obs_fr["ObservationStart"]).split(".")[0]
            resol   = str(map_res/core.G3Units.arcmin)+"' map"
            if self.smooth_map_with_gaussian:
                resol += " map smoothed with " + \
                         str(self.gaussian_fwhm) + "' Gaussian"
            title_a = source + "  " + obs_id + " (" + date + ") " + \
                      "   " + map_id_for_fig + " " + mp_ty_str
        
        elif (obs_fr is None) and (obs_id_list is not None):
            season = get_season_based_on_fields(obs_id_list.keys())
            if season == "winter":
                source = "Winter fields"
            elif season == "summer":
                source = "Summer fields"
            obs_ids = []
            for oids_one_field in obs_id_list.values():
                for oid in oids_one_field:
                    obs_ids.append(oid)
            min_id = str(numpy.min(obs_ids))
            max_id = str(numpy.max(obs_ids))
            min_dt = "-".join(str(std_processing.obsid_to_g3time(min_id)).\
                              split(":")[0].split("-")[0:2])
            max_dt = "-".join(str(std_processing.obsid_to_g3time(max_id)).\
                              split(":")[0].split("-")[0:2])
            dt_rng = "between " + min_dt + " and " + max_dt
            obs_id = min_id + " to " + max_id
            el_dic = {"ra0hdec-44.75": "el 0", "ra0hdec-52.25": "el 1",
                      "ra0hdec-59.75": "el 2", "ra0hdec-67.25": "el 3",
                      "ra5hdec-24.5" : "el 0", "ra5hdec-31.5" : "el 1",
                      "ra5hdec-38.5" : "el 2", "ra5hdec-45.5" : "el 3",
                      "ra5hdec-52.5" : "el 4", "ra5hdec-59.5" : "el 5"}
            n_obss = ",  ".join([el_dic[source] + " : " + str(len(obss)) \
                                for source, obss in obs_id_list.items()])
            resol  = str(map_res/core.G3Units.arcmin)+"'"
            if self.smooth_map_with_gaussian:
                resol += " map smoothed with " + \
                         str(self.gaussian_fwhm) + "' Gaussian"
            title_a = source + "   " + map_id_for_fig + \
                      " coadded " + mp_ty_str + "s " + "\n" + \
                      "(data taken " + dt_rng + "," + "\n" + \
                      "{" + n_obss + "})" + "\n"
        
        else:
            raise RuntimeError("Unclear how to build the fig. title!")
        
        if vals_for_stats is None:
            if show_res:
                title_b = "(" + resol + ")" + "\n"
            else:
                title_b = ""
            full_ttl = title_a + title_b
        else:
            if show_res:
                title_b = "(" + resol + ",  "
                if len(resol) > 20:
                    title_b += "\n"
            else:
                title_b = "("
            median  = str(numpy.round(numpy.median(vals_for_stats)))
            pctl_15 = str(numpy.round(numpy.percentile(vals_for_stats, 15), nr))
            pctl_85 = str(numpy.round(numpy.percentile(vals_for_stats, 85), nr))
            title_c = "15th pctl. = " +pctl_15+",  "+ \
                      "85th = " +pctl_85+ ")"
            full_ttl = title_a + title_b + title_c + "\n"
        
        file_nm = source + "-" + obs_id + self.map_id + "_" + mp_ty_str
        file_nm = (file_nm + ".png").replace(" ", "_")
        
        return full_ttl, file_nm
    
    
    def get_colorbar_limits_from_model(self, band, n_maps_ea_fld, distro_stdev):
        
        signal_variances = { "90GHz": 0.00120,
                            "150GHz": 0.00140,
                            "220GHz": 0.00250}
        
        predicted_variances = []
        for sub_field, obs_ids_used in n_maps_ea_fld.items():
            median_var_one_obs  = numpy.median(distro_stdev[sub_field].values())
            median_var_one_obs  = median_var_one_obs ** 2
            median_var_one_obs /= core.G3Units.mK * core.G3Units.mK
            signal_var = signal_variances[band]
            n_map_this_fld = len(obs_ids_used)
            predicted_var_this_fld = \
                signal_var + median_var_one_obs / n_map_this_fld
            predicted_variances.append(predicted_var_this_fld)
            """print("%s %6.4f %d %6.4f %6.4f"
                  %(sub_field, median_var_one_obs, n_map_this_fld,
                    signal_var, predicted_var_this_fld))"""
        
        max_var = numpy.max(predicted_variances)
        max_std = numpy.sqrt(max_var)
        
        return -2.5*max_std, 2.5*max_std
    
    
    def visualize_entire_map(
            self, map_to_view,
            w=12.0, h=8.5, dpi=100, aspect="equal",
            cmap="gray", custom_vmin=None, custom_vmax=None,
            cbar_label="", fig_title="",
            fig_dir="", file_name="map.png"):
        
        figure_obj, plot_obj = get_figure_and_plot_objects(w=w, h=h, dpi=dpi)
        
        if (custom_vmin is None) or (custom_vmax is None):
            pctl_lo = numpy.nanpercentile(numpy.asarray(map_to_view),  1)
            pctl_hi = numpy.nanpercentile(numpy.asarray(map_to_view), 99)
            larger  = numpy.max([numpy.abs(pctl_lo), numpy.abs(pctl_hi)])
            vmin = -1.0 * larger
            vmax =  1.0 * larger
            cbar_label += "  ( 1 to 99 pctl. range)"
        else:
            vmin = custom_vmin
            vmax = custom_vmax
        
        cax = plot_obj.imshow(
                  numpy.asarray(map_to_view),
                  origin="lower", aspect=aspect, interpolation="none",
                  cmap=cmap, vmin=vmin, vmax=vmax)
        
        cbar = figure_obj.colorbar(
                   cax, ax=plot_obj, pad=0.030, shrink=0.75, aspect=20)
        
        lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
        
        plot_obj.tick_params(
            axis="both", which="both",
            bottom=False, top=False, labelbottom=False, labeltop=False,
            left=False, right=False, labelleft=False, labelright=False)
        
        cbar.ax.tick_params(labelsize=tck_fs)
        
        plot_obj.set_title(fig_title, fontsize=self.ttl_fs)
        cbar.ax.set_ylabel(cbar_label, fontsize=lbl_fs)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name,
                        rect=(0.02, 0.0, 1.0, 0.94))
    
    
    def visualize_normalized_weight_map_cross_section(
            self, band, weight_map, n_obss,
            w=12.0, h=9.0, dpi=100,
            xlabel="", ylabel="", fig_title="",
            fig_dir="", file_name="cross_section.png"):
        
        figure_obj, plot_obj = get_figure_and_plot_objects(w=w, h=h, dpi=dpi)
        
        deg = core.G3Units.deg
        
        season = get_season_based_on_fields(n_obss.keys())
        if season  == "winter":
            center_ra   =   0.0
            cen_ra_str  = "0h"
            norm_dec_l  = -70.10
            norm_dec_r  = -69.90
            most_neg_d  = -74.75
            most_pos_d  = -37.25
            center_decs = [-67.25, -59.75, -52.25, -44.75]
            one_fld_uni =   3.25
            one_fld_zr  =   4.25
            n_flds      =   4
        elif season == "summer":
            center_ra   =  75.0
            cen_ra_str  = "5h"
            norm_dec_l  = -62.10
            norm_dec_r  = -61.90
            most_neg_d  = -66.50
            most_pos_d  = -17.50
            center_decs = [-59.5, -52.5, -45.5, -38.5, -31.5, -24.5]
            one_fld_uni =   3.00
            one_fld_zr  =   4.00
            n_flds      =   6
        
        vertical_cr_sec_values = \
            numpy.asarray(weight_map).transpose()[weight_map.shape[1]//2]
        
        norm_dec_r_idx = numpy.unravel_index(
                             weight_map.angle_to_pixel(
                                 center_ra*deg, norm_dec_r*deg),
                             weight_map.shape)[0]
        norm_dec_l_idx = numpy.unravel_index(
                             weight_map.angle_to_pixel(
                                 center_ra*deg, norm_dec_l*deg),
                                 weight_map.shape)[0]
        avg_w_hi_el = \
            numpy.mean(vertical_cr_sec_values[norm_dec_l_idx:norm_dec_r_idx+1])
        if avg_w_hi_el != 0.0:
            hard_to_normalize  = False
            normalized_weights = vertical_cr_sec_values / avg_w_hi_el
        else:
            hard_to_normalize  = True
            normalized_weights = \
                vertical_cr_sec_values / numpy.max(vertical_cr_sec_values)
            max_idx = int(numpy.argmax(vertical_cr_sec_values))
            max_w_ra, max_w_dec = \
                weight_map.pixel_to_angle(weight_map.shape[1]//2, max_idx)
            max_w_dec /= core.G3Units.deg
                
        decs = numpy.linspace(most_neg_d, most_pos_d,
                              (most_pos_d-most_neg_d)*1000+1)
        pids = [weight_map.angle_to_pixel(center_ra*deg, d*deg) for d in decs]
        idcs = [numpy.unravel_index(pid, weight_map.shape)[0] for pid in pids]
        
        focal_plane_hits = numpy.zeros(len(decs))
        for center_dec in center_decs:
            one_sub_field_hits = numpy.zeros(len(decs))
            l_zero = center_dec - one_fld_zr
            l_edge = center_dec - one_fld_uni
            r_edge = center_dec + one_fld_uni
            r_zero = center_dec + one_fld_zr
            idx_lz = numpy.argmin(numpy.absolute(decs - l_zero))
            idx_le = numpy.argmin(numpy.absolute(decs - l_edge))
            idx_re = numpy.argmin(numpy.absolute(decs - r_edge))
            idx_rz = numpy.argmin(numpy.absolute(decs - r_zero))
            for i in range(idx_lz, idx_le+1):
                one_sub_field_hits[i] = \
                    1.0 * (i - idx_lz) / (idx_le - idx_lz)
            for i in range(idx_le, idx_re+1):
                one_sub_field_hits[i] = 1.0
            for i in range(idx_re, idx_rz+1):
                one_sub_field_hits[i] = \
                    1.0 - 1.0 * (i - idx_re) / (idx_rz - idx_re)
            try:
                n_obs = len(n_obss["ra"+cen_ra_str+"dec"+str(center_dec)])
            except KeyError:
                n_obs = 0
            focal_plane_hits += one_sub_field_hits * n_obs
        
        inv_cos = 1.0 / numpy.cos(decs * numpy.pi / 180.0)
        prediction = focal_plane_hits * inv_cos
        hi_el_edge = numpy.where(decs==numpy.mean([norm_dec_l, norm_dec_r]))
        if hard_to_normalize:
            max_w = prediction[numpy.argmin(numpy.abs(decs-max_w_dec))]
        else:
            max_w = prediction[hi_el_edge]
        prediction /= max_w
        
        data_linewidth = 1.50
        color_dict = {"90GHz": "red", "150GHz": "green", "220GHz": "blue"}
        
        plot_obj.plot(normalized_weights, label="Actual data",
                      linewidth=data_linewidth,
                      color=color_dict[band], alpha=0.5)
        plot_obj.plot(idcs, prediction, label="Simple model",
                      linewidth=data_linewidth, color="black")
        
        decs = numpy.linspace(most_neg_d, most_pos_d, 2*n_flds+3)
        pids = [weight_map.angle_to_pixel(
                    center_ra*deg, dec * deg) for dec in decs]
        idcs = [numpy.unravel_index(pid, weight_map.shape)[0] for pid in pids]
        mjridcs = []
        mnridcs = []
        for i, idc in enumerate(idcs[1:-1]):
            if i % 2 == 0:
                mnridcs.append(idc)
            else:
                mjridcs.append(idc)
        
        xlim_left  = idcs[0]
        xlim_right = idcs[-1]
        
        larger_max = numpy.max([numpy.nanmax(normalized_weights),
                                numpy.nanmax(prediction)])
        
        if larger_max > 10.0:
            ytop = 5.2
        else:
            nonzero_weights = \
                normalized_weights[numpy.where(normalized_weights>0.0)]
            max_nonzero = numpy.max(nonzero_weights)
            med_nonzero = numpy.median(nonzero_weights)
            if max_nonzero / med_nonzero > 10.0:
                ytop = med_nonzero * 3.0
            else:
                ytop = larger_max * 1.23
        
        set_lims(plot_obj, xlim_left, xlim_right, -0.02, ytop)
        set_ticks(plot_obj, mjridcs, mnridcs, [str(dec) for dec in center_decs],
                  None, None, None, self.ttl_fs, data_linewidth)
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name,
                        rect=(0.02, 0.02, 0.97, 0.97))
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Observation:
            self.obs_fr = frame
        
        if "FluctuationMetricsIndividualSignalMapsTMapStandardDeviations" \
        in frame.keys():
            self.distro_stdev = \
                frame["FluctuationMetrics"+\
                      "IndividualSignalMaps"+\
                      "TMapStandardDeviations"]
        
        if frame.type == core.G3FrameType.Map:
            if frame["Id"] != self.map_id:
                return
            
            self.log("")
            self.log("* Found a Map frame with ID "+self.map_id+"!")
            self.log("* Figures will be made for this frame.")
            self.log("")
            
            if self.rebin_map_before_plotting:
                self.log("Rebinning the field map to "
                         "%s arcminute resolution ...",
                         self.res_to_rebin)
                frame = self.rebin_map_to_diff_res(frame, self.res_to_rebin)
                self.log("Done.")
                self.log("")
            
            
            if self.make_figure_for_field_map:
                if self.map_type == "T":
                    fd_mp, fd_mp_str = self.get_field_map(frame)
                    
                    self.log("Making a figure for the field map ...")
                    self.log("  Some basic properties of the map:")
                    self.log(" "*4 +"Map shape    : %s", fd_mp.shape)
                    self.log(" "*4 +"Map units    : %s", fd_mp.units)
                    self.log(" "*4 +"Map size     : %s", fd_mp.size)
                    self.log(" "*4 +"x resolution : %s arc minutes",
                             fd_mp.x_res/core.G3Units.arcmin)
                    self.log(" "*4 +"y resolution : %s arc minutes",
                             fd_mp.y_res/core.G3Units.arcmin)
                    
                    fd_mp   /= core.G3Units.mK
                    res      = fd_mp.res
                    cbar_lbl = "\n" + r"$mK_{CMB}$"
                    
                    fd_mp = numpy.asarray(fd_mp)
                    if self.smooth_map_with_gaussian:
                        sigma  = self.gaussian_fwhm * core.G3Units.arcmin
                        sigma /= res
                        sigma /= 2.0 * numpy.sqrt(2.0 * numpy.log(2.0))
                        self.log("  Smoothing the map with a gaussian "
                                 "whose sigma is %s ...", sigma)
                        fd_mp = ndimage.gaussian_filter(
                                    fd_mp, sigma, mode="nearest")
                    
                    self.log("  Preparing the figure title ...")
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                    else:
                        obs_id_list = None
                    if self.smooth_map_with_gaussian:
                        vals_stats = None
                    else:
                        vals_stats = fd_mp[numpy.where(numpy.isfinite(fd_mp))]
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                            self.obs_fr, obs_id_list, fd_mp_str,
                            vals_stats, 3, res, show_res=True)
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                fd_mp_str.replace(" ", "_") + ".png"
                    
                    if (self.custom_vmin is None) and \
                       (self.custom_vmax is None) and \
                       (self.coadded_data)        and \
                       (not self.smooth_map_with_gaussian):
                        self.log("  Deciding on the colorbar limits ...")
                        self.custom_vmin, self.custom_vmax = \
                            self.get_colorbar_limits_from_model(
                                self.map_id, obs_id_list, self.distro_stdev)
                    
                    self.log("  Actually making a figure ...")
                    self.visualize_entire_map(
                        fd_mp, dpi=150, cmap="bwr",
                        custom_vmin=self.custom_vmin,
                        custom_vmax=self.custom_vmax,
                        cbar_label=cbar_lbl, fig_title=fig_ttl, file_name=fl_nm)
                    del fd_mp
                    gc.collect()
                    
                    self.log("Done.")
                    self.log("")
            
            
            if self.make_figure_for_w_map_cr_sec:
                if self.map_type == "T":
                    wt_mp, wt_mp_str = self.get_weight_map(frame)
                    
                    self.log("Making a figure for a cross section of "
                             "the TT weight map ...")
                    
                    map_res = wt_mp.res
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                        field_names = obs_id_list.keys()
                    else:
                        obs_id_list = None
                        field_names = [self.obs_fr["SourceName"]]
                    season = get_season_based_on_fields(field_names)
                    if season == "winter":
                        center_ra = "0h"
                    if season == "summer":
                        center_ra = "5h"
                    
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                            self.obs_fr, obs_id_list, wt_mp_str,
                            None, None, map_res, show_res=False)
                    xlabel  = "\n" + "Declination [degree]"
                    ylabel  = "Normalized weight" + "\n"
                    if self.map_id == "90GHz":
                        map_id_for_fig = "95GHz"
                    else:
                        map_id_for_fig = self.map_id
                    fig_ttl = map_id_for_fig + " " + wt_mp_str + "\n" + \
                              "Cross sectional view " + \
                              "along the RA = "+center_ra+" contour" + "\n"
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                wt_mp_str.replace(" ", "_") + \
                                "_cross_sectional_view.png"
                    
                    self.visualize_normalized_weight_map_cross_section(
                        self.map_id, wt_mp, obs_id_list,
                        xlabel=xlabel, ylabel=ylabel,
                        fig_title=fig_ttl, file_name=fl_nm)
                    
                    del wt_mp
                    gc.collect()
                    
                    self.log("Done.")
                    self.log("")
            
            
            if self.make_figure_for_weight_map:
                if self.map_type == "T":
                    wt_mp, wt_mp_str = self.get_weight_map(frame)
                    
                    self.log("Making a figure for the entire weight map ...")
                    
                    map_res = wt_mp.res
                    wu      = 1.0 / (core.G3Units.mK*core.G3Units.mK)
                    wt_mp   = numpy.asarray(wt_mp/wu)
                    
                    self.log("  Preparing the figure title ...")
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                    else:
                        obs_id_list = None
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                            self.obs_fr, obs_id_list, wt_mp_str,
                            None, 0, map_res, show_res=False)
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                wt_mp_str.replace(" ", "_") + ".png"
                    
                    wt_mp[numpy.where(wt_mp==0.0)] = numpy.nan
                    
                    self.log("  Actually making a figure ...")
                    cbar_lbl = "\n" + "1 / " + r"${mK_{CMB}}^2$"
                    self.visualize_entire_map(
                        wt_mp, dpi=150,
                        custom_vmin=0.0, custom_vmax=numpy.nanmax(wt_mp),
                        cbar_label=cbar_lbl, fig_title=fig_ttl, file_name=fl_nm)
                    
                    del wt_mp
                    gc.collect()
                    
                    self.log("Done.")
                    self.log("")




class MakeFiguresForTimeVariationsOfMapRelatedQuantities(object):
    
    def __init__(self,
                 map_id=None, map_type=None,
                 fig_fs=False, fig_tc=False, fig_rc=False,
                 fig_fl=False, fig_pt=False,
                 fig_rp=False, fig_ns=False,
                 xlim_left=None, xlim_right=None,
                 figure_title_font_size=11,
                 directory_to_save_figures="",
                 logging_function=logging.info):
        
        self.map_type = map_type
        self.map_id   = map_id
        self.all_fields = ["ra0hdec-44.75", "ra0hdec-52.25",
                           "ra0hdec-59.75", "ra0hdec-67.25",
                           "ra5hdec-24.5" , "ra5hdec-31.5" ,
                           "ra5hdec-38.5" , "ra5hdec-45.5" ,
                           "ra5hdec-52.5" , "ra5hdec-59.5"]
        self.make_fig_for_flggg_stats   = fig_fs
        self.make_fig_for_temp_cal_facs = fig_tc
        self.make_fig_for_frac_cal_chng = fig_rc
        self.make_fig_for_fluc_metrics  = fig_fl
        self.make_fig_for_ptg_offsets   = fig_pt
        self.make_fig_for_xspec_ratios  = fig_rp
        self.make_fig_for_noise_levels  = fig_ns
        self.fig_fs_made = False
        self.fig_tc_made = False
        self.fig_rc_made = False
        self.fig_fl_made = False
        self.fig_pt_made = False
        self.fig_rp_made = False
        self.fig_ns_made = False
        self.fig_od_made = False   # * for observation durations
        
        self.obs_ids_of_interest  = None
        self.already_made_figures = False
        
        self.xlim_left  = xlim_left
        self.xlim_right = xlim_right
        self.n_rmarg_el = 4.8
        self.n_rmarg_wf = 6.0
        self.ttl_fs  = figure_title_font_size
        self.el_dict = {"ra0hdec-44.75": "el 0"   , "ra0hdec-52.25": "el 1",
                        "ra0hdec-59.75": "el 2"   , "ra0hdec-67.25": "el 3",
                        "ra5hdec-24.5" : "el 0"   , "ra5hdec-31.5" : "el 1",
                        "ra5hdec-38.5" : "el 2"   , "ra5hdec-45.5" : "el 3",
                        "ra5hdec-52.5" : "el 4"   , "ra5hdec-59.5" : "el 5"}
        self.cl_dict = {"ra0hdec-44.75": "#1f77b4", "ra0hdec-52.25": "#ff7f0e",
                        "ra0hdec-59.75": "#2ca02c", "ra0hdec-67.25": "#d62728",
                        "ra5hdec-24.5" : "#1f77b4", "ra5hdec-31.5" : "#ff7f0e",
                        "ra5hdec-38.5" : "#2ca02c", "ra5hdec-45.5" : "#d62728",
                        "ra5hdec-52.5" : "#9467bd", "ra5hdec-59.5" : "#8c564b"}
        self.waf_cl_dict = {"W172": "#1f77b4", "W174": "#ff7f0e",
                            "W176": "#2ca02c", "W177": "#d62728",
                            "W180": "#9467bd", "W181": "#8c564b",
                            "W188": "#e377c2", "W203": "#7f7f7f",
                            "W204": "#bcbd22", "W206": "#17becf"}
        self.typical_alpha = 0.5
        
        self.ln_wdth = 1.25
        self.ln_styl = "solid"
        self.mrkrsz  = 10.0
        self.typical_xlabel = "\n" + "Time (UTC)"
        
        self.fig_dir = directory_to_save_figures
        self.log = logging_function
    
    
    def frame_has_info(self, frame, keyword):
        
        if any([keyword in key for key in frame.keys()]):
            return True
        else:
            return False
        
    
    def get_data_points_to_plot(self, map_map_double, desired_fields, units):
        
        xs_and_ys = []
        for sub_field, oids_and_data in map_map_double.items():
            if sub_field not in desired_fields:
                continue
            if sub_field not in self.obs_ids_of_interest.keys():
                continue
            relevant_ids = self.obs_ids_of_interest[sub_field]
            for oid in relevant_ids:
                try:
                    xs_and_ys.append((int(oid), oids_and_data[str(oid)]/units))
                except:
                    xs_and_ys.append((int(oid), numpy.nan))
            """for oid, datum in oids_and_data.items():
                if int(oid) not in relevant_ids:
                    continue
                if int(oid) < self.xlim_left:
                    continue
                if int(oid) > self.xlim_right:
                    continue
                xs_and_ys.append((int(oid), datum/units))"""
        
        if len(xs_and_ys) == 0:
            xs = []
            ys = []
            return xs, ys
        
        xs_and_ys = sorted(xs_and_ys, key=itemgetter(0))
        xs = [entry[0] for entry in xs_and_ys]
        ys = [entry[1] for entry in xs_and_ys]
        
        return xs, ys
    
    
    def get_xlims_from_obs_id_range(self, min_obs_id, max_obs_id, nmarg_r):
        
        if (min_obs_id != None) and (max_obs_id != None):
            margin = (max_obs_id - min_obs_id) * 0.05
            xlims_dict = {"left" : min_obs_id - margin * 1.0,
                          "right": max_obs_id + margin * nmarg_r}
        else:
            xlims_dict = {"left": None, "right": None}
        
        return xlims_dict
    
    
    def get_xticks_and_labels_from_obs_id_range(
            self, plot_obj, min_obs_id, max_obs_id, no_hour=False):
        
        if (min_obs_id is None) or (max_obs_id is None):
            xtick_locs   = plot_obj.get_xticks()
            xtick_labels = \
                ["/".join(str(std_processing.obsid_to_g3time(obs_id)).\
                          split(":")[0].split("-")[0:2])              \
                 for obs_id in xtick_locs]
        else:
            xtick_locs_major   = []
            xtick_locs_minor   = []
            xtick_labels_major = []
            
            total_secs = max_obs_id - min_obs_id
            
            current_tick = min_obs_id
            while current_tick <= (max_obs_id + 5):
                tstr  = str(std_processing.obsid_to_g3time(current_tick))
                t_obj = time.strptime(tstr, "%d-%b-%Y:%H:%M:%S.000000000")
                month = str(t_obj.tm_mon)
                date  = str(t_obj.tm_mday)
                hour  = str(t_obj.tm_hour)
                if no_hour:
                    label = month + "/" + date
                else:
                    label = month + "/" + date + "\n" + hour + "h"
                xtick_locs_major.append(current_tick)
                xtick_labels_major.append(label)
                if total_secs < (4 * 24 * 3600 + 10):
                    current_tick += 6 * 3600
                elif total_secs < (8 * 24 * 3600):
                    for hour in [6, 12, 18]:
                        secs_hours = hour * 3600
                        xtick_locs_minor.append(current_tick + secs_hours)
                    current_tick += 24 * 3600
                elif total_secs < (33 * 24 * 3600):
                    for day in [1, 2, 3, 4, 5, 6]:
                        secs_days = day * 24 * 3600
                        xtick_locs_minor.append(current_tick + secs_days)
                    current_tick += 7 * 24 * 3600
                else:
                    for week in [1, 2, 3]:
                        secs_weeks = week * 7 * 24 * 3600
                        xtick_locs_minor.append(current_tick + secs_weeks)
                    current_tick += 4 * 7 * 24 * 3600
        
                xtick_locs_minor = [tick for tick in xtick_locs_minor \
                                    if tick < max_obs_id]
        
        return xtick_locs_major, xtick_labels_major, xtick_locs_minor
    
    
    def indicate_out_of_range_values(
            self, plot_obj, x_values, y_values, ylims_dict, color):
        
        ytop = ylims_dict["top"]
        ybot = ylims_dict["bottom"]
        near_ytop = ybot + 0.90 * (ytop - ybot)
        near_ybot = ybot + 0.15 * (ytop - ybot)
        lbl_fz, tck_fz, lgd_fz = determine_various_font_sizes(self.ttl_fs)
        ah_lw, av_lw, ga_lw, di_lw = determine_various_line_widths(self.ln_wdth)
        note_fs = lgd_fz
        note_lw = av_lw
        note_st = "dashed"
        
        for i, x_value in enumerate(x_values):
            y_value = y_values[i]
            if not numpy.isfinite(y_value):
                plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)
                plot_obj.text(x_value, near_ytop, "NaN", rotation="vertical",
                              color=color, fontsize=note_fs)
            elif y_value >= ytop:
                plot_obj.scatter(x_value, ytop, color=color,
                                 marker=6, s=10*self.mrkrsz)
                """plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)"""
            elif y_value <= ybot:
                plot_obj.scatter(x_value, ybot, color=color,
                                 marker=7, s=10*self.mrkrsz)
                """plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)"""
    
    
    def get_full_fig_title(self, additional_title, no_map_type=False):
        
        if self.map_id == "90GHz":
            map_id_for_fig = "95GHz"
        else:
            map_id_for_fig = self.map_id
        
        if no_map_type:
            prefix = map_id_for_fig + "  -  "
        else:
            prefix = map_id_for_fig + " " + self.map_type + " map" + "  -  "

        return prefix + additional_title
    
    
    def get_full_file_name(self, additional_name):
        
        prefix = self.map_id + "-" + self.map_type + "_" + "map" + "_"
        return prefix + additional_name + ".png"
    
    
    def make_figure_for_observation_durations(self, frame):
        
        if self.fig_od_made or ("ObservationDurations" not in frame):
            return
        
        self.log("")
        self.log("Making a figure for observation duration ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        obs_dur_data = frame["ObservationDurations"]
        
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_el)
        ylims_dict = {"bottom": 0.0, "top": 140.0}
        
        total_obs_time = 0
        for sub_field in sorted(obs_dur_data.keys()):
            oids, lens = self.get_data_points_to_plot(
                             obs_dur_data, [sub_field], core.G3Units.min)
            for one_obs_len in lens:
                total_obs_time += one_obs_len * 60 * 0.965
            label = self.el_dict[sub_field]
            color = self.cl_dict[sub_field]
            
            plot_obj.plot(
                oids, lens, label=label,
                linestyle="None", marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha*1.5)
            
            self.indicate_out_of_range_values(
                plot_obj, oids, lens, ylims_dict, color)
                
        set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                           ylims_dict["bottom"], ylims_dict["top"])
        
        total_elapsed_time = self.xlim_right - self.xlim_left
        obs_eff = 100 * total_obs_time / total_elapsed_time
        
        lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
        text_kargs = {"transform": plot_obj.transAxes,
                      "color"    : "black",
                      "alpha"    : 2.0*self.typical_alpha,
                      "fontsize" : 0.90*lgd_fs,
                      "horizontalalignment": "right"}
        plot_obj.text(0.98, 0.40,
                      "Observing\nefficiency:\n  {:4.1f}%".format(obs_eff)+\
                      "\n(turnarounds\nincluded)",
                      **text_kargs)
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Observation Duration [minute]" + "\n"
        fig_title = "Observation durations and cadence" + "\n"
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        file_name = "observation_durations.png"
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_od_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_flagging_statistics(self, frame):
        
        key_prefix = "FlaggingStatisticsAverageNumbersOfAllBolos"
        info_contained = self.frame_has_info(frame, key_prefix)
        if (not info_contained) or (self.fig_fs_made):
            return
        self.log("")
        self.log("Making a figure for the flagging statistics ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        cl_mk_lbl_dict  = \
            {"BadHk"                  : ["#1f77b4" , "Bad HK"      , "." , 8],
             "Latched"                : ["#ff7f0e" , "Latched"     , "." , 8],
             "Overbiased"             : ["#2ca02c" , "Saturated"   , "." , 8],
             "Oscillating"            : ["#1f77b4" , "Oscillating" , "d" , 6],
             "Glitchy"                : ["#ff7f0e" , "Glitchy"     , "d" , 6],
             "UnphysicalLowVariance"  : ["#2ca02c" , "Low Var."    , "d" , 6],
             "BadCalSn"               : ["#1f77b4" , "Low Cal S/N" , "v" , 6],
             "MissingFluxCalibration" : ["#ff7f0e" , "No Fluxcal"  , "v" , 6],
             "PostCalibrationNaNs"    : ["#2ca02c" , "Has NaNs"    , "v" , 6],
             "BadWeight"              : ["#d62728" , "Bad Weight"  , "v" , 6],
             "Others"                 : ["black"   , "Others"      , "." , 8],
             "TotalNotFlagged"        : ["green"   , "Survived"    , "*" , 9],
             "TotalRemoved"           : ["red"     , "Removed"     , "x" , 8]}
        
        for flag_type in cl_mk_lbl_dict.keys():
            key_name = key_prefix + flag_type
            averages = frame[key_name]
            obs_ids, averages = \
                self.get_data_points_to_plot(averages, self.all_fields, 1.0)
            
            plot_obj.plot(
                obs_ids, averages, label=cl_mk_lbl_dict[flag_type][1],
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=cl_mk_lbl_dict[flag_type][2],
                markersize=cl_mk_lbl_dict[flag_type][3],
                color=cl_mk_lbl_dict[flag_type][0], alpha=self.typical_alpha)
        
        plot_obj.set_yscale("symlog", linthreshy=150,
                            subsy=[2, 3, 4, 5, 6, 7, 8, 9])
        
        ylims_dict = {"bottom": -10, "top": 5000}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, 9.0)
        
        set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                           ylims_dict["bottom"], ylims_dict["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Average (over all scans of an obs.)" + "\n"
        more_title = "Average numbers related to detector flagging" + "\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=True)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "average_numbers_of_flagged_detectors"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_fs_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_calibration_factors(self, frame):
        
        key_prefix = "MediansOfTemperatureCalRelatedFactors"
        info_contained = self.frame_has_info(frame, key_prefix)
        if (not info_contained) or (self.fig_tc_made):
            return
        self.log("")
        self.log("Making figures for the pW/K factors ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        cal_uni = core.G3Units.pW / core.G3Units.K        
        ylims_dict = { "90GHz": {"bottom": -0.16, "top": -0.00},
                      "150GHz": {"bottom": -0.28, "top": -0.00},
                      "220GHz": {"bottom": -0.06, "top": -0.00}}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_wf)
        
        for wafer, color in self.waf_cl_dict.items():
            cal_data = frame[key_prefix + wafer + "PicowattsPerKelvin"]
            obs_ids, cals = self.get_data_points_to_plot(
                              cal_data, self.all_fields, cal_uni)
            
            plot_obj.plot(
                obs_ids, cals, label=wafer,
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha)
            
            self.indicate_out_of_range_values(
                plot_obj, obs_ids, cals, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Median pW/K calibration factor" + "\n"
        more_title = "Medians of temperature calibration factors by wafer"+"\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=True)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "median_temperature_calibration_factors_by_wafer"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_tc_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_responsivity_change(self, frame):
        
        key_prefix = "MedianCalibratorResponse"
        info_contained = self.frame_has_info(frame, key_prefix)
        if (not info_contained) or (self.fig_rc_made):
            return
        self.log("")
        self.log("Making figures for the responsivity changes ...")
                
        # - Make a figure showing variations among wafers
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        ylims_dict = { "90GHz": {"bottom": -5, "top": 25},
                      "150GHz": {"bottom": -5, "top": 25},
                      "220GHz": {"bottom": -2, "top": 10}}
        
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_wf)
        
        for wafer, color in self.waf_cl_dict.items():
            full_key  = key_prefix + wafer + "FractionalChangesTopToBottom"
            frac_chgs = frame[full_key]
            obs_ids, changes = self.get_data_points_to_plot(
                                 frac_chgs, self.all_fields, 0.01)
            
            plot_obj.plot(
                obs_ids, changes, label=wafer,
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha)
            
            self.indicate_out_of_range_values(
                plot_obj, obs_ids, changes, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Median percentage change" + "\n"
        more_title = "Median percentage changes in calibrator response" +"\n"+\
                     "(@ top of the field w.r.t. @ bottom) by wafer" + "\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=True)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "median_cal_resp_percentage_changes_by_wafer"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        
        # - Make a figure showing variations among fields
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        ylims_dict = { "90GHz": {"bottom": -3, "top": 15},
                      "150GHz": {"bottom": -3, "top": 15},
                      "220GHz": {"bottom": -2, "top": 10}}
        
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_el)
        
        for sub_field in self.all_fields:
            full_key  = key_prefix + "AllBolosFractionalChangesTopToBottom"
            frac_chgs = frame[full_key]
            if sub_field not in frac_chgs.keys():
                continue
            obs_ids, changes = self.get_data_points_to_plot(
                                 frac_chgs, [sub_field], 0.01)
            
            color = self.cl_dict[sub_field]
            plot_obj.plot(
                obs_ids, changes, label=self.el_dict[sub_field],
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha)
            
            self.indicate_out_of_range_values(
                plot_obj, obs_ids, changes, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Median percentage change" + "\n"
        more_title = "Median percentage changes in calibrator response" +"\n"+\
                     "(@ top of the field w.r.t. @ bottom) by sub-field" + "\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=True)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "median_cal_resp_percentage_changes_by_field"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_rc_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_map_fluctuation_metrics(self, frame):
        
        key_prefix = "FluctuationMetricsIndividualSignalMaps"
        param_lists = \
            [{"key"   : key_prefix + "MeansOfTTWeights",
              "ylims" : { "90GHz": {"bottom": 20, "top": 100},
                         "150GHz": {"bottom": 30, "top": 170},
                         "220GHz": {"bottom":  1, "top":  15}},
              "yunits": 1 / (core.G3Units.mK * core.G3Units.mK),
              "lnstyl": self.ln_styl,
              "ylabel": "Mean weight  " + \
                        "[1 / " + r"${mK_{CMB}}^2$" + "]" + "\n",
              "title" : "Means of TT weight maps" + "\n",
              "file"  : "mean_of_tt_weight_map_values"},
             {"key"   : key_prefix + "NumbersOfPixelsWithGoodTTWeights",
              "ylims" : { "90GHz": {"bottom": -0.01, "top": 1.10},
                         "150GHz": {"bottom": -0.01, "top": 1.10},
                         "220GHz": {"bottom": -0.01, "top": 1.10}},
              "yunits": 1.0,
              "lnstyl": "None",
              "ylabel": "Fraction" + "\n",
              "title" : "Fractions of pixels having nominal weights" + "\n",
              "file"  : "fractional_coverage_in_tt_weight_maps"}]
        
        
        info_contained = self.frame_has_info(frame, key_prefix)
        if (not info_contained) or (self.fig_fl_made):
            return
        self.log("")
        self.log("Making figures for some fluctuation metrics ...")
        
        for plist in param_lists:
            
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            fluc_data = frame[plist["key"]]
            
            ylims = plist["ylims"]
            if get_season_based_on_fields(fluc_data.keys()) == "summer":
                if "MeansOfTTWeights" in  plist["key"]:
                    ylims = { "90GHz": {"bottom": 20, "top": 160},
                             "150GHz": {"bottom": 30, "top": 250},
                             "220GHz": {"bottom":  1, "top":  25}}
            xlims = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_el)
            
            records_for_later = {}
            for sub_field in sorted(fluc_data.keys()):
                xvals, yvals = self.get_data_points_to_plot(
                                   fluc_data, [sub_field], plist["yunits"])
                if "NumbersOfPixelsWithGoodTTWeights" in plist["key"]:
                    tot_n_pix = {"ra0hdec-44.75": 23493916,
                                 "ra0hdec-52.25": 20249796,
                                 "ra0hdec-59.75": 16651643,
                                 "ra0hdec-67.25": 12786678,
                                 "ra5hdec-24.5" : 11458275,
                                 "ra5hdec-31.5" : 10736705,
                                 "ra5hdec-38.5" :  9854259,
                                 "ra5hdec-45.5" :  8811520,
                                 "ra5hdec-52.5" :  7635120,
                                 "ra5hdec-59.5" :  6370108}
                    yvals = [yval/tot_n_pix[sub_field] for yval in yvals]
                    n_st_0p9 = len([yval for yval in yvals if yval < 0.9])
                    n_tot    = len(yvals)
                    if n_tot == 0:
                        n_rec = 0
                    else:
                        n_rec = int(100.0 * n_st_0p9 / n_tot)
                    records_for_later[sub_field] = n_rec
                if "MeansOfTTWeights" in plist["key"]:
                    records_for_later[sub_field] = numpy.median(yvals)
                
                label = self.el_dict[sub_field]
                color = self.cl_dict[sub_field]
                plot_obj.plot(
                    xvals, yvals, label=label,
                    linestyle=plist["lnstyl"], linewidth=self.ln_wdth,
                    marker=".", markersize=self.mrkrsz,
                    color=color, alpha=self.typical_alpha)
                
                self.indicate_out_of_range_values(
                    plot_obj, xvals, yvals, ylims[self.map_id], color)
            
            set_lims(plot_obj, xlims["left"], xlims["right"],
                     ylims[self.map_id]["bottom"], ylims[self.map_id]["top"])
            
            xtick_locs_major, xtick_labels, xtick_locs_minor = \
                self.get_xticks_and_labels_from_obs_id_range(
                    plot_obj, self.xlim_left, self.xlim_right)
            
            set_ticks(
                plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                None, None, None, self.ttl_fs, self.ln_wdth)
            
            lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
            if "NumbersOfPixelsWithGoodTTWeights" in plist["key"]:
                for counter, sub_field in enumerate(sorted(fluc_data.keys())):
                    good_pc = records_for_later[sub_field]
                    plot_obj.text(0.98, 0.65-counter*0.06,
                                  "{:2d}% maps".format(good_pc),
                                  transform=plot_obj.transAxes,
                                  color=self.cl_dict[sub_field],
                                  alpha=2.0*self.typical_alpha,
                                  fontsize=0.9*lgd_fs,
                                  horizontalalignment="right")
                plot_obj.text(0.98, 0.40, "< 0.9",
                              transform=plot_obj.transAxes,
                              color="black", alpha=2.0*self.typical_alpha,
                              fontsize=0.9*lgd_fs, horizontalalignment="right")
            if "MeansOfTTWeights" in plist["key"]:
                """el3_median = records_for_later["ra0hdec-67.25"]
                el2_to_el3 = records_for_later["ra0hdec-59.75"] / el3_median
                el1_to_el3 = records_for_later["ra0hdec-52.25"] / el3_median
                el0_to_el3 = records_for_later["ra0hdec-44.75"] / el3_median
                text_kargs = {"transform": plot_obj.transAxes,
                              "color"    : "black",
                              "alpha"    : 2.0*self.typical_alpha,
                              "fontsize" : 0.90*lgd_fs,
                              "horizontalalignment": "right"}
                plot_obj.text(0.98, 0.62,
                              "w0 / w3\n= {:4.2f}".format(el0_to_el3),
                              **text_kargs)
                plot_obj.text(0.98, 0.52,
                              "w1 / w3\n= {:4.2f}".format(el1_to_el3),
                              **text_kargs)
                plot_obj.text(0.98, 0.42,
                              "w2 / w3\n= {:4.2f}".format(el2_to_el3),
                              **text_kargs)
                plot_obj.text(0.98, 0.30,
                              "cos3 / cos0\n= 0.54", **text_kargs)
                plot_obj.text(0.98, 0.20,
                              "cos3 / cos1\n= 0.63", **text_kargs)
                plot_obj.text(0.98, 0.10,
                              "cos3 / cos2\n= 0.77", **text_kargs)"""
            
            xlabel = self.typical_xlabel
            ylabel = plist["ylabel"]
            more_title = plist["title"]
            fig_title  = self.get_full_fig_title(more_title)
            
            set_ax_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
            
            more_file_name = plist["file"]
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_fl_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_pointing_discrepancies(self, frame):
        
        key_prefix = "Delta"
        info_contained = self.frame_has_info(frame, key_prefix)
        if (not info_contained) or (self.fig_pt_made):
            return
        self.log("")
        self.log("Making figures for pointing discrepancies ...")
        
        xlims_dict  = self.get_xlims_from_obs_id_range(
                          self.xlim_left, self.xlim_right, 9.1)
        ylims_dict  = {"bottom": -45, "top": 45}
        srccl_dict  = {"1st": "red",    "2nd": "blue",   "3rd": "green"}
        marker_dict = {"1st": ["*", 9], "2nd": ["p", 7], "3rd": [".", 10]}
        
        for sub_field in self.all_fields:
            for coordinate in ["Ra", "Dec"]:
                
                figure_obj, plot_obj = get_figure_and_plot_objects()
                
                all_offsets = {"1st": None, "2nd": None, "3rd": None}
                
                for source_rank in all_offsets.keys():
                    key  = "Delta" + coordinate + "s" + "Of" + \
                           source_rank + "BrightestSourceFromEachSubfield"
                    data = frame[key]
                    
                    obs_ids, offsets = \
                        self.get_data_points_to_plot(
                            data, [sub_field], core.G3Units.arcsec)
                    all_offsets[source_rank] = offsets
                
                avgs = numpy.nanmean(list(all_offsets.values()), axis=0)
                
                for source_rank, offsets in all_offsets.items():
                    plot_obj.plot(
                        obs_ids, offsets,
                        label="Source. "+source_rank[0],
                        linewidth=0.2*self.ln_wdth,
                        marker=marker_dict[source_rank][0],
                        markersize=marker_dict[source_rank][1],
                        color=srccl_dict[source_rank],
                        alpha=self.typical_alpha)
                    
                    self.indicate_out_of_range_values(
                        plot_obj, obs_ids, offsets, ylims_dict,
                        srccl_dict[source_rank])
                
                plot_obj.plot(
                    obs_ids, avgs, label="Avg. offsets",
                    linestyle=self.ln_styl, linewidth=self.ln_wdth,
                    color="black", alpha=self.typical_alpha)
                
                set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                                   ylims_dict["bottom"], ylims_dict["top"])
                
                text_tr = 0.67
                separat = 0.06
                counter = 0
                lbl_fs, tck_fs, lgd_fs = \
                    determine_various_font_sizes(self.ttl_fs)
                for source_rank, offsets in all_offsets.items():
                    avg = numpy.nanmean(offsets)
                    plot_obj.text(0.98, text_tr-counter*separat,
                                  "Src. {}, avg.: {:+05.1f}".\
                                  format(source_rank[0], avg),
                                  transform=plot_obj.transAxes,
                                  horizontalalignment="right",
                                  verticalalignment="top",
                                  color=srccl_dict[source_rank],
                                  fontsize=lgd_fs)
                    std = numpy.nanstd(offsets)
                    plot_obj.text(0.98, text_tr-(counter+4)*separat,
                                  "Src. {}, std.: {:04.1f}".\
                                  format(source_rank[0], std),
                                  transform=plot_obj.transAxes,
                                  horizontalalignment="right",
                                  verticalalignment="top",
                                  color=srccl_dict[source_rank],
                                  fontsize=lgd_fs)
                    counter += 1
                
                avg_avgs = numpy.nanmean(avgs)
                std_avgs = numpy.nanstd(avgs)
                plot_obj.text(0.98, text_tr-8*separat,
                              "Avg., avg.: {:+05.1f}".format(avg_avgs),
                              transform=plot_obj.transAxes,
                              horizontalalignment="right",
                              verticalalignment="top",
                              color="black",
                              fontsize=lgd_fs)
                plot_obj.text(0.98, text_tr-9*separat,
                              "Avg., std.: {:04.1f}".format(std_avgs),
                              transform=plot_obj.transAxes,
                              horizontalalignment="right",
                              verticalalignment="top",
                              color="black",
                              fontsize=lgd_fs)
                
                xtick_locs_major, xtick_labels, xtick_locs_minor = \
                    self.get_xticks_and_labels_from_obs_id_range(
                        plot_obj, self.xlim_left, self.xlim_right)
                
                set_ticks(plot_obj,
                          xtick_locs_major, xtick_locs_minor, xtick_labels,
                          None, None, None, self.ttl_fs, self.ln_wdth)
                
                xlabel = self.typical_xlabel
                if coordinate == "Ra":
                    ylabel = "(Measured R.A. - True R.A.) " + \
                             r"$\times$" + " cos(True Dec.)"+ ' ["]' + "\n"
                else:
                    ylabel = "Measured Dec. - True Dec." + ' ["]' + "\n"
                more_title = coordinate + " offsets of point sources " + \
                             "in " + sub_field + "\n"
                
                fig_title = self.get_full_fig_title(more_title)
                
                set_ax_labels_and_title(
                    plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
                
                more_file_name = "delta_" + coordinate + \
                                 "s_from_point_sources_in_" + sub_field
                file_name = self.get_full_file_name(more_file_name)
                
                save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_pt_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_xspec_ratios(self, frame):
        
        rat_key = "AveragesOfRatiosOfSPTxPlancktoPlckxPlckTTspectra"
        info_contained = self.frame_has_info(frame, rat_key)
        if (not info_contained) or (self.fig_rp_made):
            return
        self.log("")
        self.log("Making a figure for the ratios of the SPTxPlanck / PxP ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        ylims_dict = { "90GHz": {"bottom": 0.5, "top": 2.0},
                      "150GHz": {"bottom": 0.4, "top": 1.6},
                      "220GHz": {"bottom": 0.1, "top": 1.9}}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_el)
        
        data = frame[rat_key]
        for sub_field in data.keys():
            obs_ids, ratios = self.get_data_points_to_plot(
                                 data, [sub_field], 1.0)
            
            color = self.cl_dict[sub_field]
            plot_obj.plot(
                obs_ids, ratios, label=self.el_dict[sub_field],
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha)
            
            self.indicate_out_of_range_values(
                plot_obj, obs_ids, ratios, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = self.typical_xlabel
        ylabel = "Average of S x P / P x P" + "\n" +\
                 "in the ell range [750, 1250]" + "\n"
        more_title = "Averages of ratios of" +"\n" + \
                     "SPT x Planck to Planck x Planck spectra" + "\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=False)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "averages_of_power_spectra_ratios"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_rp_made = True
        self.log("Done.")
        self.log("")
    
    
    def make_figure_for_noise_levels(self, frame):
        
        noise_key = "NoiseLevelsFromIndividualTMaps"
        info_contained = self.frame_has_info(frame, noise_key)
        if (not info_contained) or (self.fig_ns_made):
            return
        self.log("")
        self.log("Making a figure for noise levels ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        season = get_season_based_on_fields(frame[noise_key].keys())
        if season  == "winter":
            ylims_dict = { "90GHz": {"bottom":  90, "top": 210},
                          "150GHz": {"bottom":  90, "top": 210},
                          "220GHz": {"bottom": 300, "top": 700}}
        if season == "summer":
            ylims_dict = { "90GHz": {"bottom":  60, "top": 220},
                          "150GHz": {"bottom":  60, "top": 220},
                          "220GHz": {"bottom": 200, "top": 740}}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, self.n_rmarg_el)
        
        records_for_later = {}
        data = frame[noise_key]
        for sub_field in sorted(data.keys()):
            nuni = core.G3Units.uK * core.G3Units.arcmin
            obs_ids, noises = self.get_data_points_to_plot(
                                  data, [sub_field], nuni)
            records_for_later[sub_field] = numpy.nanmedian(noises)**2
            
            color = self.cl_dict[sub_field]
            plot_obj.plot(
                obs_ids, noises, label=self.el_dict[sub_field],
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=self.typical_alpha)
            
            self.indicate_out_of_range_values(
                plot_obj, obs_ids, noises, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
        """el3_value  = records_for_later["ra0hdec-67.25"]
        el2_to_el3 = records_for_later["ra0hdec-59.75"] / el3_value
        el1_to_el3 = records_for_later["ra0hdec-52.25"] / el3_value
        el0_to_el3 = records_for_later["ra0hdec-44.75"] / el3_value
        text_kargs = {"transform": plot_obj.transAxes,
                      "color"    : "black",
                      "alpha"    : 2.0*self.typical_alpha,
                      "fontsize" : 0.90*lgd_fs,
                      "horizontalalignment": "right"}
        plot_obj.text(0.98, 0.62,
                      "(n0/n3)^2\n= {:4.2f}".format(el0_to_el3),
                      **text_kargs)
        plot_obj.text(0.98, 0.52,
                      "(n1/n3)^2\n= {:4.2f}".format(el1_to_el3),
                              **text_kargs)
        plot_obj.text(0.98, 0.42,
                      "(n2/n3)^2\n= {:4.2f}".format(el2_to_el3),
                      **text_kargs)
        plot_obj.text(0.98, 0.30,
                      "cos0/cos3\n= 1.84", **text_kargs)
        plot_obj.text(0.98, 0.20,
                      "cos1/cos3\n= 1.58", **text_kargs)
        plot_obj.text(0.98, 0.10,
                      "cos2/cos3\n= 1.30", **text_kargs)"""
        
        xlabel = self.typical_xlabel
        ylabel = "Average  " + r"$\sqrt{C_l}$" + "  " + "in the " + "\n" + \
                 "ell range [3000, 5000]  " + r"$[{\mu}K \cdot arcmin]$" + "\n"
        more_title = "Noise (+ a little signal) of individual maps" + "\n"
        fig_title  = self.get_full_fig_title(more_title, no_map_type=False)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "noise_levels_from_individual_maps"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.fig_ns_made = True
        self.log("Done.")
        self.log("")
    
    
    def __call__(self, frame):
        
        if "Id" not in frame:
            return
        
        if frame["Id"] != self.map_id:
            return
                
        if "CoaddedObservationIDs" in frame:
            self.obs_ids_of_interest = frame["CoaddedObservationIDs"]
        
        self.make_figure_for_observation_durations(frame)
        
        if self.make_fig_for_flggg_stats:
            self.make_figure_for_flagging_statistics(frame)
        
        if self.make_fig_for_temp_cal_facs:
            self.make_figure_for_calibration_factors(frame)
        
        if self.make_fig_for_frac_cal_chng:
            self.make_figure_for_responsivity_change(frame)
        
        if self.make_fig_for_fluc_metrics:
            self.make_figure_for_map_fluctuation_metrics(frame)
        
        if self.make_fig_for_ptg_offsets:
            self.make_figure_for_pointing_discrepancies(frame)
        
        if self.make_fig_for_xspec_ratios:
            self.make_figure_for_xspec_ratios(frame)
        
        if self.make_fig_for_noise_levels:
            self.make_figure_for_noise_levels(frame)




class MakeFiguresForTimeEvolutionOfMapRelatedQuantities(object):
    
    def __init__(self,
                 map_id=None, map_type=None,
                 fig_fl=False, fig_ns=False,
                 figure_title_font_size=11,
                 directory_to_save_figures=".",
                 logging_function=logging.info):
        
        self.log = logging_function
        self.map_id   = map_id
        self.map_type = map_type
        self.obs_id_list      = None
        self.stddevs_individs = None
        
        self.make_figs_for_fluc_metrics = fig_fl
        self.make_figs_for_noise_levels = fig_ns
        self.fig_an_made = False
        self.fig_rn_made = False
        self.fig_oa_made = False
        self.fig_fl_made = False
        
        self.ttl_fs  = figure_title_font_size
        self.el_dict = {"ra0hdec-44.75": "el 0"   , "ra0hdec-52.25": "el 1",
                        "ra0hdec-59.75": "el 2"   , "ra0hdec-67.25": "el 3",
                        "ra5hdec-24.5" : "el 0"   , "ra5hdec-31.5" : "el 1",
                        "ra5hdec-38.5" : "el 2"   , "ra5hdec-45.5" : "el 3",
                        "ra5hdec-52.5" : "el 4"   , "ra5hdec-59.5" : "el 5"}
        self.cl_dict = {"ra0hdec-44.75": "#1f77b4", "ra0hdec-52.25": "#ff7f0e",
                        "ra0hdec-59.75": "#2ca02c", "ra0hdec-67.25": "#d62728",
                        "ra5hdec-24.5" : "#1f77b4", "ra5hdec-31.5" : "#ff7f0e",
                        "ra5hdec-38.5" : "#2ca02c", "ra5hdec-45.5" : "#d62728",
                        "ra5hdec-52.5" : "#9467bd", "ra5hdec-59.5" : "#8c564b"}
        self.ln_wdth = 1.25
        self.ln_styl = "solid"
        self.mrkrsz  = 10.0
        self.fig_dir = directory_to_save_figures
        self.lbl_fs, self.tck_fs, self.lgd_fs = \
            determine_various_font_sizes(self.ttl_fs)
    
    
    def get_full_fig_title(self, additional_title, no_map_type=False):
        
        if self.map_id == "90GHz":
            map_id_for_fig = "95GHz"
        else:
            map_id_for_fig = self.map_id
        
        if no_map_type:
            prefix = map_id_for_fig + "  -  "
        else:
            prefix = map_id_for_fig + " " + self.map_type + " map" + "  -  "

        return prefix + additional_title
    
    
    def get_full_file_name(self, additional_name):
        
        prefix = self.map_id + "-" + self.map_type + "_" + "map" + "_"
        return prefix + additional_name + ".png"
    
    
    def make_figure_for_all_noise_levels(self, frame):
        
        self.log("")
        self.log("Making a figure for noise levels of coadded maps ...")
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        noise_data = frame["NoiseLevelsFromCoaddedTMaps"]
        
        plot_obj.set_xscale("log")
        plot_obj.set_yscale("log")
        season = get_season_based_on_fields(noise_data.keys())
        if season == "winter":
            ylims_dict = { "90GHz": {"bottom":  5, "top": 200},
                          "150GHz": {"bottom":  5, "top": 200},
                          "220GHz": {"bottom": 15, "top": 600}}
        elif season == "summer":
            ylims_dict = { "90GHz": {"bottom":  5, "top": 200},
                          "150GHz": {"bottom":  5, "top": 200},
                          "220GHz": {"bottom": 15, "top": 600}}
        
        max_n_obss = 0
        for sub_field, noise_dict in noise_data.items():
            n_obss = len(noise_dict.keys())
            if n_obss > max_n_obss:
                max_n_obss = n_obss
        xlims_dict = {"left" : 0.9,
                      "right": 10**numpy.ceil(numpy.log10(max_n_obss))}
        
        n_excluded = {}
        bad_ids = frame["IgnoredObservationIDs"]
        for sub_field in self.el_dict.keys():
            if sub_field in bad_ids.keys():
                n_excluded[sub_field] = len(bad_ids[sub_field])
            """else:
                n_excluded[sub_field] = 0"""
        
        for sub_field, obs_ids in self.obs_id_list.items():
            noise_units  = core.G3Units.uK * core.G3Units.arcmin
            noise_levels = [noise_data[sub_field][str(obs_id)] / noise_units \
                            for obs_id in obs_ids \
                            if str(obs_id) in noise_data[sub_field].keys()]
            x_data = range(1, len(noise_levels)+1)
            
            label = self.el_dict[sub_field]
            color = self.cl_dict[sub_field]
            
            plot_obj.plot(
                x_data, noise_levels, label=label,
                linestyle="None",
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=0.5)
            
            curr_nois = noise_levels[-1]
            ylocs_dict  = {"ra0hdec-44.75": 0.93,
                           "ra0hdec-52.25": 0.86,
                           "ra0hdec-59.75": 0.79,
                           "ra0hdec-67.25": 0.72,
                           "ra5hdec-24.5" : 0.93,
                           "ra5hdec-31.5" : 0.86,
                           "ra5hdec-38.5" : 0.79,
                           "ra5hdec-45.5" : 0.72,
                           "ra5hdec-52.5" : 0.65,
                           "ra5hdec-59.5" : 0.58}
            plot_obj.text(0.50, ylocs_dict[sub_field],
                          "YTD noise: {:05.2f}".format(curr_nois),
                          transform=plot_obj.transAxes,
                          horizontalalignment="left",
                          color=color, fontsize=self.lgd_fs)
            
            if len(x_data) > 20:
                x_for_lls = numpy.log(numpy.asarray(x_data[11:]))
                y_for_lls = numpy.log(numpy.asarray(noise_levels[11:]))
            else:
                x_for_lls = numpy.log(numpy.asarray(x_data))
                y_for_lls = numpy.log(numpy.asarray(noise_levels))
            design_matrix = numpy.asarray(x_for_lls)\
                            [:, numpy.newaxis]**[1.0, 0.0]
            parameters, residues, rank, singular_values = \
                linalg.lstsq(design_matrix, numpy.asarray(y_for_lls))
            power = parameters[0]
            const = parameters[1]
            x_for_plot = numpy.arange(1, xlims_dict["right"], 2)
            y_for_plot = numpy.exp(const) * numpy.power(x_for_plot, power)
            plot_obj.plot(x_for_plot, y_for_plot,
                          color=color, alpha=0.8, linewidth=0.4*self.ln_wdth)
            explanation = "Slope: " + str(power)[0:6]
            ylocs_dict  = {"ra0hdec-44.75": 0.25,
                           "ra0hdec-52.25": 0.18,
                           "ra0hdec-59.75": 0.11,
                           "ra0hdec-67.25": 0.04,
                           "ra5hdec-24.5" : 0.39,
                           "ra5hdec-31.5" : 0.32,
                           "ra5hdec-38.5" : 0.25,
                           "ra5hdec-45.5" : 0.18,
                           "ra5hdec-52.5" : 0.11,
                           "ra5hdec-59.5" : 0.04}
            plot_obj.text(0.03, ylocs_dict[sub_field],
                          explanation,
                          transform=plot_obj.transAxes,
                          horizontalalignment="left",
                          color=color, fontsize=self.lgd_fs)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
       
        xtick_locs_major = None
        xtick_labels     = None
        xtick_locs_minor = None
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs, self.ln_wdth)
        
        xlabel = "\n" + "Number of pairs of difference map added"
        ylabel = "Average  " + r"$\sqrt{C_l}$" + "  " + "in the " + "\n" + \
                 "ell range [3000, 5000]  " + r"$[{\mu}K \cdot arcmin]$" + "\n"
        
        more_title = "Noise in the running coadded maps\n"
        more_info  = []
        for sub_field in sorted(n_excluded.keys()):
            alias = self.el_dict[sub_field]
            n_bad = n_excluded[sub_field]
            more_info.append(alias + ": " + str(n_bad))
        more_info   = "(excluded maps: {" + \
                      ",  ".join(more_info) + "})\n"
        more_title += more_info
        fig_title = self.get_full_fig_title(more_title)
        
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "noise_levels_from_running_coadds"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.log("Done.")
        self.log("")
    
    
    def make_figures_for_recent_noise_levels(self, frame):
        
        self.log("")
        self.log("Making figures for the time evolution of the noise in the ")
        self.log("coadded noise maps due to recently added maps ...")
        
        noise_data = frame["NoiseLevelsFromCoaddedTMaps"]
        
        for sub_field, obs_ids_this_sf in self.obs_id_list.items():
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            l_ids = []
            r_ids = []
            for obs_id in obs_ids_this_sf:
                if str(obs_id) in noise_data[sub_field].keys():
                    r_ids.append(obs_id)
                else:
                    l_ids.append(obs_id)
            
            all_obs_ids = []
            for i in range(numpy.min([len(l_ids), len(r_ids)])):
                all_obs_ids.append(l_ids[i])
                all_obs_ids.append(r_ids[i])
            
            all_noises = []
            for obs_id in all_obs_ids:
                if str(obs_id) in noise_data[sub_field].keys():
                    noise  = noise_data[sub_field][str(obs_id)]
                    noise /= core.G3Units.uK * core.G3Units.arcmin
                    all_noises.append(noise)
                else:
                    all_noises.append(numpy.nan)
            
            if len(all_obs_ids) > 50:
                all_obs_ids = all_obs_ids[-50:]
                all_noises  = all_noises[-50:]
            dummy_x_data = numpy.arange(len(all_obs_ids)) + 1
            
            color = self.cl_dict[sub_field]
            
            plot_obj.plot(dummy_x_data, all_noises,
                          linestyle=self.ln_styl, linewidth=self.ln_wdth,
                          marker=".", markersize=self.mrkrsz,
                          color=color)
            
            max_noise = numpy.nanmax(all_noises)
            min_noise = numpy.nanmin(all_noises)
            diff      = max_noise - min_noise
            ylim_top  = max_noise + 0.05 * diff
            near_top  = max_noise + 0.01 * diff
            
            set_lims(plot_obj, 0, dummy_x_data[-1]+1, None, ylim_top)
            
            lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
            for index, noise in enumerate(all_noises):
                if not numpy.isfinite(noise):
                    plot_obj.text(dummy_x_data[index], near_top,
                                  "N/A", color=color, fontsize=lgd_fs*0.7,
                                  rotation="vertical",
                                  verticalalignment="center",
                                  horizontalalignment="center")
            
            def get_month_day_hour(obs_id):
                full_tstr = str(std_processing.obsid_to_g3time(obs_id))
                month     = full_tstr.split("-")[1]
                date      = full_tstr.split("-")[0]
                hour      = full_tstr.split(":")[1]
                return "{}/{} {}h".format(month, date, hour)
            
            xtick_labels = ["{} ({})".format(
                            str(obs_id), get_month_day_hour(obs_id)) \
                            for obs_id in all_obs_ids]
            
            set_ticks(plot_obj, dummy_x_data, [], xtick_labels,
                      None, None, None, self.ttl_fs, self.ln_wdth,
                      xtrot="vertical")
            
            xlabel = ""
            ylabel = "Noise " + r"$[{\mu}K \cdot arcmin]$"
            more_title = "Noise due to 50 recently added maps of " + sub_field 
            fig_title  = self.get_full_fig_title(more_title)
            
            set_ax_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs,
                add_legend=False)
            
            more_file_name = "recent_noise_levels_from_observations_of_{}".\
                             format(sub_field)
            file_name = self.get_full_file_name(more_file_name)
            save_figure_etc(figure_obj, self.fig_dir, file_name)
            
        self.log("Done.")
        self.log("")
    
    
    def make_figures_for_order_of_addition(self, frame):
        
        self.log("")
        self.log("Making figures for the order in which maps were added ...")
        
        noise_data = frame["NoiseLevelsFromCoaddedTMaps"]
        
        for sub_field, obs_ids_this_sf in self.obs_id_list.items():
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            l_ids = []
            r_ids = []
            for obs_id in obs_ids_this_sf:
                if str(obs_id) in noise_data[sub_field].keys():
                    r_ids.append(obs_id)
                else:
                    l_ids.append(obs_id)
            
            all_obs_ids = []
            for i in range(numpy.min([len(l_ids), len(r_ids)])):
                all_obs_ids.append(l_ids[i])
                all_obs_ids.append(r_ids[i])
            
            plot_obj.plot(all_obs_ids,
                          linestyle="None",
                          marker=".", markersize=self.mrkrsz,
                          color=self.cl_dict[sub_field])
            
            min_oid = numpy.min(all_obs_ids)
            max_oid = numpy.max(all_obs_ids)
            
            ylims_dict = \
                MakeFiguresForTimeVariationsOfMapRelatedQuantities.\
                get_xlims_from_obs_id_range(
                    "dummy", min_oid, max_oid, 0.0)
            ylims_dict["bottom"] = ylims_dict.pop("left")
            ylims_dict["top"]    = ylims_dict.pop("right")
            
            set_lims(plot_obj, 0, len(obs_ids_this_sf)+1,
                               ylims_dict["bottom"], ylims_dict["top"])
            
            ytick_locs_major, ytick_labels, ytick_locs_minor = \
                MakeFiguresForTimeVariationsOfMapRelatedQuantities.\
                get_xticks_and_labels_from_obs_id_range(
                    "dummy", plot_obj, min_oid, max_oid, no_hour=True)
            
            set_ticks(plot_obj, None, None, None,
                      ytick_locs_major, ytick_locs_minor, ytick_labels,
                      self.ttl_fs, self.ln_wdth)
            
            xlabel = "\n" + "Number of observations added"
            ylabel = "Date of observation" + "\n"
            
            more_title = "Order in which observations of " + \
                         sub_field + " were added" + "\n"
            fig_title = self.get_full_fig_title(more_title)
            
            set_ax_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs,
                add_legend=False)
            
            more_file_name = "order_of_addition_of_maps_of_" + sub_field
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
            
        self.log("Done.")
        self.log("")
    
    
    def make_figures_for_variances_of_maps(self, frame):
        
        self.log("")
        self.log("Making figures for variances of coadded maps ...")
        
        stddevs_sigs = \
            frame["FluctuationMetricsCoaddedSignalMapsTMapStandardDeviations"]
        stddevs_nois = \
            frame["FluctuationMetricsCoaddedNoiseMapsTMapStandardDeviations"]
        vu = core.G3Units.uK * core.G3Units.uK
        
        for sub_field, obs_ids_this_sf in self.obs_id_list:
            figure_obj, plot_obj = get_figure_and_plot_objects(w=15.0, h=9.0)
            
            median_individ_vars = numpy.nanmedian(
                numpy.asarray(self.stddevs_individs[sub_field].values())**2)/vu
            median_individ_vars /= 2
            
            y_data = []
            for obs_id in obs_ids_this_sf:
                if str(obs_id) in stddevs_sigs[sub_field].keys():
                    y_data.append(stddevs_sigs[sub_field][str(obs_id)]**2/vu)
            x_data = numpy.arange(1, len(y_data)+1, 1)
            plot_obj.plot(x_data, y_data,
                          label="Variances of coadded signal maps",
                          linestyle="None", marker=".", markersize=8,
                          color="#1f77b4", alpha=0.5)
            
            if len(y_data) > 40:
                x_for_fit = x_data[20:]
                y_for_fit = y_data[20:]
            else:
                x_for_fit = x_data
                y_for_fit = y_data
            def sig_plu_noi_mod(n_mp, sig):
                return sig + median_individ_vars/(n_mp)
            popt, pcov = optimize.curve_fit(
                sig_plu_noi_mod, x_for_fit, y_for_fit, p0=[y_data[-1]])
            
            x_for_show = numpy.arange(1, 250, 1)
            y_for_show = sig_plu_noi_mod(x_for_show, popt[0])
            plot_obj.plot(x_for_show, y_for_show,
                          label="Curve fit y = (med. blue distro.) "+\
                                 "/ x + b"+\
                                 "\n(b = {:4.2e})".format(popt[0]),
                                 linewidth=1.0, color="#1f77b4", alpha=0.5)
            
            
            y_data = []
            for obs_id in obs_ids_this_sf:
                if str(obs_id) in stddevs_nois[sub_field].keys():
                    y_data.append(stddevs_nois[sub_field][str(obs_id)]**2/vu)
            x_data = numpy.arange(1, len(y_data)+1, 1)
            plot_obj.plot(x_data, y_data,
                          label="Variances of coadded noise maps",
                          linestyle="None", marker=".", markersize=8,
                          color="#ff7f0e", alpha=0.5)
            
            if len(y_data) > 40:
                x_for_fit = x_data[20:]
                y_for_fit = y_data[20:]
            else:
                x_for_fit = x_data
                y_for_fit = y_data
            def noi_mod(n_mp, power):
                return median_individ_vars/n_mp**power
            popt, pcov = optimize.curve_fit(
                noi_mod, x_for_fit, y_for_fit, p0=[1.0])
            
            x_for_show = numpy.arange(1, 250, 1)
            y_for_show = noi_mod(x_for_show, popt[0])
            plot_obj.plot(x_for_show, y_for_show,
                          label="Curve fit y = (med. blue distro.) "+\
                                "/ "+r"$x^p$"+\
                                "\n(p = {:4.2e})".format(popt[0]),
                                linewidth=1.0, color="#ff7f0e", alpha=0.5)
            
            
            y_data = numpy.asarray(
                         self.stddevs_individs[sub_field].values())**2/vu
            y_data /= 2
            
            bps = plot_obj.boxplot(
                      y_data, whis=[0, 100], sym="", positions=[1], widths=0.20,
                      patch_artist=True,
                      boxprops=dict(facecolor="blue", alpha=0.1),
                      medianprops=dict(linestyle="none"),
                      capprops=dict(color="blue", alpha=0.1, linewidth=2.0),
                      whiskerprops=dict(color="blue", alpha=0.1))
            
            plot_obj.set_xscale("log")
            plot_obj.set_yscale("log")
            plot_obj.set_xlim(left=0.8, right=260)
            if self.map_id in ["90GHz", "150GHz"]:
                plot_obj.set_ylim(bottom=2e2, top=2.5e5)
            if self.map_id in ["220GHz"]:
                plot_obj.set_ylim(bottom=2e3, top=2.5e6)
            
            plot_obj.tick_params(axis="both", which="major",
                                 direction="in", labelsize=14)
            
            handles, labels = plot_obj.get_legend_handles_labels()
            handles += [bps["boxes"][0]]
            labels  += ["Variances of individual maps / 2"]
            plot_obj.legend(handles, labels, loc="upper right",
                            fontsize=15, framealpha=0.8)
            
            plot_obj.grid(axis="both", which="both",
                          linestyle="dashed", linewidth=0.05, color="black")
            plot_obj.set_axisbelow(True)
            
            plot_obj.set_xlabel("\nNumber of pairs of maps used", fontsize=16)
            plot_obj.set_ylabel("Variance ["+r"$uK^2$"+"]\n", fontsize=16)
            plot_obj.set_title("Variances of {} maps of {} sub-field\n".\
                               format(self.map_id, sub_field), fontsize=18)
            
            
            inset_plot = inset_axes(
                             plot_obj, width="30%", height="35%",
                             loc=3, borderpad=3)
            
            y_data = numpy.asarray(
                         self.stddevs_individs[sub_field].values())**2/vu
            y_data /= 2 * 1e3
            
            range_lo = numpy.percentile(y_data,  1)
            range_hi = numpy.percentile(y_data, 90)
            inset_plot.hist(y_data, bins=50, range=(range_lo, range_hi),
                            histtype="stepfilled", color="blue", alpha=0.1)
            
            inset_plot.tick_params(axis="both", which="both",
                                   direction="in", labelsize=12, left=False)
            inset_plot.set_yticklabels([])
            inset_plot.set_title("Blue Distribution  "+\
                                 "[1000 "+r"${\mu}K^2$"+"]",
                                 fontsize=13)
            
            
            more_file_name = "variances_of_signals_and_noises_in_" + sub_field
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        self.log("Done.")
        self.log("")
    
    
    def __call__(self, frame):
        
        if "Id" not in frame.keys():
            return
        
        if frame["Id"] != self.map_id:
            return
        
        if "CoaddedObservationIDs" in frame.keys():
            self.obs_id_list = frame["CoaddedObservationIDs"]
        
        try:
            k = "FluctuationMetricsIndividualSignalMapsTMapStandardDeviations"
            self.stddevs_individs = frame[k]
        except:
            pass
        
        if self.make_figs_for_noise_levels:
            if "NoiseLevelsFromCoaddedTMaps" in frame.keys():
                if not self.fig_an_made:
                    self.make_figure_for_all_noise_levels(frame)
                    self.fig_an_made = True
                if not self.fig_rn_made:
                    self.make_figures_for_recent_noise_levels(frame)
                    self.fig_rn_made = True
                if not self.fig_oa_made:
                    self.make_figures_for_order_of_addition(frame)
                    self.fig_oa_made = True
        
        if self.make_figs_for_fluc_metrics:
            k1 = "FluctuationMetricsCoaddedSignalMapsTMapStandardDeviations"
            k2 = "FluctuationMetricsCoaddedNoiseMapsTMapStandardDeviations"
            if (k1 in frame.keys())   and \
               (k2 in frame.keys())   and \
               (not self.fig_fl_made) and \
               (self.stddevs_individs is not None):
                self.make_figures_for_variances_of_maps(frame)
                self.fig_fl_made = True




class MakeFiguresForDistributionsOfMapRelatedQuantities(object):
    
    def __init__(self,
                 map_id=None, map_type=None,
                 fig_rc=False, fig_fl=False, fig_pt=False,
                 fig_rp=False, fig_ns=False,
                 figure_title_font_size=11,
                 directory_to_save_figures=".",
                 logging_function=logging.info):
        
        self.log = logging_function
        self.map_id   = map_id
        self.map_type = map_type
        
        self.make_fig_for_frac_cal_chng = fig_rc
        self.make_fig_for_fluc_metrics  = fig_fl
        self.make_fig_for_ptg_offsets   = fig_pt
        self.make_fig_for_xspec_ratios  = fig_rp
        self.make_fig_for_noise_levels  = fig_ns
        
        self.fig_rc_made = False
        self.fig_fl_made = False
        self.fig_pt_made = False
        self.fig_rp_made = False
        self.fig_ns_made = False
        
        self.ttl_fs  = figure_title_font_size
        font_sizes = determine_various_font_sizes(self.ttl_fs)
        self.lbl_fs = font_sizes[0]
        self.tck_fs = font_sizes[1]
        self.lgd_fs = 0.9*font_sizes[2]
        
        self.fields  = ["ra0hdec-44.75", "ra0hdec-52.25",
                        "ra0hdec-59.75", "ra0hdec-67.25",
                        "ra5hdec-24.5" , "ra5hdec-31.5",
                        "ra5hdec-38.5" , "ra5hdec-45.5",
                        "ra5hdec-52.5" , "ra5hdec-59.5"]
        self.el_dict = {"ra0hdec-44.75": "el 0"   , "ra0hdec-52.25": "el 1",
                        "ra0hdec-59.75": "el 2"   , "ra0hdec-67.25": "el 3",
                        "ra5hdec-24.5" : "el 0"   , "ra0hdec-31.5" : "el 1",
                        "ra5hdec-38.5" : "el 2"   , "ra5hdec-45.5" : "el 3",
                        "ra5hdec-52.5" : "el 4"   , "ra5hdec-59.5" : "el 5"}
        self.cl_dict = {"ra0hdec-44.75": "#1f77b4", "ra0hdec-52.25": "#ff7f0e",
                        "ra0hdec-59.75": "#2ca02c", "ra0hdec-67.25": "#d62728",
                        "ra5hdec-24.5" : "#1f77b4", "ra5hdec-31.5" : "#ff7f0e",
                        "ra5hdec-38.5" : "#2ca02c", "ra5hdec-45.5" : "#d62728",
                        "ra5hdec-52.5" : "#9467bd", "ra5hdec-59.5" : "#8c564b"}
        
        self.alpha   = 0.5
        self.ln_wdth = 1.25
        self.ln_styl = "solid"
        
        self.fig_dir = directory_to_save_figures
    
    
    def get_full_fig_title(self, additional_title, no_map_type=False):
        
        if self.map_id == "90GHz":
            map_id_for_fig = "95GHz"
        else:
            map_id_for_fig = self.map_id
        
        if no_map_type:
            prefix = map_id_for_fig + "  -  "
        else:
            prefix = map_id_for_fig + " " + self.map_type + " map" + "  -  "
        return prefix + additional_title
    
    
    def get_full_file_name(self, additional_name):
        
        prefix = self.map_id + "-" + self.map_type + "_" + "map" + "_"
        return prefix + additional_name + ".png"
    
    
    def determine_histogram_range(self, map_vector):
        
        pctl_los = []
        pctl_his = []
        for sub_field, data in map_vector.items():
            pctl_los.append(numpy.nanpercentile(data,  1))
            pctl_his.append(numpy.nanpercentile(data, 99))
        lo_of_lo = numpy.nanmin(pctl_los)
        hi_of_hi = numpy.nanmax(pctl_his)
        
        return lo_of_lo, hi_of_hi
    
    
    def draw_histograms(self, plot_obj, map_vector, xlim_left, xlim_right):
        
        xlim_right += 0.2 * (xlim_right - xlim_left)
        xlim_left  -= 0.2 * (xlim_right - xlim_left)
        max_occurs = []
        
        for sub_field, data in map_vector.items():
            n, bins, patches = \
                plot_obj.hist(data, bins=50, range=(xlim_left, xlim_right),
                              label=sub_field,
                              histtype="step",
                              linewidth=self.ln_wdth, linestyle=self.ln_styl,
                              color=self.cl_dict[sub_field], alpha=self.alpha)
            mode = (bins[numpy.argmax(n)] + bins[numpy.argmax(n)+1]) / 2.0
            plot_obj.axvline(mode, linestyle="dotted", linewidth=self.ln_wdth,
                             label="Mode = {:4.2e}".format(mode),
                             color=self.cl_dict[sub_field], alpha=self.alpha)
            max_occurs.append(numpy.max(n))
        
        max_of_max_n = numpy.max(max_occurs)
        plot_obj.set_ylim(bottom=-0.02*max_of_max_n, top=1.08*max_of_max_n)
        plot_obj.set_xlim(left=xlim_left, right=xlim_right)
                
        bin_size = bins[1] - bins[0]
        plot_obj.text(0.02, 0.97, s="bin size\n = {:4.2e}".format(bin_size),
                      horizontalalignment="left", verticalalignment="top",
                      transform=plot_obj.transAxes,
                      color="black", fontsize=self.lgd_fs)

        plot_obj.tick_params(axis="x", labelsize=self.tck_fs)
        plot_obj.tick_params(axis="y", left=False, labelleft=False)
        plot_obj.legend(loc="upper right",
                            fontsize=self.lgd_fs, framealpha=0.2)
    
    
    def indicate_histogram_statistics(self, plot_obj, map_vector):
        
        start_ht = 0.80
        counter  = 0
        n_fld    = len(map_vector.keys())
        for sub_field, data in map_vector.items():
            counter += 1
            med = numpy.nanmedian(data)
            mad = numpy.nanmedian(numpy.absolute(data-med))
            opkwargs = {"transform": plot_obj.transAxes,
                        "horizontalalignment": "left",
                        "verticalalignment": "top",
                        "color": self.cl_dict[sub_field],
                        "fontsize": self.lgd_fs*0.9}
            plot_obj.text(0.01, start_ht-counter*0.11,
                          s="MED = {:4.2e}".format(med),
                          **opkwargs)
            plot_obj.text(0.01, start_ht-(counter*0.11+0.05),
                          s="MAD = {:4.2e}".format(mad),
                          **opkwargs)
    
    
    def make_figure_for_distributions_of_cal_vs_el(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        original_frac_changes_data = \
            frame["MedianCalibratorResponseAllBolos"+\
                  "FractionalChangesTopToBottom"]
        data_reorganized = {}
        for sub_field, map_double in original_frac_changes_data.items():
            data_reorganized[sub_field] = \
                numpy.asarray(map_double.values()) * 100
        
        xlim_lefts  = {"90GHz": -3,  "150GHz": -3,  "220GHz": -2}
        xlim_rights = {"90GHz":  15, "150GHz":  15, "220GHz":  10}
        
        self.draw_histograms(
            plot_obj, data_reorganized,
            xlim_lefts[self.map_id], xlim_rights[self.map_id])
        
        self.indicate_histogram_statistics(plot_obj, data_reorganized)
        
        xlabel = "\n" + "Median percentage change in CalibratorResponse"
        ylabel = ""
        title  = self.get_full_fig_title(
                     "Distributions of the medians of the" + "\n" +\
                     "percentage changes in CalibratorResponse "  +\
                     "from bottom to top\n")
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, title, self.ttl_fs, add_legend=False)
        
        file_name = self.get_full_file_name(
                        "distributions_of_cal_response_changes")
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_distributions_of_mean_weights(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        original_mean_weights_data = \
            frame["FluctuationMetricsIndividualSignalMapsMeansOfTTWeights"]
        wu = 1 / (core.G3Units.mK * core.G3Units.mK)
        data_reorganized = {}
        for sub_field, map_double in original_mean_weights_data.items():
            data_reorganized[sub_field] = \
                numpy.asarray(map_double.values()) / wu
        
        season = get_season_based_on_fields(
                     original_mean_weights_data.keys())
        if season == "winter":
            xlim_lefts  = {"90GHz":  20.0, "150GHz":  30.0, "220GHz":  1.0}
            xlim_rights = {"90GHz": 120.0, "150GHz": 180.0, "220GHz": 15.0}
        elif season == "summer":
            xlim_lefts  = {"90GHz":  20.0, "150GHz":  30.0, "220GHz":  1.0}
            xlim_rights = {"90GHz": 160.0, "150GHz": 250.0, "220GHz": 25.0}
        
        self.draw_histograms(
            plot_obj, data_reorganized,
            xlim_lefts[self.map_id], xlim_rights[self.map_id])
        
        self.indicate_histogram_statistics(plot_obj, data_reorganized)
        
        xlabel = "\n" + "Mean TT weight  " + "[1 / " + r"${mK_{CMB}}^2$" + "]"
        ylabel = ""
        title  = self.get_full_fig_title(
                     "Distributions of the" + "\n" +\
                     "mean TT weights of individual maps\n")
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, title, self.ttl_fs, add_legend=False)
        
        file_name = self.get_full_file_name(
                        "distributions_of_individual_mean_tt_weights")
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_distributions_of_pntng_offsets(self, frame):
        
        for coordinate in ["Ra", "Dec"]:
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            original_offsets_data_sets = \
                [frame["Delta"+coordinate+"sOf"+\
                       rank+"BrightestSourceFromEachSubfield"] \
                 for rank in ["1st", "2nd", "3rd"]]
            
            data_reorganized = {}
            for sub_field, map_double in original_offsets_data_sets[0].items():
                three_offsets_all_obss = \
                    [[ds[sub_field][obs_id] \
                      for ds in original_offsets_data_sets] \
                     for obs_id in map_double.keys()]
                avgs_offsets = [numpy.nanmean(three_offsets_one_obs) \
                                for three_offsets_one_obs            \
                                in  three_offsets_all_obss]
                data_reorganized[sub_field] = \
                    numpy.asarray(avgs_offsets)/core.G3Units.arcsec
            
            self.draw_histograms(
                plot_obj, data_reorganized, -45, 45)
            
            self.indicate_histogram_statistics(plot_obj, data_reorganized)
            
            xlabel = "\n" + "Average " + coordinate + " offset [arcsec.]"
            ylabel = ""
            title  = self.get_full_fig_title(
                         "Distributions of the averages of the" + "\n" +\
                         coordinate+" offsets of the point sources\n")
            set_ax_labels_and_title(
                plot_obj, xlabel, ylabel, title, self.ttl_fs, add_legend=False)
            
            file_name = self.get_full_file_name(
                            "distribution_of_average_"+coordinate+"_offsets")
            save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_distributions_of_avg_xspc_ratios(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        original_noise_data = \
            frame["AveragesOfRatiosOfSPTxPlancktoPlckxPlckTTspectra"]
        data_reorganized = {}
        for sub_field, map_double in original_noise_data.items():
            data_reorganized[sub_field] = \
                numpy.asarray(map_double.values())
        
        season = get_season_based_on_fields(data_reorganized.keys())
        if season  == "winter":
            xlim_lefts  = {"90GHz": 0.5, "150GHz": 0.5, "220GHz": 0.2}
            xlim_rights = {"90GHz": 2.0, "150GHz": 2.0, "220GHz": 5.0}
        if season  == "summer":
            xlim_lefts  = {"90GHz": 0.5, "150GHz": 0.5, "220GHz": 0.2}
            xlim_rights = {"90GHz": 2.0, "150GHz": 2.0, "220GHz": 5.0}
        
        self.draw_histograms(
            plot_obj, data_reorganized,
            xlim_lefts[self.map_id], xlim_rights[self.map_id])
        
        self.indicate_histogram_statistics(plot_obj, data_reorganized)
        
        xlabel = "\n" + "Average of SPT x Planck / Planck x Planck"
        ylabel = ""
        title  = self.get_full_fig_title(
                     "Distributions of the averages of the ratios of\n"+\
                     "SPT x Planck to Planck x Planck spectra\n")
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, title, self.ttl_fs, add_legend=False)
        
        file_name = self.get_full_file_name(
                        "distributions_of_crass_spectra_ratios")
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_distributions_of_noise_levels(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        original_noise_data = frame["NoiseLevelsFromIndividualTMaps"]
        nu = core.G3Units.uK * core.G3Units.arcmin
        data_reorganized = {}
        for sub_field, map_double in original_noise_data.items():
            data_reorganized[sub_field] = \
                numpy.asarray(map_double.values())/nu
        
        season = get_season_based_on_fields(data_reorganized.keys())
        if season  == "winter":
            xlim_lefts  = {"90GHz": 100, "150GHz":  80, "220GHz": 280}
            xlim_rights = {"90GHz": 220, "150GHz": 180, "220GHz": 650}
        elif season == "summer":
            xlim_lefts  = {"90GHz":  60, "150GHz":  60, "220GHz": 200}
            xlim_rights = {"90GHz": 220, "150GHz": 220, "220GHz": 740}

        self.draw_histograms(
            plot_obj, data_reorganized,
            xlim_lefts[self.map_id], xlim_rights[self.map_id])
        
        self.indicate_histogram_statistics(plot_obj, data_reorganized)
        
        xlabel = "\n" + "Noise level  " + r"$[{\mu}K \cdot arcmin]$"
        ylabel = ""
        title  = self.get_full_fig_title(
                     "Distributions of the" + "\n" +\
                     "noise levels of individual maps\n")
        set_ax_labels_and_title(
            plot_obj, xlabel, ylabel, title, self.ttl_fs, add_legend=False)
        
        file_name = self.get_full_file_name(
                        "distributions_of_individual_noise_levels")
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def __call__(self, frame):
        
        if "Id" not in frame.keys():
            return
        
        if frame["Id"] != self.map_id:
            return
        
        if self.make_fig_for_frac_cal_chng:
            if not self.fig_rc_made:
                if "MedianCalibratorResponseAllBolos"+\
                   "FractionalChangesTopToBottom" in frame.keys():
                    self.log("")
                    self.log("Makding a figure for "
                             "distributions of responsivity change ...")
                    self.make_figure_for_distributions_of_cal_vs_el(frame)
                    self.fig_rc_made = True
                    self.log("Done.")
                    self.log("")
        
        if self.make_fig_for_fluc_metrics:
            if not self.fig_fl_made:
                if "FluctuationMetricsIndividualSignalMapsMeansOfTTWeights" \
                in frame.keys():
                    self.log("")
                    self.log("Making a figure for "
                             "distributions of mean TT weights ...")
                    self.make_figure_for_distributions_of_mean_weights(frame)
                    self.fig_fl_made = True
                    self.log("Done.")
                    self.log("")
        
        if self.make_fig_for_ptg_offsets:
            if not self.fig_pt_made:
                if "DeltaRasOf1stBrightestSourceFromEachSubfield" \
                in frame.keys():
                    self.log("")
                    self.log("Making a figure for "
                             "distributions of pointing offsets ...")
                    self.make_figure_for_distributions_of_pntng_offsets(frame)
                    self.fig_pt_made = True
                    self.log("Done.")
                    self.log("")
        
        if self.make_fig_for_xspec_ratios:
            if not self.fig_rp_made:
                if "AveragesOfRatiosOfSPTxPlancktoPlckxPlckTTspectra" \
                in frame.keys():
                    self.log("")
                    self.log("Making a figure for "
                             "distributions of avgs. of X spectra ratios ...")
                    self.make_figure_for_distributions_of_avg_xspc_ratios(frame)
                    self.fig_rp_made = True
                    self.log("Done.")
                    self.log("")
        
        if self.make_fig_for_noise_levels:
            if not self.fig_ns_made:
                if "NoiseLevelsFromIndividualTMaps" in frame.keys():
                    self.log("")
                    self.log("Making a figure for "
                             "distributions of noise levels ...")
                    self.make_figure_for_distributions_of_noise_levels(frame)
                    self.fig_ns_made = True
                    self.log("Done.")
                    self.log("")
        
        return




class RecordIDsUsedForPlotting(object):
    
    def __init__(self, logging_function=print):
        
        self.already_recorded = False
        self.log = logging_function
    
    
    def __call__(self, frame):
        
        if self.already_recorded:
            return
        
        if "CoaddedObservationIDs" in frame:
            self.log("")
            self.log("Recording what observation IDs "
                     "were used to make figures...")
            all_ids = core.G3VectorInt()
            for sub_field, obs_ids in frame["CoaddedObservationIDs"].items():
                for obs_id in obs_ids:
                    all_ids.append(obs_id)
            
            new_frame = core.G3Frame()
            new_frame["IDsUsedForMakingFiguresLastTime"] = all_ids
            self.log("Done.")
            
            self.already_recorded = True
            
            return [new_frame, frame]




# ==============================================================================




# ==============================================================================
# Define functions needed to start the pipeline
# ------------------------------------------------------------------------------


def run(input_files=[], decide_whether_to_make_figures_at_all=False,
        map_id=None, map_type=None, coadded_data=True,
        make_figure_for_field_map=False,
        rebin_map_before_plotting=False, new_map_resolution=None,
        smooth_map_with_gaussian=False, gaussian_fwhm=None,
        color_bar_upper_limit=None, color_bar_lower_limit=None,
        make_figure_for_entire_weight_map=False,
        make_figure_for_weight_map_cross_section=False,
        make_figures_showing_time_variations=False,
        make_figures_showing_time_evolution=False,
        make_figures_showing_distributions=False,
        make_figure_for_flagging_statistics=False,
        make_figure_for_pW_to_K_factors=False,
        make_figures_for_responsivity_changes=False,
        make_figures_for_fluctuation_metrics=False,
        make_figures_for_pointing_discrepancies=False,
        make_figure_for_ratios_of_power_spectra=False,
        make_figure_for_noise_levels=False,
        left_xlimit_for_time_variations=None,
        right_xlimit_for_time_variations=None,
        figure_title_font_size=11,
        directory_to_save_figures='.',
        simpler_file_names=False,
        logger_name="", log_file=None):
    
    
    # - Define global variables
    
    good_input_files = [input_file for input_file in input_files 
                        if os.path.isfile(input_file)]
    
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter('%(message)s')
    
    if log_file is None:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.DEBUG)
        stream_handler.setFormatter(log_format)
        logger.addHandler(stream_handler)
    else:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(log_format)
        logger.addHandler(file_handler)
        logging.captureWarnings(True)
    
    log = logger.info
    
    
    # - Time to start!
    
    log("")
    log("# ======================= #")
    log("#  Start making figures!  #")
    log("# ======================= #")
    log("")
    
    
    if len(good_input_files) == 0:
        log("Actually, no applicable input files!")
        log("")
        return
    else:
        log("")
        log("- These are the input files to be used:")
        log("")
        for input_file in good_input_files:
            log(input_file.split("/")[-1])
        log("")
        log("")
    
    
    if decide_whether_to_make_figures_at_all == True:
        
        bookkeeping_file = \
            "bookkeeping_for_plotting_" + \
            os.path.basename(good_input_files[-1]).replace(".gz", "")
        bookkeeping_file = \
            os.path.join(directory_to_save_figures, bookkeeping_file)
        
        if os.path.isfile(bookkeeping_file):
            try:
                previous_ids = list(list(core.G3File(bookkeeping_file))[0] \
                                        ["IDsUsedForMakingFiguresLastTime"])
                iterator  = core.G3File(good_input_files[0])
                frame     = iterator.next()
                newer_ids = []
                for sub_field, obs_ids \
                in  frame["CoaddedObservationIDs"].items():
                    for obs_id in obs_ids:
                        newer_ids.append(obs_id)
                common_ids = set(list(previous_ids)) & set(list(newer_ids))
                if len(common_ids) == len(newer_ids):
                    log("")
                    log("* The input files do not seem to contain any")
                    log("* new information, so the same set of figures")
                    log("* will not be generated again!")
                    log("")
                    return
            except:
                pass
    
    
    pipeline = core.G3Pipeline()
    
    pipeline.Add(core.G3Reader,
                 filename=good_input_files)
    
    pipeline.Add(MakeFiguresForFieldMapsAndWeightMaps,
                 fig_f=make_figure_for_field_map,
                 fig_w=make_figure_for_entire_weight_map,
                 fig_cr=make_figure_for_weight_map_cross_section,
                 rebin_map_before_plotting=rebin_map_before_plotting,
                 new_map_resolution=new_map_resolution,
                 smooth_map_with_gaussian=smooth_map_with_gaussian,
                 gaussian_fwhm=gaussian_fwhm,
                 map_type=map_type,
                 map_id=map_id,
                 coadded_data=coadded_data,
                 custom_vmin_field_map=color_bar_lower_limit,
                 custom_vmax_field_map=color_bar_upper_limit,
                 figure_title_font_size=figure_title_font_size,
                 simpler_file_names=simpler_file_names,
                 directory_to_save_figures=directory_to_save_figures,
                 logging_function=log)
    
    if make_figures_showing_time_variations:
        pipeline.Add(MakeFiguresForTimeVariationsOfMapRelatedQuantities,
                     map_id=map_id,
                     map_type=map_type,
                     fig_fs=make_figure_for_flagging_statistics,
                     fig_tc=make_figure_for_pW_to_K_factors,
                     fig_rc=make_figures_for_responsivity_changes,
                     fig_fl=make_figures_for_fluctuation_metrics,
                     fig_pt=make_figures_for_pointing_discrepancies,
                     fig_rp=make_figure_for_ratios_of_power_spectra,
                     fig_ns=make_figure_for_noise_levels,
                     xlim_left=left_xlimit_for_time_variations,
                     xlim_right=right_xlimit_for_time_variations,
                     figure_title_font_size=figure_title_font_size,
                     directory_to_save_figures=directory_to_save_figures,
                     logging_function=log)
    
    if make_figures_showing_time_evolution:
        pipeline.Add(MakeFiguresForTimeEvolutionOfMapRelatedQuantities,
                     map_id=map_id,
                     map_type=map_type,
                     fig_fl=make_figures_for_fluctuation_metrics,
                     fig_ns=make_figure_for_noise_levels,
                     figure_title_font_size=figure_title_font_size,
                     directory_to_save_figures=directory_to_save_figures,
                     logging_function=log)
    
    if make_figures_showing_distributions:
        pipeline.Add(MakeFiguresForDistributionsOfMapRelatedQuantities,
                     map_id=map_id,
                     map_type=map_type,
                     fig_rc=make_figures_for_responsivity_changes,
                     fig_fl=make_figures_for_fluctuation_metrics,
                     fig_pt=make_figures_for_pointing_discrepancies,
                     fig_rp=make_figure_for_ratios_of_power_spectra,
                     fig_ns=make_figure_for_noise_levels,
                     figure_title_font_size=figure_title_font_size,
                     directory_to_save_figures=directory_to_save_figures,
                     logging_function=log)
    
    if decide_whether_to_make_figures_at_all == True:
        pipeline.Add(RecordIDsUsedForPlotting,
                     logging_function=log)
        pipeline.Add(lambda frame: "IDsUsedForMakingFiguresLastTime" in frame)
        pipeline.Add(core.G3Writer,
                     filename=bookkeeping_file)
    
    if log_file is None:
        profile = True
    else:
        profile = False
    
    pipeline.Run(profile=profile)
    
    log("\n")


# ==============================================================================




# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                 description="This script is mainly intended to be used to "
                             "make figures showing coadded maps and relevant "
                             "quantities such as noise levels obtained by the "
                             "script fields_coadding.py.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                 epilog="Hopefully you will like the figures generated!")
    
    parser.add_argument("-I", "--input_files",
                        action="store", type=str, nargs="+", default=[],
                        help="G3 files that contain maps and map-related "
                             "analysis results.")
    
    parser.add_argument("-A", "--decide_whether_to_make_figures_at_all",
                        action="store_true", default=False,
                        help="Whether to check if it is actually necessary to "
                             "make any figures at all. If the input files do "
                             "not have new data (new maps added, new analysis "
                             "results) that were not present in them when "
                             "figures were made last time, then there is no "
                             "need to generate figures again.")
    
    parser.add_argument("-m", "--make_figure_for_field_map",
                        action="store_true", default=False,
                        help="Whether a figure showing the coadded field map "
                             "is to be made.")
    
    parser.add_argument("-i", "--map_id",
                        action="store", type=str, default="90GHz",
                        help="A string indicating figures for what type of "
                             "maps (mainly what band) are supposed to be made. "
                             "Each data frame in the input files need to have "
                             "an 'Id' key that stores a value that matches "
                             "this argument in order for the data in the "
                             "frame to be used to generate figures.")
    
    parser.add_argument("-y", "--map_type",
                        action="store", type=str, choices=["T"], default="T",
                        help="Whether figures of T, Q, or U maps are desired. "
                             "Currently, the script is still under development "
                             "and ignores options other than 'T'.")
    
    parser.add_argument("-c", "--coadded_data",
                        action="store_true", default=False,
                        help="Whether the frames contain maps resulting from "
                             "the coaddition of multiple observations' maps. "
                             "The script uses different titles and file names "
                             "for the figures depending on whether this option "
                             "is True or False. The default is False.")
    
    parser.add_argument("-u", "--color_bar_upper_limit",
                        action="store", type=float, default=None,
                        help="The upper limit of the values that will be used "
                             "in the colorbar for showing a field map, "
                             "if a value automatically determined by the "
                             "script is not desired.")
    
    parser.add_argument("-l", "--color_bar_lower_limit",
                        action="store", type=float, default=None,
                        help="The lower limit of the values that will be used "
                             "in the colorbar for showing a field map.")
    
    parser.add_argument("-b", "--rebin_map_before_plotting",
                        action="store_true", default=False,
                        help="Whether to rebin maps to a different resolution "
                             "before making figures for them.")
    
    parser.add_argument("-B", "--new_map_resolution",
                        action="store", type=float, default=None,
                        help="The new resolution that will be used during the "
                             "rebinning mentioned above. This number should be "
                             "expressed in the units of arcminute.")
    
    parser.add_argument("-t", "--smooth_map_with_gaussian",
                        action="store_true", default=False,
                        help="Whether to smooth a map by convolving it with "
                             "a 2D Gaussian.")
    
    parser.add_argument("-T", "--gaussian_fwhm",
                        action="store", type=float, default=None,
                        help="The full width at half maximum of the Gaussian "
                             "expressed in the units of arcminute.")
    
    parser.add_argument("-w", "--make_figure_for_entire_weight_map",
                        action="store_true", default=False,
                        help="Whether to make a figure for entire coadded TT "
                             "weight map.")
    
    parser.add_argument("-W", "--make_figure_for_weight_map_cross_section",
                        action="store_true", default=False,
                        help="Whether to make a figure only showing the cross "
                             "section of the weight map along the RA = 0h "
                             "contour.")
    
    parser.add_argument("-V", "--make_figures_showing_time_variations",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of several map-related quantities.")
    
    parser.add_argument("-E", "--make_figures_showing_time_evolution",
                        action="store_true", default=False,
                        help="Whether to make figures showing time evolution "
                             "of a few map-related quantities.")
    
    parser.add_argument("-D", "--make_figures_showing_distributions",
                        action="store_true", default=False,
                        help="Whether to make figures showing distributions "
                             "of several map-related quantities.")
    
    parser.add_argument("-F", "--make_figure_for_flagging_statistics",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of average numbers of flagged detectors "
                             "over time.")
    
    parser.add_argument("-C", "--make_figure_for_pW_to_K_factors",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of the temperature calibration factors pW/K or "
                             "their distributions.")
    
    parser.add_argument("-r", "--make_figures_for_responsivity_changes",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of the fractional changes in detectors' response "
                             "to the calibrator at the top of each sub-field "
                             "with respect to the bottom or their "
                             "distributions.")
    
    parser.add_argument("-U", "--make_figures_for_fluctuation_metrics",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations, "
                             "time evolution, and/or distributions of some "
                             "basic metrics of fluctuations of map values, "
                             "such as variances and mean weights.")
    
    parser.add_argument("-p", "--make_figures_for_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of pointing discrepancies and/or their "
                             "distributions.")
    
    parser.add_argument("-x", "--make_figure_for_ratios_of_power_spectra",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of the ratios of the SPT x Planck power spectra "
                             "to Planck x Planck ones and/or their "
                             "distributions.")
    
    parser.add_argument("-n", "--make_figure_for_noise_levels",
                        action="store_true", default=False,
                        help="Whether to make figures showing time variations "
                             "of noise levels from individual observations or "
                             "time evolution of noise levels of coadded maps, "
                             "depending on which type of data is available. "
                             "In the former case, figures showing "
                             "distributions noise levels will also be made.")
    
    parser.add_argument("-L", "--left_xlimit_for_time_variations",
                        action="store", type=int, default=None,
                        help="The observation ID that will be used as the "
                             "left limit of the x-axis of a figure that shows "
                             "time variations of certain quantity.")
    
    parser.add_argument("-R", "--right_xlimit_for_time_variations",
                        action="store", type=int, default=None,
                        help="The observation ID that will be used as the "
                             "right limit of the x-axis of a figure that shows "
                             "time variations of certain quantity.")
    
    parser.add_argument("-z", "--figure_title_font_size",
                        action="store", type=float, default=11,
                        help="The font size to be used for figure titles. "
                             "The font sizes of other elements such as "
                             "axis labels will be determined "
                             "based on this value.")
    
    parser.add_argument("-d", "--directory_to_save_figures",
                        action="store", type=str, default=".",
                        help="The directory where figures will be saved.")
    
    parser.add_argument("-f", "--simpler_file_names",
                        action="store_true", default=False,
                        help="Whether to omit infomration on sub-field, "
                             "observation ID, and so on in the names of the "
                             "PNG files that show field maps and weight maps.")
    
    parser.add_argument("-g", "--logger_name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "
                             "record log messages.")
    
    parser.add_argument("-G", "--log_file",
                        type=str, action="store", default=None,
                        help="The file to which the logger will send messages.")
    
    
    arguments = parser.parse_args()
    
    run(**vars(arguments))

    
# ==============================================================================

