## ******* Always need to import modules! ******* ##

import  matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["mathtext.default"] = "regular"
from matplotlib import pyplot

import  os
import  gc
import  sys
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




# ==============================================================================
# Define functions and pipeline modules needed for the pipeline
# ------------------------------------------------------------------------------


# Define functions related to general plotting utilities
# ------------------------------------------------------------------------------

def get_figure_and_plot_objects(w=12.0, h=7.5, dpi=None):
    
    figure_obj = pyplot.figure(figsize=(w, h), dpi=dpi)
    plot_obj   = figure_obj.add_subplot(1, 1, 1)
    return figure_obj, plot_obj


def determine_various_font_sizes(title_fontsize):

    labels_fs = title_fontsize - 1.0
    ticks_fs  = title_fontsize - 1.0
    legend_fs = title_fontsize - 2.0
    return labels_fs, ticks_fs, legend_fs


def set_lims(plot_obj, l, r, b, t):
    
    if l is not None: plot_obj.set_xlim(left=l)
    if r is not None: plot_obj.set_xlim(right=r)
    if b is not None: plot_obj.set_ylim(bottom=b)
    if t is not None: plot_obj.set_ylim(top=t)


def set_ticks(plot_obj, xta, xti, xtl, yta, yti, ytl, ttl_fs,
              xtrot="horizontal"):

    lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(ttl_fs)
    if xta is not None: plot_obj.set_xticks(xta, minor=False)
    if xti is not None: plot_obj.set_xticks(xti, minor=True)
    if yta is not None: plot_obj.set_yticks(yta, minor=False)
    if yti is not None: plot_obj.set_yticks(yti, minor=True)
    if xtl is not None: plot_obj.set_xticklabels(xtl, fontsize=tck_fs, rotation=xtrot)
    if ytl is not None: plot_obj.set_yticklabels(ytl, fontsize=tck_fs)
    plot_obj.grid(which="major", linestyle="dotted", linewidth=0.20)
    plot_obj.grid(which="minor", linestyle="dotted", linewidth=0.10)


def set_labels_and_title(plot_obj, xlabel, ylabel, title, ttl_fs):
    
    lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(ttl_fs)
    plot_obj.tick_params(axis="both", which="major",
                         direction="in", labelsize=tck_fs)
    plot_obj.tick_params(axis="both", which="minor",
                         direction="in")
    plot_obj.legend(loc="upper right", fontsize=lgd_fs, framealpha=0.2)
    plot_obj.set_xlabel(xlabel, fontsize=lbl_fs)
    plot_obj.set_ylabel(ylabel, fontsize=lbl_fs)
    plot_obj.set_title(title, fontsize=ttl_fs)


def save_figure_etc(figure_obj, fig_dir, file_name):
    
    figure_obj.tight_layout()
    fig_path = os.path.join(fig_dir, file_name)
    figure_obj.savefig(fig_path, bbox_inches="tight")
    pyplot.close(figure_obj)




# Define pipeline modules and more specific plotting functions
# ------------------------------------------------------------------------------


class PossiblyMakeFiguresForFieldMapsAndWeightMaps(object):
    
    def __init__(self,
                 fig_f=False, fig_w=False, fig_cr=False,
                 rebin_map_before_plotting=False, new_map_resolution=None,
                 smooth_map_with_gaussian=False, gaussian_fwhm=None,
                 map_type=None, map_id=None, coadded_data=True,
                 custom_vmin_field_map=None,
                 custom_vmax_field_map=None,
                 figure_title_font_size=11,
                 simpler_file_names=True,
                 directory_to_save_figures="",
                 logging_function=logging.info):
        
        self.make_figure_for_field_map    = fig_f
        self.make_figure_for_weight_map   = fig_w
        self.make_figure_for_w_map_cr_sec = fig_cr
        self.rebin_map_before_plotting    = rebin_map_before_plotting
        self.new_res  = new_map_resolution
        self.smooth_map_with_gaussian = smooth_map_with_gaussian
        self.gaussian_fwhm = gaussian_fwhm
        self.map_type = map_type
        self.map_id   = map_id
        self.coadded_data = coadded_data
        self.custom_vmin  = custom_vmin_field_map
        self.custom_vmax  = custom_vmax_field_map
        self.simpler_file_names = simpler_file_names
        self.ttl_fs = figure_title_font_size
        self.fig_dir = directory_to_save_figures
        self.log = logging_function
        self.obs_fr = None
        
        
    def rebin_map_to_diff_res(self, frame, new_res):
        
        map_parameters = {"x_len"       : 1200,
                          "y_len"       : 1200,
                          "res"         : new_res*core.G3Units.arcmin,
                          "proj"        : frame["T"].proj,
                          "alpha_center": frame["T"].alpha_center,
                          "delta_center": frame["T"].delta_center,
                          "pol_type"    : frame["T"].pol_type,
                          "coord_ref"   : frame["T"].coord_ref}
        rebin_factor = int(new_res/(frame["T"].res/core.G3Units.arcmin))
        print("The rebin factor is", rebin_factor)

        if self.map_type == "T":
            if "Wunpol" in frame.keys():
                new_field_map  = coordinateutils.FlatSkyMap(**map_parameters)
                new_weight_map = coordinateutils.FlatSkyMap(**map_parameters)
                coordinateutils.reproj_map(frame["T"],
                                           new_field_map,
                                           rebin_factor)
                coordinateutils.reproj_map(frame["Wunpol"].TT,
                                           new_weight_map,
                                           rebin_factor)
                new_frame = core.G3Frame(core.G3FrameType.Map)
                new_frame["T"]         = new_field_map
                new_frame["Wunpol"]    = core.G3SkyMapWeights()
                new_frame["Wunpol"].TT = new_weight_map
        
        return new_frame
    
    
    def get_field_map(self, frame, map_type):
        
        if map_type == "T":
            if "Wunpol" in frame:
                return frame["T"]/frame["Wunpol"].TT, "T map"
    
    
    def get_weight_map(self, frame, map_type):
        
        if map_type == "T":
            return frame["Wunpol"].TT, "TT weight map"
    
    
    def get_title_and_file_name_of_figure_for_map(
            self, obs_fr, obs_id_list, vals_for_stats, res, mp_ty_str):
        
        if (obs_fr is not None) and (obs_id_list is None):
            source  = obs_fr["SourceName"]
            obs_id  = str(obs_fr["ObservationID"])
            date    = str(obs_fr["ObservationStart"]).split(".")[0]
            resol   = str(res/core.G3Units.arcmin)+"' map"
            if self.smooth_map_with_gaussian:
                resol += " smoothed w/ " + \
                         str(self.gaussian_fwhm) + "' Gaussian"
            title_a = source + "  " + obs_id + " (" + date + ") " + \
                      "   " + self.map_id + " " + mp_ty_str #+ " (" + resol + ")"
        
        elif (obs_fr is None) and (obs_id_list is not None):
            source  = "1500 square degrees field"
            obs_ids = []
            for oids_one_field in obs_id_list.values():
                for oid in oids_one_field:
                    obs_ids.append(oid)
            min_id  = str(numpy.min(obs_ids))
            max_id  = str(numpy.max(obs_ids))
            min_dt  = str(std_processing.obsid_to_g3time(min_id)).split(":")[0]
            max_dt  = str(std_processing.obsid_to_g3time(max_id)).split(":")[0]
            obs_id  = min_id + " to " + max_id
            dt_rng  = "from " + min_dt + " to " + max_dt
            el_dic  = {"ra0hdec-44.75": "el 0", "ra0hdec-52.25": "el 1",
                       "ra0hdec-59.75": "el 2", "ra0hdec-67.25": "el 3"}
            n_obss  = ["{"+str(len(obss))+" "+el_dic[source]+"}" \
                       for source, obss in obs_id_list.items()]
            n_obss  = " ".join(n_obss)
            resol   = str(res/core.G3Units.arcmin)+"' map"
            if self.smooth_map_with_gaussian:
                resol += " smoothed w/ " + \
                         str(self.gaussian_fwhm) + "' Gaussian"
            title_a = source + "   " + self.map_id + " coadded " + mp_ty_str + \
                      "s " + "\n" + \
                      "(data from observations taken " + dt_rng + ":" + "\n" + \
                      n_obss + ", " + obs_id + ")"
        else:
            raise RuntimeError("Unclear how to build the fig. title!")
        
        if vals_for_stats is None:
            title_b = "(No information on statistics)"
        else:
            median  = str(numpy.round(numpy.median(vals_for_stats)))
            pctl_10 = str(numpy.round(numpy.percentile(vals_for_stats, 15), 3))
            pctl_90 = str(numpy.round(numpy.percentile(vals_for_stats, 85), 3))
            title_b = "(15th pctl. = " +pctl_10+",  "+\
                      "Median = " +median+ ",  " + "85th pctl. = " +pctl_90+ ")"
        
        full_ttl = title_a + "\n" + title_b + "\n"
        file_nm  = source + "-" + obs_id + self.map_id + "_" + mp_ty_str
        file_nm += (file_nm + ".png").replace(" ", "_")
        
        return full_ttl, file_nm
    
    
    def visualize_entire_map(
            self, map_to_view,
            w=12.0, h=8.0, dpi=None, aspect="equal",
            cmap="gray", custom_vmin=None, custom_vmax=None,
            cbar_label="", fig_title="",
            fig_dir="", file_name="map.png"):
        
        figure_obj, plot_obj = get_figure_and_plot_objects(w=w, h=h, dpi=100)
        
        if (custom_vmin == None) or (custom_vmax == None):
            pctl_lo = numpy.nanpercentile(numpy.asarray(map_to_view),  1)
            pctl_hi = numpy.nanpercentile(numpy.asarray(map_to_view), 99)
            larger  = numpy.max([numpy.abs(pctl_lo), numpy.abs(pctl_hi)])
            vmin = -1.0 * larger
            vmax =  1.0 * larger
            cbar_label += "  ( 1 to 99 percentile range )"
        else:
            vmin = custom_vmin
            vmax = custom_vmax
        
        cax = plot_obj.imshow(
                  numpy.asarray(map_to_view),
                  origin="lower", aspect=aspect, interpolation="none",
                  cmap=cmap, vmin=vmin, vmax=vmax)
        
        cbar = figure_obj.colorbar(
                   cax, ax=plot_obj, pad=0.01, shrink=0.75, aspect=30)
        cbar.ax.tick_params(labelsize=self.ttl_fs-2.0)
        
        plot_obj.tick_params(
            axis="both", which="both",
            bottom=False, top=False, labelbottom=False, labeltop=False,
            left=False, right=False, labelleft=False, labelright=False)
        
        lbl_fs, tck_fs, lgd_fs = determine_various_font_sizes(self.ttl_fs)
        cbar.ax.set_ylabel(cbar_label, fontsize=lbl_fs)
        plot_obj.set_title(fig_title, fontsize=self.ttl_fs)
                
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def visualize_normalized_map_cross_section(
            self, map_data,
            fig_w=12.0, fig_h=7.5,
            xlabel="", ylabel="", fig_title="",
            fig_dir="", file_name="cross_section.png"):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        vertical_cr_sec_values = \
            numpy.asarray(map_data).transpose()[map_data.shape[1]//2]
        vertical_cr_sec_values = \
            vertical_cr_sec_values / numpy.max(vertical_cr_sec_values)        
        
        plot_obj.plot(vertical_cr_sec_values, linewidth=0.45, color="black")
        
        deg  = core.G3Units.deg
        ra   = 0.0
        decs = numpy.linspace(-74.75, -37.25, 11)
        pids = [map_data.angle_to_pixel(ra*deg, dec*deg) for dec in decs]
        idcs = [pid//map_data.shape[1]-1 for pid in pids]
        
        xlim_left  = idcs[0]
        xlim_right = idcs[-1]
                
        set_lims(plot_obj, xlim_left, xlim_right, -0.02, 1.05)
        set_ticks(plot_obj, idcs[1:-1], None, [str(dec) for dec in decs[1:-1]],
                  None, None, None, self.ttl_fs)
        set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Observation:
            self.obs_fr = frame
        
        if frame.type == core.G3FrameType.Map:
            if self.map_id == frame["Id"]:
                
                self.log("")
                self.log("* Found a Map frame with ID "+self.map_id+"!")
                self.log("* Figures will possibly be made for this frame.")
                self.log("")
                
                if self.make_figure_for_field_map:                    
                    if self.rebin_map_before_plotting:
                        self.log("Rebinning the field map to "
                                 "%s arcminute resolution...", self.new_res)
                        new_frame = \
                            self.rebin_map_to_diff_res(frame, self.new_res)
                        fd_mp, fd_mp_str = \
                            self.get_field_map(new_frame, self.map_type)
                        self.log("Done.")
                        self.log("")
                    else:
                        fd_mp, fd_mp_str = \
                            self.get_field_map(frame, self.map_type)
                    
                    self.log("Making a figure for the field map...")
                    self.log("  Some basic properties of the map:")
                    self.log(" "*4 +"Map shape    : %s", fd_mp.shape)
                    self.log(" "*4 +"Map units    : %s", fd_mp.units)
                    self.log(" "*4 +"Map size     : %s", fd_mp.size)
                    self.log(" "*4 +"x resolution : %s arc minutes",
                             fd_mp.x_res/core.G3Units.arcmin)
                    self.log(" "*4 +"y resolution : %s arc minutes",
                             fd_mp.y_res/core.G3Units.arcmin)
                    
                    fd_mp   /= core.G3Units.mK
                    res      = fd_mp.x_res
                    cbar_lbl = "\n" + r"$mK_{CMB}$"
                    
                    fd_mp = numpy.asarray(fd_mp)
                    if self.smooth_map_with_gaussian:
                        sigma  = self.gaussian_fwhm * core.G3Units.arcmin
                        sigma /= res
                        sigma /= 2.0 * numpy.sqrt(2.0 * numpy.log(2.0))
                        self.log("  Smoothing the map with a gaussian "
                                 "whose sigma is %s...", sigma)
                        fd_mp = ndimage.gaussian_filter(fd_mp, sigma, mode="nearest")
                    
                    self.log("  Preparing the figure title...")
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                    else:
                        obs_id_list = None
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                            self.obs_fr, obs_id_list,
                            fd_mp[numpy.where(numpy.isfinite(fd_mp))],
                            res, fd_mp_str)
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                    fd_mp_str.replace(" ", "_") + ".png"
                    
                    self.log("  Actually making a figure...")
                    self.visualize_entire_map(
                        fd_mp, dpi=100,
                        custom_vmin=self.custom_vmin,
                        custom_vmax=self.custom_vmax,
                        cbar_label=cbar_lbl, fig_title=fig_ttl, file_name=fl_nm)
                    self.log("Done.")
                    self.log("")
                    del fd_mp
                    gc.collect()
                
                if self.make_figure_for_w_map_cr_sec:
                    wt_mp, wt_mp_str = self.get_weight_map(frame, self.map_type)
                    
                    res = wt_mp.x_res
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                    else:
                        obs_id_list = None
                    
                    xlabel  = "\n" + "Declination [degree]"
                    ylabel  = "Normalized Weight [unitless]" + "\n"
                    
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                                    self.obs_fr, obs_id_list,
                                    None, res, wt_mp_str)
                    
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                wt_mp_str.replace(" ", "_") + \
                                "_cross_sectional_view.png"
                    fig_ttl = self.map_id + " " + wt_mp_str + "   " + "\n" + \
                              "Cross sectional view " + \
                              "along the RA = 0h contour" + "\n"
                    
                    self.log("Making a figure for a cross section of "
                             "the weight map...")
                    self.visualize_normalized_map_cross_section(
                        wt_mp, xlabel=xlabel, ylabel=ylabel,
                        fig_title=fig_ttl, file_name=fl_nm)
                    self.log("Done.")
                    self.log("")
                    del wt_mp
                    gc.collect()

                if self.make_figure_for_weight_map:
                    if self.rebin_map_before_plotting:
                        self.log("Rebinning the weight map to "
                                 "%s arcminute resolution...", self.new_res)
                        new_frame = \
                            self.rebin_map_to_diff_res(frame, self.new_res)
                        wt_mp, wt_mp_str = \
                            self.get_weight_map(frame, self.map_type)
                        self.log("Done.")
                        self.log("")
                    else:
                        wt_mp, wt_mp_str = \
                            self.get_weight_map(frame, self.map_type)
                    
                    self.log("Making a figure for the entire weight map...")
                    
                    res = wt_mp.x_res
                    
                    wt_mp  = numpy.asarray(wt_mp)
                    wt_mp /= 1.0 / (core.G3Units.mK*core.G3Units.mK)
                    cbar_lbl = "\n" + "1 / " + r"${mK_{CMB}}^2$"
                    
                    self.log("  Preparing the figure title...")
                    if self.coadded_data:
                        obs_id_list = frame["CoaddedObservationIDs"]
                    else:
                        obs_id_list = None
                    non_zero_wghts = wt_mp[numpy.where(wt_mp!=0.0)]
    
                    fig_ttl, fl_nm = \
                        self.get_title_and_file_name_of_figure_for_map(
                            self.obs_fr, obs_id_list,
                            non_zero_wghts, res, wt_mp_str)
                    if self.simpler_file_names:
                        fl_nm = self.map_id + "-" + \
                                    wt_mp_str.replace(" ", "_") + ".png"
                    
                    wt_mp[numpy.where(wt_mp==0.0)] = numpy.nan
                    
                    self.log("  Actually making a figure...")
                    self.visualize_entire_map(
                        wt_mp, dpi=100,
                        custom_vmin=0.0,
                        custom_vmax=numpy.nanmax(wt_mp),
                        cbar_label=cbar_lbl, fig_title=fig_ttl, file_name=fl_nm)
                    
                    if self.simpler_file_names:
                        file_name = self.map_id + "-" + \
                                    wt_mp_str.replace(" ", "_") + ".png"
                    self.log("Done.")
                    self.log("")
                    del wt_mp
                    gc.collect()
                




class PossiblyMakeFiguresForTimeVariationsOfMapRelatedQuantities(object):
    
    def __init__(self,
                 fig_fs=False, fig_tc=False, fig_fl=False,
                 fig_pt=False, fig_ns=False,
                 map_type=None, map_id=None,
                 xlim_left=None, xlim_right=None,
                 figure_title_font_size=11,
                 directory_to_save_figures="",
                 logging_function=logging.info):
        
        self.make_fig_for_flggg_stats   = fig_fs
        self.make_fig_for_temp_cal_facs = fig_tc
        self.make_fig_for_fluc_metrics  = fig_fl
        self.make_fig_for_ptg_offsets   = fig_pt
        self.make_fig_for_noise_levels  = fig_ns
        self.map_type = map_type
        self.map_id   = map_id
        
        self.already_made_figures = False
                
        self.xlim_left  = xlim_left
        self.xlim_right = xlim_right
        self.ttl_fs     = figure_title_font_size
        self.el_dict = {"ra0hdec-44.75": "el 0"   , "ra0hdec-52.25": "el 1",
                        "ra0hdec-59.75": "el 2"   , "ra0hdec-67.25": "el 3"}
        self.cl_dict = {"ra0hdec-44.75": "#1f77b4", "ra0hdec-52.25": "#ff7f0e",
                        "ra0hdec-59.75": "#2ca02c", "ra0hdec-67.25": "#d62728"}
        self.ln_wdth = 0.25
        self.ln_styl = "dotted"
        self.mrkrsz  = 8.0
        
        self.fig_title_prefix = map_id + " " + map_type + " " + "map" + "  -  "
        self.file_name_prefix = map_id + "-" + map_type + "_" + "map" + "_"
        
        self.fig_dir = directory_to_save_figures
        self.log = logging_function
    
    
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
                month = tstr.split("-")[1]
                date  = tstr.split("-")[0]
                hour  = tstr.split(":")[1]
                if no_hour:
                    label = month + " " + date
                else:
                    label = month + " " + date + "\n" + hour + "h"
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
        
        return xtick_locs_major, xtick_labels_major, xtick_locs_minor
    
    
    def indicate_out_of_range_values(
            self, plot_obj, x_values, y_values, ylims_dict, color):
        
        ytop = ylims_dict["top"]
        ybot = ylims_dict["bottom"]
        near_ytop = ybot + 0.90 * (ytop - ybot)
        near_ybot = ybot + 0.15 * (ytop - ybot)
        note_fs = self.ttl_fs-4
        note_lw = 1.5*self.ln_wdth
        note_st = "--"
        
        for i, x_value in enumerate(x_values):
            y_value = y_values[i]
            if not numpy.isfinite(y_value):
                plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)
                plot_obj.text(x_value, near_ytop, "NaN", rotation="vertical",
                              color=color, fontsize=note_fs)
            elif y_value >= ytop:
                plot_obj.scatter(x_value, ytop, color=color, marker=6)
                plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)
                """plot_obj.text(x_value, near_ytop, "O\nu\nt\nl\ni\ne\nr",
                                 color=color, fontsize=note_fs)"""
            elif y_value <= ybot:
                plot_obj.scatter(x_value, ybot, color=color, marker=7)
                plot_obj.axvline(x_value, color=color,
                                 linewidth=note_lw, linestyle=note_st)
                """plot_obj.text(x_value, near_ybot, "O\nu\nt\nl\ni\ne\nr",
                                 color=color, fontsize=note_fs)"""
                              
    
    def get_full_fig_title(self, additional_title):
        
        return self.fig_title_prefix + additional_title
    
    
    def get_full_file_name(self, additional_name):
        
        return self.file_name_prefix + additional_name + ".png"
    
    
    def make_figure_for_flagging_statistics(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        cl_mk_lbl_dict  = \
            {"BadHk"                  : ["#1f77b4" , "Bad HK"      , "." , 8],
             "Latched"                : ["#1f77b4" , "Latched"     , "p" , 5],
             "Overbiased"             : ["#1f77b4" , "Saturated"   , "s" , 5],
             "Oscillating"            : ["#9467bd" , "Oscillating" , "v" , 5],
             "Glitchy"                : ["#9467bd" , "Glitchy"     , "^" , 6],
             "UnphysicalLowVariance"  : ["#9467bd" , "Low Var."    , "D" , 4],
             "BadCalSn"               : ["#1fbecf" , "Low Cal S/N" , "1" , 8],
             "MissingFluxCalibration" : ["#1fbecf" , "No Fluxcal"  , "2" , 8],
             "PostCalibrationNaNs"    : ["#1fbecf" , "Has NaNs"    , "3" , 8],
             "BadWeight"              : ["#1fbecf" , "Bad Weight"  , "4" , 8],
             "Others"                 : ["black"   , "Others"      , "." , 6],
             "TotalNotFlagged"        : ["green"   , "Survived"    , "*" , 8],
             "TotalRemoved"           : ["red"     , "Removed"     , "x" , 6]}
        
        for flag_type in cl_mk_lbl_dict.keys():
            key_name = "FlaggingStatisticsAverageNumberOf" + flag_type
            averages = frame[key_name]
            
            times_and_averages = []
            for sub_field, oids_and_avgs in averages.items():
                for oid, avg in oids_and_avgs.items():
                    times_and_averages.append((int(oid), avg))
            times_and_averages = sorted(times_and_averages, key=itemgetter(0))
            times    = [entry[0] for entry in times_and_averages]
            averages = [entry[1] for entry in times_and_averages]
            
            plot_obj.plot(
                times, averages, label=cl_mk_lbl_dict[flag_type][1],
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=cl_mk_lbl_dict[flag_type][2],
                markersize=cl_mk_lbl_dict[flag_type][3],
                color=cl_mk_lbl_dict[flag_type][0], alpha=0.8)
        
        plot_obj.set_yscale("symlog", linthreshy=150,
                            subsy=[2, 3, 4, 5, 6, 7, 8, 9])
        
        ylims_dict = {"bottom": -10, "top": 5000}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, 6.2)
        
        set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                           ylims_dict["bottom"], ylims_dict["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs)
        
        xlabel = "\n" + "Time (UTC)"
        ylabel = "Averages (over all scans of an obs.)" + "\n"
        more_title = "Average number of flagged detectors" + "\n"
        fig_title  = self.get_full_fig_title(more_title)
        
        set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "average_numbers_of_flagged_detectors"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_observation_durations(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        obs_dur_data = frame["ObservationDurations"]
        
        for sub_field, oids_and_olens in obs_dur_data.items():
            oids = [int(oid) for oid in oids_and_olens.keys()]
            lens = [oids_and_olens[str(oid)]/core.G3Units.min \
                    for oid in oids_and_olens.keys()]
            label = self.el_dict[sub_field]
            color = self.cl_dict[sub_field]

            plot_obj.plot(
                oids, lens, label=label,
                linestyle="None", marker=".", markersize=self.mrkrsz,
                color=color, alpha=0.8)
        
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, 3.5)
        ylims_dict = {"bottom": 0.0, "top": 140.0}
        
        set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                           ylims_dict["bottom"], ylims_dict["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs)
        
        xlabel = "\n" + "Time (UTC)"
        ylabel = "Duration [minute]" + "\n"
        fig_title = "Durations of observations" + "\n"
        
        set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        file_name = "observation_durations.png"
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_calibration_factors(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        waf_cl_dict = {"W172": "#1f77b4", "W174": "#ff7f0e",
                       "W176": "#2ca02c", "W177": "#d62728",
                       "W180": "#9467bd", "W181": "#8c564b",
                       "W188": "#e377c2", "W203": "#7f7f7f",
                       "W204": "#bcbd22", "W206": "#17becf"}
        cal_uni = core.G3Units.pW / core.G3Units.K
        
        ylims_dict = { "90GHz": {"bottom": -0.14, "top": -0.03},
                      "150GHz": {"bottom": -0.21, "top": -0.02},
                      "220GHz": {"bottom": -0.04, "top": -0.00}}
        xlims_dict = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, 4.0)

        for wafer, color in waf_cl_dict.items():
            cal_data = frame["MeansOfTemperatureCalibrationFactors" + wafer]
            
            times_and_cals = []
            for sub_field, oids_and_cals in cal_data.items():
                for oid, cal in oids_and_cals.items():
                    times_and_cals.append((int(oid), cal/cal_uni))
            
            times_and_cals = sorted(times_and_cals, key=itemgetter(0))
            times = [entry[0] for entry in times_and_cals]
            cals  = [entry[1] for entry in times_and_cals]
            
            plot_obj.plot(
                times, cals, label=wafer,
                linestyle=self.ln_styl, linewidth=self.ln_wdth,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=0.8)
            
            self.indicate_out_of_range_values(
                        plot_obj, times, cals, ylims_dict[self.map_id], color)
        
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        xtick_locs_major, xtick_labels, xtick_locs_minor = \
            self.get_xticks_and_labels_from_obs_id_range(
                plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs)

        xlabel = "\n" + "Time (UTC)"
        ylabel = "Mean from each wafer [pW/K]" + "\n"
        more_title = "Temperature calibration factors" + "\n"
        fig_title  = self.get_full_fig_title(more_title)
        
        set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        more_file_name = "mean_temperature_calibration_factors"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_map_fluctuation_metrics(self, frame):
        
        key_prefix = "BasicMapFluctuationMetricsIndividualMaps"
        param_lists = \
            [{"key"   : key_prefix+ "MapStdDevs",
              "ylims" : { "90GHz": {"bottom": 0.15, "top": 0.60},
                         "150GHz": {"bottom": 0.15, "top": 0.60},
                         "220GHz": {"bottom": 0.75, "top": 1.65}},
              "yunits": core.G3Units.mK,
              "ylabel": "Standard deviation  " + r"$[mK_{CMB}]$" + "\n",
              "title" : "Standard deviation of map values" + "\n",
              "file"  : "standard_deviation_of_map_values"},
             {"key"   : key_prefix + "MeansOfWeights",
              "ylims" : { "90GHz": {"bottom": 20, "top": 100},
                         "150GHz": {"bottom": 10, "top": 150},
                         "220GHz": {"bottom":  0, "top":  12}},
              "yunits": 1 / (core.G3Units.mK * core.G3Units.mK),
              "ylabel": "Mean weight  " + \
                        "[1 / " + r"${mK_{CMB}}^2$" + "]" + "\n",
              "title" : "Mean of TT weight map" + "\n",
              "file"  : "mean_of_tt_weight_map_values"}]
        
        
        for plist in param_lists:
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            data = frame[plist["key"]]
            
            ylims = plist["ylims"]
            xlims = self.get_xlims_from_obs_id_range(
                         self.xlim_left, self.xlim_right, 3.5)
            
            for sub_field, xs_and_ys in data.items():
                xvals = sorted([int(oid) for oid in xs_and_ys.keys()])
                yvals = [xs_and_ys[str(oid)]/plist["yunits"] for oid in xvals]
                label = self.el_dict[sub_field]
                color = self.cl_dict[sub_field]
                
                plot_obj.plot(
                    xvals, yvals, label=label,
                    linestyle=self.ln_styl, linewidth=self.ln_wdth*2,
                    marker=".", markersize=self.mrkrsz,
                    color=color, alpha=0.8)
                
                self.indicate_out_of_range_values(
                    plot_obj, xvals, yvals, ylims[self.map_id], color)
            
            set_lims(plot_obj, xlims["left"], xlims["right"],
                     ylims[self.map_id]["bottom"], ylims[self.map_id]["top"])
            
            xtick_locs_major, xtick_labels, xtick_locs_minor = \
                self.get_xticks_and_labels_from_obs_id_range(
                    plot_obj, self.xlim_left, self.xlim_right)
            
            set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                      None, None, None, self.ttl_fs)
            
            xlabel = "\n" + "Time (UTC)"
            ylabel = plist["ylabel"]
            more_title = plist["title"]
            fig_title  = self.get_full_fig_title(more_title)
            
            set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
            
            more_file_name = plist["file"]
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_all_noise_levels(self, frame):
        
        figure_obj, plot_obj = get_figure_and_plot_objects()
        
        obs_data = frame["CoaddedObservationIDs"]
        
        if "NoiseFromIndividualMaps" in frame.keys():
            noise_from_running_coadds = False
            noise_data = frame["NoiseFromIndividualMaps"]
        elif "NoiseFromCoaddedMaps" in frame.keys():
            noise_from_running_coadds = True
            noise_data = frame["NoiseFromCoaddedMaps"]
            try:
                oprts_hist   = frame["NoiseCalculationsOperationsDoneToMaps"]
                ops_recorded = True
            except KeyError:
                ops_recorded = False
        else:
            raise NameError("Not clear where noise data are!")
        
        if noise_from_running_coadds:
            plot_obj.set_xscale("log")
            plot_obj.set_yscale("log")
            ylims_dict = { "90GHz": {"bottom":  5, "top": 250},
                          "150GHz": {"bottom":  5, "top": 250},
                          "220GHz": {"bottom": 15, "top": 700}}
            max_n_obss = 0
            for sub_field, noise_dict in noise_data.items():
                n_obss = len(noise_dict.keys())
                if n_obss > max_n_obss:
                    max_n_obss = n_obss
            xlims_dict = {"left" : 0.9,
                          "right": 10**numpy.ceil(numpy.log10(max_n_obss))}
        else:
            ylims_dict = { "90GHz": {"bottom":  90, "top": 210},
                          "150GHz": {"bottom":  90, "top": 210},
                          "220GHz": {"bottom": 300, "top": 700}}
            xlims_dict = self.get_xlims_from_obs_id_range(
                             self.xlim_left, self.xlim_right, 3.5)
        
        if noise_from_running_coadds:
            n_excluded = {}
        
        for sub_field, obs_ids in obs_data.items():
            if not noise_from_running_coadds:
                obs_ids = sorted(obs_ids)
            noise_units  = core.G3Units.uK * core.G3Units.arcmin
            noise_levels = [noise_data[sub_field][str(obs_id)] / noise_units \
                            for obs_id in obs_ids \
                            if str(obs_id) in noise_data[sub_field].keys()]
            
            label = self.el_dict[sub_field]
            color = self.cl_dict[sub_field]
            
            if noise_from_running_coadds:
                if ops_recorded:
                    operations_performed = \
                        numpy.array([oprts_hist[sub_field][str(obs_id)] \
                                     for obs_id in obs_ids])
                    valid_indices = numpy.where(operations_performed==-1.0)[0]
                    bad_indices   = numpy.where(operations_performed== 0.0)[0]
                    n_excluded[sub_field] = len(bad_indices)
                    x_data = range(1, len(valid_indices)+1)
                    noise_levels = numpy.asarray(noise_levels)[valid_indices]
                else:
                    if sub_field in frame["IgnoredObservationIDs"].keys():
                        n_bad_obs = frame["IgnoredObservationIDs"][sub_field]
                    else:
                        n_bad_obs = 0
                    n_excluded[sub_field] = n_bad_obs
                    x_data = [obs_id for obs_id in obs_ids \
                              if str(obs_id) in noise_data[sub_field].keys()]
                    x_data = range(1, len(x_data)+1)
            else:
                x_data = obs_ids
            plot_obj.plot(
                x_data, noise_levels, label=label,
                linestyle=self.ln_styl, linewidth=self.ln_wdth*2,
                marker=".", markersize=self.mrkrsz,
                color=color, alpha=0.5)
            
            if noise_from_running_coadds:
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
                              color=color, alpha=0.8, linewidth=self.ln_wdth)
                explanation = "Slope from LLS fit: " + str(power)[0:6]
                ylocs_dict  = {"ra0hdec-44.75": 0.19,
                               "ra0hdec-52.25": 0.14,
                               "ra0hdec-59.75": 0.09,
                               "ra0hdec-67.25": 0.04}
                plot_obj.text(0.03, ylocs_dict[sub_field], explanation,
                              transform=plot_obj.transAxes,
                              horizontalalignment="left",
                              color=color, fontsize=self.ttl_fs-4)
            else:
                self.indicate_out_of_range_values(plot_obj, x_data, noise_levels,
                                                  ylims_dict[self.map_id], color)
                
        set_lims(plot_obj, xlims_dict["left"], xlims_dict["right"],
                 ylims_dict[self.map_id]["bottom"],
                 ylims_dict[self.map_id]["top"])
        
        
        if noise_from_running_coadds:
            xtick_locs_major = None
            xtick_labels     = None
            xtick_locs_minor = None
        else:
            xtick_locs_major, xtick_labels, xtick_locs_minor = \
                self.get_xticks_and_labels_from_obs_id_range(
                    plot_obj, self.xlim_left, self.xlim_right)
        
        set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                  None, None, None, self.ttl_fs)

        if noise_from_running_coadds:
            xlabel = "\n" + "Number of pairs of difference map added"
        else:
            xlabel = "\n" + "Time (UTC)"
        ylabel = "Average  " + r"$\sqrt{C_l}$" + "  " + "in the " + "\n" + \
                 "ell range [3000, 5000]  " + r"$[{\mu}K \cdot arcmin]$" + "\n"
        
        if noise_from_running_coadds:
            more_title = "Noise of year-to-date coadded noise maps\n"
            more_info  = []
            for sub_field in sorted(n_excluded.keys()):
                alias = self.el_dict[sub_field]
                n_bad = n_excluded[sub_field]
                more_info.append("{" + str(n_bad) + " " + alias + "}")
            more_info   = "(numbers of excluded maps: " + \
                          " ".join(more_info) + ")\n"
            more_title += more_info
        else:
            more_title = "Noise (+some signal) of individual maps\n"
        fig_title = self.get_full_fig_title(more_title)
        
        set_labels_and_title(plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
        
        if noise_from_running_coadds:
            more_file_name = "noise_levels_from_running_coadds"
        else:
            more_file_name = "noise_levels_from_individual_maps"
        file_name = self.get_full_file_name(more_file_name)
        
        save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_recent_noise_levels(self, frame):
        
        noise_key = "NoiseFromCoaddedMaps"
        if noise_key not in frame.keys():
            return
        
        all_noise_calculated = frame[noise_key]
        all_obs_ids_used     = frame["CoaddedObservationIDs"]
        try:
            all_operations_performed = frame["NoiseCalculationsOperationsDoneToMaps"]
            ops_recorded = True
        except KeyError:
            ops_recorded = False
        
        for sub_field, obs_ids_this_sf in all_obs_ids_used.items():
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            if ops_recorded:
                operations_this_sf = \
                    numpy.array([all_operations_performed[sub_field][str(obs_id)] \
                                 for obs_id in obs_ids_this_sf])
                noise_this_sf = \
                    numpy.array([all_noise_calculated[sub_field][str(obs_id)] \
                                 for obs_id in obs_ids_this_sf])
                
                ok_indices   = numpy.where(operations_this_sf!=0.0)[0]
                ok_obs_ids   = numpy.asarray(obs_ids_this_sf)[ok_indices][-50:]
                ok_noises    = numpy.asarray(noise_this_sf)[ok_indices][-50:]
                ok_noises   /= core.G3Units.uK * core.G3Units.arcmin
                dummy_x_data = numpy.arange(len(ok_obs_ids)) + 1
            else:
                l_ids = []
                r_ids = []
                for obs_id in frame["CoaddedObservationIDs"][sub_field]:
                    if str(obs_id) in frame["NoiseFromCoaddedMaps"][sub_field].keys():
                        r_ids.append(obs_id)
                    else:
                        l_ids.append(obs_id)
                ok_obs_ids = []
                for i in range(numpy.min([len(l_ids), len(r_ids)])):
                    ok_obs_ids.append(l_ids[i])
                    ok_obs_ids.append(r_ids[i])
                ok_noises = []
                for obs_id in ok_obs_ids:
                    if str(obs_id) in frame["NoiseFromCoaddedMaps"][sub_field].keys():
                        noise  = frame["NoiseFromCoaddedMaps"][sub_field][str(obs_id)]
                        noise /= core.G3Units.uK * core.G3Units.arcmin
                        ok_noises.append(noise)
                    else:
                        ok_noises.append(numpy.nan)
                ok_obs_ids = ok_obs_ids[-50:]
                ok_noises  = ok_noises[-50:]
                dummy_x_data = numpy.arange(len(ok_obs_ids)) + 1
            
            color = self.cl_dict[sub_field]
            
            plot_obj.plot(dummy_x_data, ok_noises,
                          linestyle=self.ln_styl, linewidth=self.ln_wdth*2,
                          marker=".", markersize=self.mrkrsz,
                          color=color)
            
            max_noise = numpy.nanmax(ok_noises)
            min_noise = numpy.nanmin(ok_noises)
            diff      = max_noise - min_noise
            ylim_top  = max_noise + 0.05 * diff
            near_top  = max_noise + 0.02 * diff
            
            set_lims(plot_obj, 0, dummy_x_data[-1]+1, None, ylim_top)
            
            for index, noise in enumerate(ok_noises):
                if not numpy.isfinite(noise):
                    plot_obj.text(dummy_x_data[index], near_top,
                                  "N/A", color=color, fontsize=7,
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
                            for obs_id in ok_obs_ids]
            
            set_ticks(plot_obj, dummy_x_data, [], xtick_labels,
                      None, None, None, 8, xtrot="vertical")
            
            xlabel = ""
            ylabel = "Noise " + r"$[{\mu}K \cdot arcmin]$"
            more_title = "Noise due to 50 recently added maps of " + sub_field 
            fig_title  = self.get_full_fig_title(more_title)
            
            set_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
            
            more_file_name = "noise_levels_from_observations_of_{}".\
                             format(sub_field)
            file_name = self.get_full_file_name(more_file_name)
            save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_order_of_addition(self, frame):
        
        noise_key = "NoiseFromCoaddedMaps"
        
        if noise_key not in frame.keys():
            return
        if "NoiseCalculationsOperationsDoneToMaps" in frame.keys():
            easy_stuff = True
        else:
            easy_stuff = False
        
        all_obs_ids_used = frame["CoaddedObservationIDs"]
        
        for sub_field, obs_ids_this_sf in all_obs_ids_used.items():
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            if not easy_stuff:
                l_ids = []
                r_ids = []
                for obs_id in obs_ids_this_sf:
                    if str(obs_id) in frame["NoiseFromCoaddedMaps"][sub_field].keys():
                        r_ids.append(obs_id)
                    else:
                        l_ids.append(obs_id)
                obs_ids_this_sf = []
                for i in range(numpy.min([len(l_ids), len(r_ids)])):
                    obs_ids_this_sf.append(l_ids[i])
                    obs_ids_this_sf.append(r_ids[i])
            
            plot_obj.plot(obs_ids_this_sf,
                          linestyle="None",
                          marker=".", markersize=self.mrkrsz,
                          color=self.cl_dict[sub_field])
            
            min_oid = numpy.min(obs_ids_this_sf)
            max_oid = numpy.max(obs_ids_this_sf)
            
            ylims_dict = self.get_xlims_from_obs_id_range(
                min_oid, max_oid, 0.0)
            ylims_dict["bottom"] = ylims_dict.pop("left")
            ylims_dict["top"]    = ylims_dict.pop("right")
            
            set_lims(plot_obj, 0, len(obs_ids_this_sf)+1,
                               ylims_dict["bottom"], ylims_dict["top"])
            
            ytick_locs_major, ytick_labels, ytick_locs_minor = \
                self.get_xticks_and_labels_from_obs_id_range(
                    plot_obj, min_oid, max_oid, no_hour=True)
            
            set_ticks(plot_obj, None, None, None,
                      ytick_locs_major, ytick_locs_minor, ytick_labels,
                      self.ttl_fs)
            
            xlabel = "\n" + "Number of observations added"
            ylabel = "Date of observation" + "\n"
            
            more_title = "Order in which observations of " + \
                         sub_field + " were added" + "\n"
            fig_title = self.get_full_fig_title(more_title)
            
            set_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
            
            more_file_name = "order_of_addition_of_maps_of_" + sub_field
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def make_figure_for_pointing_discrepancies(self, frame):
        
        marker_dict = {"1": ["*", 8], "2": ["p", 6], "3": [".", 9]}
        

        for coordinate in ["Ra", "Dec"]:
            figure_obj, plot_obj = get_figure_and_plot_objects()
            
            xlims_dict = self.get_xlims_from_obs_id_range(
                             self.xlim_left, self.xlim_right, 5.4)
            ylims_dict = {"bottom": -45, "top": 45}            

            for source_rank in ["1", "2", "3"]:
                key  = "Delta" + coordinate + "s" + "FromSources" + source_rank
                data = frame[key]
                
                for sub_field, ids_and_offsets in data.items():
                    obs_ids = sorted(ids_and_offsets.keys())
                    offsets = [ids_and_offsets[oid]/core.G3Units.arcsec \
                               for oid in obs_ids]
                    obs_ids = [int(oid) for oid in obs_ids]
                    
                    plot_obj.plot(
                        obs_ids, offsets,
                        label=self.el_dict[sub_field]+" src "+source_rank,
                        linestyle="", linewidth=self.ln_wdth,
                        marker=marker_dict[source_rank][0],
                        markersize=marker_dict[source_rank][1],
                        color=self.cl_dict[sub_field], alpha=0.8)
                    
                    self.indicate_out_of_range_values(
                        plot_obj, obs_ids, offsets, ylims_dict,
                        self.cl_dict[sub_field])
                        
            set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                               ylims_dict["bottom"], ylims_dict["top"])
            
            xtick_locs_major, xtick_labels, xtick_locs_minor = \
                self.get_xticks_and_labels_from_obs_id_range(
                    plot_obj, self.xlim_left, self.xlim_right)
            
            set_ticks(plot_obj, xtick_locs_major, xtick_locs_minor, xtick_labels,
                      None, None, None, self.ttl_fs)
            
            xlabel = "\n" + "Time (UTC)"
            if coordinate == "Ra":
                ylabel = "(Measured R.A. - True R.A.) " + \
                         r"$\times$" + " cos(True Dec.)"+ ' ["]' + "\n"
            else:
                ylabel = "Measured Dec. - True Dec." + ' ["]' + "\n"
            more_title = "Difference between measured and true " + \
                         coordinate + " of point sources" + "\n"
            
            fig_title = self.get_full_fig_title(more_title)
            
            set_labels_and_title(
                plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
            
            more_file_name = "delta_" + coordinate + "s_from_point_sources"
            file_name = self.get_full_file_name(more_file_name)
            
            save_figure_etc(figure_obj, self.fig_dir, file_name)
        
        
        for sub_field in self.el_dict.keys():
            xlims_dict = self.get_xlims_from_obs_id_range(
                             self.xlim_left, self.xlim_right, 6.4)
            ylims_dict = {"bottom": -45, "top": 45}     
            srccl_dict = {"1": "red", "2": "blue", "3": "green"}

            for coordinate in ["Ra", "Dec"]:
                figure_obj, plot_obj = get_figure_and_plot_objects()
                
                all_offsets = {"1": None, "2": None, "3": None}
                
                for source_rank in ["1", "2", "3"]:
                    key  = "Delta"+coordinate+"s"+"FromSources"+source_rank
                    data = frame[key][sub_field]
                    
                    obs_ids = sorted(data.keys())
                    offsets = numpy.array([data[oid]/core.G3Units.arcsec \
                                           for oid in obs_ids])
                    obs_ids = [int(oid) for oid in obs_ids]
                    all_offsets[source_rank] = offsets
                
                avgs = numpy.nanmean([all_offsets["1"],
                                      all_offsets["2"],
                                      all_offsets["3"]], axis=0)
                
                for source_rank, offsets in all_offsets.items():
                    plot_obj.plot(
                        obs_ids, offsets,
                        label="Src. "+source_rank+" offsets",
                        linestyle="None",
                        marker=marker_dict[source_rank][0],
                        markersize=marker_dict[source_rank][1],
                        color=srccl_dict[source_rank], alpha=0.8)
                    
                    self.indicate_out_of_range_values(
                        plot_obj, obs_ids, offsets, ylims_dict,
                        self.cl_dict[sub_field])
                
                plot_obj.plot(
                    obs_ids, avgs, label="Avg. offsets",
                    linestyle=self.ln_styl, linewidth=self.ln_wdth,
                    color="black", alpha=0.8)
                
                set_lims(plot_obj, xlims_dict["left"],   xlims_dict["right"],
                                   ylims_dict["bottom"], ylims_dict["top"])
                
                text_tr = 0.70
                separat = 0.05
                counter = 0
                for source_rank, offsets in all_offsets.items():
                    avg = numpy.nanmean(offsets)
                    plot_obj.text(0.98, text_tr-counter*separat,
                                  "Src. {} - avg.: {:+05.1f}".format(source_rank, avg),
                                  transform=plot_obj.transAxes,
                                  horizontalalignment="right",
                                  verticalalignment="top",
                                  color=srccl_dict[source_rank],
                                  fontsize=self.ttl_fs-4)
                    std = numpy.nanstd(offsets)
                    plot_obj.text(0.98, text_tr-(counter+4)*separat,
                                  "Src. {} - std.: {:04.1f}".format(source_rank, std),
                                  transform=plot_obj.transAxes,
                                  horizontalalignment="right",
                                  verticalalignment="top",
                                  color=srccl_dict[source_rank],
                                  fontsize=self.ttl_fs-4)
                    counter += 1
                
                avg_avgs = numpy.nanmean(avgs)
                std_avgs = numpy.nanstd(avgs)
                plot_obj.text(0.98, text_tr-8*separat,
                              "Avg. - avg.: {:+05.1f}".format(avg_avgs),
                              transform=plot_obj.transAxes,
                              horizontalalignment="right",
                              verticalalignment="top",
                              color="black",
                              fontsize=self.ttl_fs-4)
                plot_obj.text(0.98, text_tr-9*separat,
                              "Avg. - std.: {:04.1f}".format(std_avgs),
                              transform=plot_obj.transAxes,
                              horizontalalignment="right",
                              verticalalignment="top",
                              color="black",
                              fontsize=self.ttl_fs-4)                
                
                xtick_locs_major, xtick_labels, xtick_locs_minor = \
                    self.get_xticks_and_labels_from_obs_id_range(
                        plot_obj, self.xlim_left, self.xlim_right)
                
                set_ticks(plot_obj,
                          xtick_locs_major, xtick_locs_minor, xtick_labels,
                          None, None, None, self.ttl_fs)
                
                xlabel = "\n" + "Time (UTC)"
                if coordinate == "Ra":
                    ylabel = "(Measured R.A. - True R.A.) " + \
                             r"$\times$" + " cos(True Dec.)"+ ' ["]' + "\n"
                else:
                    ylabel = "Measured Dec. - True Dec." + ' ["]' + "\n"
                more_title = coordinate + " offsets of point sources " + \
                             "in " + sub_field + "\n"
                
                fig_title = self.get_full_fig_title(more_title)
                
                set_labels_and_title(
                    plot_obj, xlabel, ylabel, fig_title, self.ttl_fs)
                
                more_file_name = "delta_" + coordinate + \
                                 "s_from_point_sources_in_" + sub_field
                file_name = self.get_full_file_name(more_file_name)
                
                save_figure_etc(figure_obj, self.fig_dir, file_name)
    
    
    def __call__(self, frame):
        
        if self.already_made_figures:
            return
        
        elif (frame.type  == core.G3FrameType.Map) and \
             (self.map_id == frame["Id"]):
            
            if self.make_fig_for_flggg_stats:
                self.log("Making a figure for flagging statistics...")
                self.make_figure_for_flagging_statistics(frame)
                self.make_figure_for_observation_durations(frame)
                self.log("Done.\n")
            
            if self.make_fig_for_temp_cal_facs:
                self.log("Making a figure for calibration factors...")
                self.make_figure_for_calibration_factors(frame)
                self.log("Done.\n")
            
            if self.make_fig_for_fluc_metrics:
                self.log("Making figures for map fluctuation metrics...")
                self.make_figure_for_map_fluctuation_metrics(frame)
                self.log("Done.\n")
            
            if self.make_fig_for_noise_levels:
                self.log("Making figures for noise levels...")
                self.make_figure_for_all_noise_levels(frame)
                self.make_figure_for_recent_noise_levels(frame)
                self.make_figure_for_order_of_addition(frame)
                self.log("Done.\n")
            
            if self.make_fig_for_ptg_offsets:
                self.log("Making figures for pointing discrepancies...")
                self.make_figure_for_pointing_discrepancies(frame)
                self.log("Done.\n")
            
            self.already_made_figures = True





class RecordIDsUsedForPlotting(object):
    
    def __init__(self, logging_function=logging.info):
        
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


def run(input_files, decide_whether_to_make_figures_at_all=False,
        make_figures_for_field_maps=False, map_id="90GHz", map_type='T',
        coadded_data=False, color_bar_upper_limit=None,
        color_bar_lower_limit=None, make_figures_for_entire_weight_maps=False,
        make_figures_for_weight_maps_cross_section=False,
        rebin_map_before_plotting=False,
        new_map_resolution=None,
        smooth_map_with_gaussian=False,
        gaussian_fwhm=None,
        make_figures_for_flagging_statistics=False,
        make_figures_for_calibration_factors=False,
        make_figures_for_fluctuation_metrics=False,
        make_figures_for_pointing_discrepancies=False,
        make_figures_for_noise_levels=False,
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


    for input_file in good_input_files:

        log("-"*80)
        log("Taking actions on %s ...", input_file)
        log("-"*80)

        if decide_whether_to_make_figures_at_all == True:

            bookkeeping_file = "bookkeeping_for_plotting_" + \
                               os.path.basename(input_file)
            bookkeeping_file = os.path.join(directory_to_save_figures,
                                            bookkeeping_file)

            if os.path.isfile(bookkeeping_file):
                try:
                    previous_ids = list(list(core.G3File(bookkeeping_file))[0] \
                                        ["IDsUsedForMakingFiguresLastTime"])
                except:
                    previous_ids = []
                iterator  = core.G3File(input_file)
                frame     = iterator.next()
                newer_ids = []
                for sub_field, obs_ids in frame["CoaddedObservationIDs"].items():
                    for obs_id in obs_ids:
                        newer_ids.append(obs_id)

                common_ids = set(list(previous_ids)) & set(list(newer_ids))
                if len(common_ids) == len(newer_ids):
                    log("")
                    log("* This input file does not seem to contain any")
                    log("* new information, so the same set of figures")
                    log("* will not be generated again!")
                    log("")
                    continue



        pipeline = core.G3Pipeline()

        pipeline.Add(core.G3Reader,
                     filename=input_file)

        pipeline.Add(PossiblyMakeFiguresForFieldMapsAndWeightMaps,
                     fig_f=make_figures_for_field_maps,
                     fig_w=make_figures_for_entire_weight_maps,
                     fig_cr=make_figures_for_weight_maps_cross_section,
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

        pipeline.Add(PossiblyMakeFiguresForTimeVariationsOfMapRelatedQuantities,
                     fig_fs=make_figures_for_flagging_statistics,
                     fig_tc=make_figures_for_calibration_factors,
                     fig_fl=make_figures_for_fluctuation_metrics,
                     fig_pt=make_figures_for_pointing_discrepancies,
                     fig_ns=make_figures_for_noise_levels,
                     map_id=map_id,
                     map_type=map_type,
                     xlim_left=left_xlimit_for_time_variations,
                     xlim_right=right_xlimit_for_time_variations,
                     figure_title_font_size=figure_title_font_size,
                     directory_to_save_figures=directory_to_save_figures,
                     logging_function=log)

        pipeline.Add(RecordIDsUsedForPlotting,
                     logging_function=log)
        pipeline.Add(lambda frame: "IDsUsedForMakingFiguresLastTime" in frame)
        
        if decide_whether_to_make_figures_at_all == True:
            pipeline.Add(core.G3Writer,
                     filename=bookkeeping_file)
        
        if log_file is None:
            profile = True
        else:
            profile = False

        pipeline.Run(profile=profile)

    log("")

    
# ==============================================================================




# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                 description="This script makes figures for maps made from "
                             "field observations. It can handle both a map from "
                             "a single field observation and a map from "
                             "the coaddition of mutltple observations.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                 epilog="Hopefully you will like the figures generated!")

    parser.add_argument("input_files",
                        action="store", type=str, nargs="+",
                        help="G3 files that contain Map frames. The data stored "
                             "in these frames will be plotted in the figures.")

    parser.add_argument("-D", "--decide_whether_to_make_figures_at_all",
                        action="store_true", default=False,
                        help="Whether to check if it is actually necessary to "
                             "make any figures at all. If figures of coadded maps "
                             "are to be made, if the corresponding g3 file that "
                             "contains the data also records what observations "
                             "were used to make the coadded maps, and if there is "
                             "another g3 file that records what observations were "
                             "contained in the coadded maps when figures were "
                             "made last time, then this check can be made.")

    parser.add_argument("-m", "--make_figures_for_field_maps",
                        action="store_true", default=False,
                        help="Whether figures showing field maps are to be made.")

    parser.add_argument("-i", "--map_id",
                        action="store", type=str, default="90GHz",
                        help="The ID a map frame needs to have in order for "
                             "a figure of the maps  stored to be made.")

    parser.add_argument("-y", "--map_type",
                        action="store", type=str, choices=["T"], default="T",
                        help="Whether figures of T, Q, or U maps are desired. "
                             "Currently, the script is still under development "
                             "and ignores options other than 'T'.")

    parser.add_argument("-c", "--coadded_data",
                        action="store_true", default=False,
                        help="Whether the frames contain maps resulting from "
                             "the coaddition of multiple observations. "
                             "The script uses different titles and file names "
                             "for the figures depending on whether this option "
                             "is True or False. The default is False.")

    parser.add_argument("-u", "--color_bar_upper_limit",
                        action="store", type=float, default=None,
                        help="The upper limit of the values that will be used "
                             "in the colorbar for showing a field map.")

    parser.add_argument("-l", "--color_bar_lower_limit",
                        action="store", type=float, default=None,
                        help="The lower limit of the values that will be used "
                             "in the colorbar for showing a field map.")

    parser.add_argument("-w", "--make_figures_for_entire_weight_maps",
                        action="store_true", default=False,
                        help="Whether to make figures for entire weight maps, "
                             "i.e. color maps showing weights at every location "
                             "of the field.")
    
    parser.add_argument("-W", "--make_figures_for_weight_maps_cross_section",
                        action="store_true", default=False,
                        help="Whether to make figures only showing the cross "
                             "section of a weight map along certain RA contour.")
    
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
                             "a Gaussian.")
    
    parser.add_argument("-T", "--gaussian_fwhm",
                        action="store", type=float, default=None,
                        help="The full width at half maximum of the Gaussian "
                             "expressed in the units of arcminute.")

    parser.add_argument("-F", "--make_figures_for_flagging_statistics",
                        action="store_true", default=False,
                        help="Whether to make figures showing average number "
                             "of flagged detectors over time.")
    
    parser.add_argument("-C", "--make_figures_for_calibration_factors",
                        action="store_true", default=False,
                        help="Whether to make figures showing temperature "
                             "calibration factors over time.")
    
    parser.add_argument("-U", "--make_figures_for_fluctuation_metrics",
                        action="store_true", default=False,
                        help="Whether to make figures showing some basic "
                             "metrics of fluctuations of map values.")

    parser.add_argument("-p", "--make_figures_for_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to make figures showing pointing "
                             "discrepancies over time.")

    parser.add_argument("-n", "--make_figures_for_noise_levels",
                        action="store_true", default=False,
                        help="Whether to make figures showing noise levels "
                             "over time.")

    parser.add_argument("-L", "--left_xlimit_for_time_variations",
                        action="store", type=int, default=None,
                        help="The observation ID that will be used as the "
                             "lower limit of the x-axis of a figure that shows "
                             "time variations of certain quantity.")

    parser.add_argument("-R", "--right_xlimit_for_time_variations",
                        action="store", type=int, default=None,
                        help="The observation ID that will be used as the "
                             "upper limit of the x-axis of a figure that shows "
                             "time variations of certain quantity.")

    parser.add_argument("-z", "--figure_title_font_size",
                        action="store", type=float, default=11,
                        help="The font size to be used for figure titles. "
                             "Then, the font sizes of other things such as "
                             "axis labels will be determined "
                             "based on this value.")

    parser.add_argument("-d", "--directory_to_save_figures",
                        action="store", type=str, default=".",
                        help="The directory where generated figures "
                             "will be saved.")

    parser.add_argument("-f", "--simpler_file_names",
                        action="store_true", default=False,
                        help="Whether to omit infomration on source, obs. ID, "
                             "and so on in the names of the PNG files that "
                             "show field maps and weight maps.")
    
    parser.add_argument("-g", "--logger_name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "+\
                             "record log messages.")
    
    parser.add_argument("-G", "--log_file",
                        type=str, action="store", default=None,
                        help="The file to which the logger will send messages.")


    arguments = parser.parse_args()

    run(**vars(arguments))

    
# ==============================================================================

