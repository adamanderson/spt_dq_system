## ******* Always need to import modules! ******* ##

import  matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

import  os
import  gc
import  sys
import  glob
import  numpy 
import  scipy
import  argparse
from    operator    import  itemgetter
from    spt3g       import  core
from    spt3g       import  mapmaker
from    spt3g       import  std_processing




# ==============================================================================
# Process command line arguments and define global variables
# ------------------------------------------------------------------------------


# Process command line arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
             description="This script makes figures for maps made from "+\
                         "field observations. It can handle both a map from "+\
                         "a single field observation and a map from "+\
                         "the coaddition of mutltple observations.",
             epilog="Hopefully you will like the figures generated!")

parser.add_argument("input_files",
                    type=str, action="store", nargs="+",
                    help="G3 files that contain Map frames. The data stored "+\
                         "in these frames will be plotted in the figures.")

parser.add_argument("-c", "--coadded_data",
                    action="store_true", default=False,
                    help="Whether the frames contain maps resulting from "+\
                         "the coaddition of multiple observations. "+\
                         "The script uses different titles and file names "+\
                         "for the figures depending on whether this option "+\
                         "is True or False. The default is False.")

parser.add_argument("-y", "--map_type",
                    type=str, action="store", default="T",
                    help="Whether figures of T, Q, or U maps are desired. "+\
                         "The default is T. Currently, the script is still "+\
                         "under development and ignores options other than "+\
                         "'T'.")

parser.add_argument("-U", "--map_values_upper_limit",
                    type=float, action="store", default=None,
                    help="The upper limit of the values that will be used "+\
                         "in the colorbar for showing a map.")

parser.add_argument("-L", "--map_values_lower_limit",
                    type=float, action="store", default=None,
                    help="The lower limit of the values that will be used "+\
                         "in the colorbar for showing a map.")

parser.add_argument("-i", "--map_ids",
                    type=str, action="store", nargs="+", default=["90GHz"],
                    help="What ID a map needs to have in order for "+\
                         "a figure to be made.")

parser.add_argument("-w", "--make_figures_for_weights",
                    action="store_true", default=False,
                    help="Whether figures of weight maps are desired. "+\
                         "The default is False.")

parser.add_argument("-F", "--make_figures_for_flagging_statistics",
                    action="store_true", default=False,
                    help="Whether figures showing average number of flagged "+\
                         "detectors over time are desired. The default "+\
                         "is False.")

parser.add_argument("-p", "--make_figures_for_pointing_discrepancies",
                    action="store_true", default=False,
                    help="Whether figures showing pointing discrepancies "+\
                         "over time are desired. The default is False.")

parser.add_argument("-n", "--make_figures_for_noise_levels",
                    action="store_true", default=False,
                    help="Whether figures showing noise levels over time "+\
                         "are desired. The default is False.")

parser.add_argument("-N", "--integrated_noise",
                    action="store_true", default=False,
                    help="Whether the noise data correspond to integrated "+\
                         "noise over time. The default is False.")

parser.add_argument("-l", "--left_xlimit_for_figures",
                    type=int, action="store", default=None,
                    help="The observation ID that will be used as the "+\
                         "lower limit of the x-axis of a figure that shows "+\
                         "time variations of certain quantity such as noise.")

parser.add_argument("-r", "--right_xlimit_for_figures",
                    type=int, action="store", default=None,
                    help="The observation ID that will be used as the "+\
                         "upper limit of the x-axis of a figure that shows "+\
                         "time variations of certain quantity such as noise.")

parser.add_argument("-z", "--figure_title_font_size",
                    type=float, action="store", default=11,
                    help="The font size to be used for figure titles. "+\
                         "The font size of axis labels, etc. will be "+\
                         "determined based on this value.")

parser.add_argument("-j", "--make_figures_for_one_frame_only",
                    action="store_true", default=False,
                    help="Whether to make figures for one Map frame only, "+\
                         "no matter how many frames the input files have. "+\
                         "This option can be used when one just wants to "+\
                         "have an idea of what the figures look like. The "+\
                         "default is False.")

parser.add_argument("-s", "--only_show_figures",
                    action="store_true", default=False,
                    help="Whether to only show figures and not save them. "+\
                         "The default is False.")

parser.add_argument("-d", "--directory_to_save_figures",
                    type=str, action="store", default=".",
                    help="The directory where generated figures will be "+\
                         "saved. The default is '.'.")

parser.add_argument("-f", "--simpler_file_names",
                    action="store_true", default=False,
                    help="Whether to omit infomration on source, obs. ID, "+\
                         "and so on in the names of the PNG files.")

arguments = parser.parse_args()


print()
print("# ======================= #")
print("#  Start making figures!  #")
print("# ======================= #")
print()

print()
print("These are the input files:")
print()
for input_file in arguments.input_files:
    print(input_file.split("/")[-1])
print()
print()


if arguments.map_type == "T":
    def get_sky_map(frame):
        return frame["T"]
    def get_weight_map(frame):
        return frame["Wunpol"].TT
    map_type_str    = "T map"
    weight_type_str = "TT weight map"
else:
    print("Well, the argument 'map_type' is unrecognized,")
    print("so, the script will do nothing...")
    sys.exit()


if arguments.directory_to_save_figures[-1] == "/":
    arguments.directory_to_save_figures = \
        arguments.directory_to_save_figures[0:-1]



# Define global variables
# -----------------------------------------------------------------------------

valid_map_ids = arguments.map_ids
ttl_fs = arguments.figure_title_font_size
lbl_fs = ttl_fs - 1.0
tck_fs = ttl_fs - 1.0
lgd_fs = ttl_fs - 2.0
fig_wd = 12
fig_ht = 7.5


# ==============================================================================





# ==============================================================================
# Define plotting function
# ------------------------------------------------------------------------------

def visualize_a_field_map(
        twod_array, x_res, y_res, fig_w, fig_h, aspect,
        cmap, vmin=None, vmax=None,
        cbar_label="", title="",
        no_xy_labels=True, image_xlabel="", image_ylabel="",
        save_fig=False, filename=None):
    
    figure_for_2d_array = pyplot.figure(figsize=(fig_w, fig_h))
    image_of_2d_array   = figure_for_2d_array.add_subplot(1, 1, 1)
    
    if (vmin == None) or (vmax == None):
        cax = image_of_2d_array.imshow(numpy.asarray(twod_array),
                                       origin="lower", aspect=aspect,
                                       interpolation="none",
                                       cmap=cmap,)
    else:
        cax = image_of_2d_array.imshow(numpy.asarray(twod_array),
                                       origin="lower", aspect=aspect,
                                       interpolation="none",
                                       cmap=cmap, vmin=vmin, vmax=vmax)
    
    cbar = figure_for_2d_array.colorbar(
               cax, ax=image_of_2d_array, pad=0.01, shrink=0.9, aspect=30)
    
    cbar.ax.set_ylabel(cbar_label, fontsize=lbl_fs)
    
    if not no_xy_labels:
        image_of_2d_array.set_xticklabels(
            numpy.round(image_of_2d_array.get_xticks()*x_res/core.G3Units.deg),
            fontsize=lbl_fs)
        image_of_2d_array.set_yticklabels(
            numpy.round(image_of_2d_array.get_yticks()*y_res/core.G3Units.deg),
            fontsize=lbl_fs)
        image_of_2d_array.set_xlabel(image_xlabel, fontsize=lbl_fs)
        image_of_2d_array.set_ylabel(image_ylabel, fontsize=lbl_fs)
    else:
        image_of_2d_array.tick_params(axis="x", which="both",
                                      bottom=False, top=False, 
                                      labelbottom=False, labeltop=False)
        image_of_2d_array.tick_params(axis="y", which="both",
                                      left=False, right=False, 
                                      labelleft=False, labelright=False)
    
    image_of_2d_array.set_title(title, fontsize=ttl_fs)
    
    figure_for_2d_array.tight_layout()
    
    if not save_fig:
        pyplot.show(figure_for_2d_array)
    else:
        pyplot.savefig(filename, bbox_inches="tight")
    pyplot.close(figure_for_2d_array)




def make_histogram_and_cross_section_plot_of_a_field_map(
        twod_array, is_a_weight_map=True, also_histogram=False,
        hist_xlabel="", hist_title="",
        cr_sec_xticks=[], cr_sec_xticklabels=[],
        cr_sec_xlabel="", cr_sec_title="", fig_title="",
        save_fig=False, filename=""):
    
    figure_for_the_array = pyplot.figure(figsize=(12, 9))
    if also_histogram:
        plot_for_histogram   = figure_for_the_array.add_subplot(1, 2, 1)
        plot_for_cr_sec_view = figure_for_the_array.add_subplot(1, 2, 2)
    else:
        plot_for_cr_sec_view = figure_for_the_array.add_subplot(1, 1, 1)
    
    if also_histogram:
        if is_a_weight_map:
            values_for_histogram = numpy.asarray(twod_array)\
                                   [numpy.where(numpy.asarray(twod_array)!=0.0)]
            values_for_histogram = \
                values_for_histogram / numpy.max(values_for_histogram)
        else:
            values_for_histogram = numpy.asarray(twod_array).flatten()
        
        plot_for_histogram.hist(
            values_for_histogram,
            bins=200, histtype="step", color="black")
        plot_for_histogram.tick_params(
            axis="y", which="both",
            left=False, labelleft=False)
        plot_for_histogram.tick_params(
            axis="x", labelsize=tck_fs)
        plot_for_histogram.set_xlabel("\n"+hist_xlabel, fontsize=lbl_fs)
        plot_for_histogram.set_title(hist_title, fontsize=ttl_fs)
    
    
    vertical_cr_sec_values = \
        numpy.asarray(twod_array).transpose()[twod_array.shape[-1]//2]
    
    if is_a_weight_map:
        vertical_cr_sec_values = \
            vertical_cr_sec_values / numpy.max(vertical_cr_sec_values)
    
    ind_diff_for_n_deg = (cr_sec_xticks[0] - cr_sec_xticks[1]) * \
                         1 / (cr_sec_xticklabels[0] - cr_sec_xticklabels[1]) * \
                         3.0
    xlim_left  = numpy.max([0, cr_sec_xticks[-1]-int(ind_diff_for_n_deg)])
    xlim_right = cr_sec_xticks[0] + int(ind_diff_for_n_deg)
    
    cr_sec_xticklabels = [str(label) for label in cr_sec_xticklabels]
    
    plot_for_cr_sec_view.plot(
        vertical_cr_sec_values, linewidth=0.3,
        color="black")
    plot_for_cr_sec_view.tick_params(
        axis="both", labelsize=tck_fs)
    
    plot_for_cr_sec_view.set_xlim(left=xlim_left, right=xlim_right)
    plot_for_cr_sec_view.set_ylim(top=1.05, bottom=-0.025)
    
    plot_for_cr_sec_view.grid(linestyle="dotted", linewidth=0.3)
    
    plot_for_cr_sec_view.set_xticks(cr_sec_xticks)
    plot_for_cr_sec_view.set_xticklabels(cr_sec_xticklabels, fontsize=tck_fs)
    
    plot_for_cr_sec_view.set_xlabel("\n"+cr_sec_xlabel, fontsize=lbl_fs)
    plot_for_cr_sec_view.set_ylabel(hist_xlabel+"\n", fontsize=lbl_fs)
    plot_for_cr_sec_view.set_title(cr_sec_title, fontsize=ttl_fs)
    
    if also_histogram:
        figure_for_the_array.suptitle(fig_title, fontsize=ttl_fs, y=0.99)
    
    figure_for_the_array.tight_layout(pad=3.0, w_pad=2.0)
    
    if not save_fig:
        pyplot.show(figure_for_the_array)
    else:
        figure_for_the_array.savefig(filename, bbox_inches="tight")
    pyplot.close(figure_for_the_array)




def decide_on_xticks_from_obs_id_range(min_obs_id, max_obs_id):
    
    ticks  = []
    labels = []
    
    total_secs   = max_obs_id - min_obs_id
    current_tick = min_obs_id
    
    while current_tick <= (max_obs_id + 5):
        tstr  = str(std_processing.obsid_to_g3time(current_tick))
        month = tstr.split("-")[1]
        date  = tstr.split("-")[0]
        hour  = tstr.split(":")[1]
        label = month + " " + date + "\n" + hour + "h"
        ticks.append(current_tick)
        labels.append(label)
        if total_secs < (4 * 24 * 3600 + 10):
            current_tick += 12 * 3600
        elif total_secs < (8 * 24 * 3600):
            current_tick += 24 * 3600
        else:
            current_tick += 7 * 24 * 3600
    
    return ticks, labels




def make_figure_for_noise_levels(
        noise_data, is_integrated_noise=False,
        xlim_left=None, xlim_right=None,
        fig_title="", save_fig=False, file_name=None):
    
    figure_for_noise = pyplot.figure(figsize=(fig_wd, fig_ht))
    plot_for_noise   = figure_for_noise.add_subplot(1, 1, 1)
    
    max_n_obss   = 0
    median_noise = []
    for source, noise_dict in noise_data.items():
        obs_ids = sorted(noise_dict.keys())
        units   = core.G3Units
        noise_levels = [noise_dict[obs_id] / (units.uK*units.arcmin) \
                        for obs_id in obs_ids]
        median_noise.append(numpy.nanmedian(noise_levels))
        obs_ids = [int(obs_id) for obs_id in obs_ids]
        if len(obs_ids) > max_n_obss:
            max_n_obss = len(obs_ids)
        
        el_dict = {"ra0hdec-44.75": "El 0", "ra0hdec-52.25": "El 1",
                   "ra0hdec-59.75": "El 2", "ra0hdec-67.25": "El 3"}
        cl_dict = {"ra0hdec-44.75": "#1f77b4", "ra0hdec-52.25": "#ff7f0e",
                   "ra0hdec-59.75": "#2ca02c", "ra0hdec-67.25": "#d62728"}
        label = el_dict[source]
        
        if is_integrated_noise:
            plot_for_noise.plot(
                range(1, len(obs_ids)+1), noise_levels, label=label,
                linewidth=0.4, marker=".", markersize=5.0,
                color=cl_dict[source], alpha=0.5)
        else:
            plot_for_noise.plot(
                obs_ids, noise_levels, label=label,
                linestyle="dotted", linewidth=0.4,
                marker=".", markersize=8.0,
                color=cl_dict[source], alpha=0.5)
    median_noise = numpy.median(median_noise)
    
    if is_integrated_noise:
        plot_for_noise.set_xscale("log")
        plot_for_noise.set_yscale("log")
        if median_noise > 300:
            plot_for_noise.set_ylim(bottom=5e0, top=1e3)
        else:
            plot_for_noise.set_ylim(bottom=5e0, top=1e3)
        plot_for_noise.set_xlim(
            left=0.9, right=10**numpy.ceil(numpy.log10(max_n_obss)))
    else:
        plot_for_noise.set_yscale("log")
        if median_noise > 300:
            plot_for_noise.set_ylim(bottom=8e1, top=1e3)
        else:
            plot_for_noise.set_ylim(bottom=8e1, top=1e3)
        if (xlim_left is not None) and (xlim_right is not None):
            margin = (xlim_right - xlim_left) * 0.05
            plot_for_noise.set_xlim(
                left=xlim_left-margin, right=xlim_right+6.0*margin)
    
    if not is_integrated_noise:
        if (xlim_left is None) or (xlim_right is None):
            xtick_locs   = plot_for_noise.get_xticks()
            xtick_labels = ["/".join(
                                str(std_processing.obsid_to_g3time(obs_id)).\
                                split(":")[0].split("-")[0:2]) \
                            for obs_id in xtick_locs]
            plot_for_noise.set_xticklabels(xtick_labels)
        else:
            xtick_locs, xtick_labels = \
                decide_on_xticks_from_obs_id_range(xlim_left, xlim_right)
            plot_for_noise.set_xticks(xtick_locs)
            plot_for_noise.set_xticklabels(xtick_labels)
    
    plot_for_noise.tick_params(direction="in", which="both", labelsize=tck_fs)
    plot_for_noise.grid(which="both", linestyle="dotted", linewidth=0.2)
    plot_for_noise.legend(loc="upper right", fontsize=tck_fs, framealpha=0.0)
        
    plot_for_noise.set_ylabel("Average  "+r"$\sqrt{C_l}$"+"  "+\
                              "in the "+"\n"+"ell range [3000, 5000]  "+\
                              r"$[{\mu}K \cdot arcmin]$"+"\n", fontsize=lbl_fs)
    if is_integrated_noise:
        plot_for_noise.set_xlabel(
            "\nNumber of observations", fontsize=lbl_fs)
    else:
        plot_for_noise.set_xlabel(
            "\nTime (UTC)", fontsize=lbl_fs)
    plot_for_noise.set_title(fig_title, fontsize=ttl_fs)
    
    figure_for_noise.tight_layout()
    
    if not save_fig:
        pyplot.show(figure_for_noise)
    else:
        pyplot.savefig(file_name, bbox_inches="tight")




def make_figure_for_pointing_discrepancies(
        discrepancy_data, ra_or_dec="ra",
        xlim_left=None, xlim_right=None,
        fig_title="", save_fig=False, file_name=None):
    
    figure_for_discrepancies = pyplot.figure(figsize=(fig_wd, fig_ht))
    plot_for_discrepancies   = figure_for_discrepancies.add_subplot(1, 1, 1)
    
    el_dict = \
        {"ra0hdec-44.75": "el0",
         "ra0hdec-52.25": "el1",
         "ra0hdec-59.75": "el2",
         "ra0hdec-67.25": "el3"}
    
    cl_dict = \
        {"ra0hdec-44.75": {"1": "aqua",   "2": "deepskyblue", "3": "darkturquoise"},
         "ra0hdec-52.25": {"1": "orange", "2": "gold",        "3": "khaki"},
         "ra0hdec-59.75": {"1": "lime",   "2": "greenyellow", "3": "forestgreen"},
         "ra0hdec-67.25": {"1": "red",    "2": "lightsalmon", "3": "maroon"}}
    
    for sub_field, sources_and_deltas in discrepancy_data.items():
        for source_ranking, deltas in sources_and_deltas.items():
            plot_for_discrepancies.plot(
                deltas["obs_ids"], deltas["delta_"+ra_or_dec+"s"],
                label=el_dict[sub_field]+" "+"src"+" "+source_ranking,
                linestyle="dotted", linewidth=0.4,
                marker=".", markersize=8.0,
                color=cl_dict[sub_field][source_ranking], alpha=0.9)
    
    plot_for_discrepancies.set_ylim(bottom=-45, top=45)
    
    if (xlim_left is None) or (xlim_right is None):
        xtick_locs   = plot_for_noise.get_xticks()
        xtick_labels = ["/".join(
                            str(std_processing.obsid_to_g3time(obs_id)).\
                            split(":")[0].split("-")[0:2]) \
                        for obs_id in xtick_locs]
        plot_for_discrepancies.set_xticklabels(xtick_labels)
    else:
        x_margin = (xlim_right - xlim_left) * 0.05
        plot_for_discrepancies.set_xlim(
            left=xlim_left-x_margin, right=xlim_right+6.0*x_margin)
        xtick_locs, xtick_labels = \
            decide_on_xticks_from_obs_id_range(xlim_left, xlim_right)
        plot_for_discrepancies.set_xticks(xtick_locs)
        plot_for_discrepancies.set_xticklabels(xtick_labels)

    plot_for_discrepancies.tick_params(
        labelsize=tck_fs, direction="in", which="both")
    plot_for_discrepancies.grid(
        which="both", linestyle="dotted", linewidth=0.2)
    plot_for_discrepancies.legend(
        loc="upper right", fontsize=tck_fs, framealpha=0.9)
    
    if ra_or_dec == "ra":
        plot_for_discrepancies.set_ylabel(
            "(Measured R.A. - True R.A.) "+"\n"+r"$\times$"+" cos(True Dec.)"+\
            " [arc second]\n", fontsize=lbl_fs)
    else:
        plot_for_discrepancies.set_ylabel(
            "(Measured Dec. - True Dec.) [arc second]\n", fontsize=lbl_fs)
    
    plot_for_discrepancies.set_xlabel("\nTime (UTC)", fontsize=lbl_fs)
    
    plot_for_discrepancies.set_title(fig_title, fontsize=ttl_fs)
    
    figure_for_discrepancies.tight_layout()
    
    if not save_fig:
        pyplot.show(figure_for_discrepancies)
    else:
        pyplot.savefig(file_name, bbox_inches="tight")




def make_figure_for_flagging_statistics(
        map_frame,
        xlim_left, xlim_right,
        fig_title, save_fig, file_name):
    
    figure_for_statistics = pyplot.figure(figsize=(fig_wd, fig_ht))
    plot_for_statistics   = figure_for_statistics.add_subplot(1, 1, 1)
    
    cl_dict  = {"BadHk"                  : "#1f77b4",
                "PostCalibrationNaNs"    : "#ff7f0e",
                "Latched"                : "#2ca02c",
                "Overbiased"             : "#d62728",
                "Oscillating"            : "#9467bd",
                "Glitchy"                : "#8c564b",
                "BadCalSn"               : "#e377c2",
                "MissingFluxCalibration" : "#7f7f7f",
                "UnphysicalLowVariance"  : "#bcbd22",
                "BadWeight"              : "#1fbecf",
                "Others"                 : "#000000",
                "TotalNotFlagged"        : "black",
                "TotalRemoved"           : "black"}
    
    lbl_dict = {"BadHk"                  : "Bad HK"     ,
                "PostCalibrationNaNs"    : "NaNs"       ,
                "Latched"                : "Latched"    ,
                "Overbiased"             : "Saturated"  ,
                "Oscillating"            : "Oscillating",
                "Glitchy"                : "Glitchy"    ,
                "BadCalSn"               : "Bad CalSN"  ,
                "MissingFluxCalibration" : "No Fluxcal" ,
                "UnphysicalLowVariance"  : "Low Var."   ,
                "BadWeight"              : "Bad Weight" ,
                "Others"                 : "Others"     ,
                "TotalNotFlagged"        : "Survived"   ,
                "TotalRemoved"           : "Removed"}
    
    mkr_dict = {"TotalNotFlagged"        : "$\u266B$",
                "TotalRemoved"           : "x"}
    
    
    for flag_reason in ["BadHk",
                        "PostCalibrationNaNs",
                        "Latched",
                        "Overbiased",
                        "Oscillating",
                        "Glitchy",
                        "BadCalSn",
                        "MissingFluxCalibration",
                        "UnphysicalLowVariance",
                        "BadWeight",
                        "Others",
                        "TotalNotFlagged",
                        "TotalRemoved"]:
        
        key_name = "FlaggingStatisticsAverageNumberOf"+flag_reason
        averages = map_frame[key_name]
        
        times_and_averages = []
        for sub_field, oids_and_avgs in averages.items():
            for oid, avg in oids_and_avgs.items():
                times_and_averages.append((int(oid), avg))
        
        times_and_averages = sorted(times_and_averages, key=itemgetter(0))
        times    = [entry[0] for entry in times_and_averages]
        averages = [entry[1] for entry in times_and_averages]
        
        if flag_reason in mkr_dict.keys():
            marker = mkr_dict[flag_reason]
        else:
            marker = "."
        
        plot_for_statistics.plot(times, averages,
                                 label=lbl_dict[flag_reason],
                                 linestyle="dotted", linewidth=0.4,
                                 marker=marker, markersize=5.0,
                                 color=cl_dict[flag_reason], alpha=0.9)
    
    plot_for_statistics.set_yscale("symlog", linthreshy=50,
                                   subsy=[2, 3, 4, 5, 6, 7, 8, 9])
    plot_for_statistics.set_ylim(bottom=-2, top=5e3)
    
    if (xlim_left is None) or (xlim_right is None):
        xtick_locs   = plot_for_statistics.get_xticks()
        xtick_labels = ["/".join(
                            str(std_processing.obsid_to_g3time(obs_id)).\
                            split(":")[0].split("-")[0:2]) \
                        for obs_id in xtick_locs]
        plot_for_statistics.set_xticklabels(xtick_labels)
    else:
        x_margin = (xlim_right - xlim_left) * 0.05
        plot_for_statistics.set_xlim(
            left=xlim_left-x_margin, right=xlim_right+6.0*x_margin)
        xtick_locs, xtick_labels = \
            decide_on_xticks_from_obs_id_range(xlim_left, xlim_right)
        plot_for_statistics.set_xticks(xtick_locs)
        plot_for_statistics.set_xticklabels(xtick_labels)
    
    plot_for_statistics.grid(
        axis="x", which="major", linestyle="dotted", linewidth=0.2)
    plot_for_statistics.grid(
        axis="y", which="both", linestyle="dotted", linewidth=0.2)
    plot_for_statistics.tick_params(
        labelsize=tck_fs, direction="in", which="both")
    plot_for_statistics.legend(
        loc="upper right", fontsize=tck_fs, framealpha=0.9)
    
    plot_for_statistics.set_ylabel(
        "Avg. (over all scans) number"+"\n", fontsize=lbl_fs)
    plot_for_statistics.set_xlabel(
        "\nTime (UTC)", fontsize=lbl_fs)
    plot_for_statistics.set_title(
        fig_title, fontsize=ttl_fs)
    
    figure_for_statistics.tight_layout()
    
    if not save_fig:
        pyplot.show(figure_for_statistics)
    else:
        pyplot.savefig(file_name, bbox_inches="tight")


# ==============================================================================





# ==============================================================================
# Load data from each G3 file and make figures
# ------------------------------------------------------------------------------


for g3_file in arguments.input_files:
    
    print()
    print("# ------------------------------------------------- ")
    print("# Loading data from", g3_file, "..."  )
    print("# ------------------------------------------------- ")
    print()
    
    try:
        frames = list(core.G3File(g3_file))
    except RuntimeError:
        print("There seems to be a problem with the G3 file.")
        print("This observation will be skipped!")
        print()
        continue
    
    print("These are the frames stored in the file:")
    print()
    for frame in frames:
        core.Dump(frame)
    print()
    
    
    if not arguments.coadded_data:    
        for frame in frames:
            if frame.type == core.G3FrameType.Observation:
                obs_frame = frame
    
        
    for frame in frames:
        if frame.type != core.G3FrameType.Map:
            continue
        elif frame["Id"] not in valid_map_ids:
            continue
        else:
            map_frame = frame
            map_id    = frame["Id"]
            print()
            print("----------------------------------------------------------")
            print(" Found a Map frame with ID", map_id+"!")
            print(" Figures will be made for this frame.")
            print("----------------------------------------------------------")
            print()
            
            sky_map    = get_sky_map(map_frame)
            weight_map = get_weight_map(map_frame)
            cmb_map_mK = (sky_map / weight_map) / core.G3Units.mK
            
            print("* Some basic properties of the map:")
            print("Map shape    :", sky_map.shape)
            print("Map units    :", sky_map.units)
            print("Map size     :", sky_map.size)
            print("x resolution :", sky_map.x_res/core.G3Units.arcmin,
                  "arc minutes")
            print("y resolution :", sky_map.y_res/core.G3Units.arcmin,
                  "arc minutes")
            print()
            
            
            print("Preparing figure titles and file names...")
            print()
            
            
            # - Prepare the first (and maybe also the second) line of the title
            
            if not arguments.coadded_data:
                source = obs_frame["SourceName"]
                obs_id = str(obs_frame["ObservationID"])
                date   = str(obs_frame["ObservationStart"]).split(".")[0]
                res    = str(cmb_map_mK.x_res/core.G3Units.arcmin)+"' resolution"
                title_a = source+"  "+obs_id+" ("+date+") "+"   "+map_id+\
                          " "+map_type_str+" ("+resoluti+")"
                if arguments.make_figures_for_weights:
                    ttl_a_w  = source+"  "+obs_id+" ("+date+") "+"   "+map_id+\
                               " "+weight_type_str+" ("+resoluti+")"
            else:
                source  = "1500 square degrees field"
                obs_ids = []
                for obs_ids_one_field \
                in  map_frame["CoaddedObservationIDs"].values():
                    for obs_id in obs_ids_one_field:
                        obs_ids.append(obs_id)
                start_dt = str(std_processing.obsid_to_g3time(
                               numpy.min(obs_ids))).split(":")[0]
                end_dt   = str(std_processing.obsid_to_g3time(
                               numpy.max(obs_ids))).split(":")[0]
                obs_id   = str(numpy.min(obs_ids))+"_to_"+\
                           str(numpy.max(obs_ids))
                el_dict  = {"ra0hdec-44.75": "El 0", "ra0hdec-52.25": "El 1",
                            "ra0hdec-59.75": "El 2", "ra0hdec-67.25": "El 3"}
                n_obss   = ["{"+str(len(obss))+" "+el_dict[source]+"}" \
                            for source, obss \
                            in  map_frame["CoaddedObservationIDs"].items()]
                n_obss   = " ".join(n_obss)
                resolu   = str(cmb_map_mK.x_res/core.G3Units.arcmin)+"' res."
                title_a  = source+"   "+map_id+" coadded "+map_type_str+"s  "+\
                           "("+resolu+")"+"\n"+\
                           "(data from observations taken "+\
                           "between "+start_dt+" and "+end_dt+":"+"\n"+n_obss+")"
                if arguments.make_figures_for_weights:
                    ttl_a_w  = source+"   "+map_id+" coadded "+\
                               weight_type_str+"  "+\
                               "("+resolu+")"+"\n"+\
                               "(data from observations taken "+\
                               "between "+start_dt+" and "+end_dt+":"+"\n"+n_obss+")"
            
            
            # - Prepare the last line in the title
            
            median  = str(numpy.round(numpy.nanmedian(cmb_map_mK)))
            pctl_10 = str(numpy.round(numpy.nanpercentile(cmb_map_mK, 10), 3))
            pctl_90 = str(numpy.round(numpy.nanpercentile(cmb_map_mK, 90), 3))
            title_b = "(Some statistics:  "+"10th pctl. = "+pctl_10+" mK,  "+\
                      "Median = "+median+" mK,  "+"90th pctl. = "+pctl_90+" mK)"
            if arguments.make_figures_for_weights:
                non_0_wts = numpy.asarray(weight_map)\
                            [numpy.where(numpy.asarray(weight_map)!=0.0)]
                median_w  = str(numpy.round(numpy.nanmedian(non_0_wts)))
                pctl_10_w = str(numpy.round(numpy.nanpercentile(
                                non_0_wts, 10), 0))
                pctl_90_w = str(numpy.round(numpy.nanpercentile(
                                non_0_wts, 90), 0))
                ttl_b_w   = "(Some statistics:  "+"10th pctl. = "+pctl_10_w+\
                            ",  "+"Median = "+median_w+",  "+\
                            "90th pctl. = "+pctl_90_w+")"
            
            
            # - Assemble the pieces for the title
            
            fig_title = title_a + "\n" + title_b + "\n"
            if arguments.make_figures_for_weights:
                fig_ttl_w = ttl_a_w + "\n" + ttl_b_w + "\n"
            
            
            # - Prepare the colorbar label
            
            arb_perc   = 1
            low_v      = numpy.nanpercentile(cmb_map_mK,     arb_perc)
            high_v     = numpy.nanpercentile(cmb_map_mK, 100-arb_perc)
            larger_abs = numpy.max([numpy.absolute(low_v),
                                    numpy.absolute(high_v)])
            if arguments.make_figures_for_weights:
                low_v_w   = numpy.nanpercentile(weight_map,     arb_perc)
                high_v_w  = numpy.nanpercentile(weight_map, 100-arb_perc)
                lar_abs_w = numpy.max([numpy.absolute(low_v_w),
                                       numpy.absolute(high_v_w)])
            
            cbar_title = "\n"+r"$mK_{CMB}$"+"  ( percentile range covered is "+\
                         "["+str(arb_perc)+", "+str(100-arb_perc)+"] )"
            if arguments.make_figures_for_weights:
                cbar_ttl_w = "\n"+"Arb."+" (percentile range covered is "+\
                             "["+str(0)+", "+str(100-arb_perc)+"])"
            
            
            # - Prepare the file name
            
            if "90" in map_id:
                filename = source+"-"+obs_id+"-0"+map_id+"_"+\
                           map_type_str.replace(" ", "_")
            else:
                filename = source+"-"+obs_id+"-"+map_id+"_"+\
                           map_type_str.replace(" ", "_")
            filename = filename + ".png"
            filename = arguments.directory_to_save_figures+"/"+filename
            filename = filename.replace(" ", "_")
            if arguments.make_figures_for_weights:
                if map_id == "90GHz":
                    fina_wts = source+"-"+obs_id+"-0"+map_id+"_"+\
                               weight_type_str.replace(" ", "_")
                else:
                    fina_wts = source+"-"+obs_id+"-"+map_id+"_"+\
                               weight_type_str.replace(" ", "_")
                fina_wts = fina_wts+ ".png"
                fina_wts = arguments.directory_to_save_figures+"/"+fina_wts
                fina_wts = fina_wts.replace(" ", "_")
            
            if arguments.simpler_file_names:
                filename = map_id+"-"+map_type_str.replace(" ", "_")+".png"
                fina_wts = map_id+"-"+map_type_str.replace(" ", "_")+\
                           "_weights.png"           
                filename = arguments.directory_to_save_figures+"/"+filename
                fina_wts = arguments.directory_to_save_figures+"/"+fina_wts
            
            
            # - Finally make figures
            
            print("Now actually making figures...")
            print()
            
            if (arguments.map_values_upper_limit != None) and \
               (arguments.map_values_lower_limit != None):
                vmin_for_map = arguments.map_values_lower_limit
                vmax_for_map = arguments.map_values_upper_limit
            else:
                vmin_for_map = -1.0 * larger_abs
                vmax_for_map =  1.0 * larger_abs
            
            visualize_a_field_map(
                cmb_map_mK,
                cmb_map_mK.x_res, cmb_map_mK.y_res, fig_wd, fig_ht, "auto",
                "gray", cbar_label=cbar_title, title=fig_title,
                vmin=vmin_for_map, vmax=vmax_for_map,
                no_xy_labels=True,
                save_fig=(not arguments.only_show_figures),
                filename=filename)
            
            if arguments.make_figures_for_weights:
                visualize_a_field_map(
                    weight_map,
                    cmb_map_mK.x_res, cmb_map_mK.y_res, fig_wd, fig_ht, "auto",
                    "gray", vmin=0.0, vmax=1.0*lar_abs_w,
                    cbar_label=cbar_ttl_w, title=fig_ttl_w,
                    no_xy_labels=True,
                    save_fig=(not arguments.only_show_figures),
                    filename=fina_wts)
                
                dec_values = numpy.linspace(-41.0, -71.0, 9)
                dec_labels = dec_values
                pixel_ids  = [cmb_map_mK.angle_to_pixel(
                                  0.0*core.G3Units.rahr, dec_value) \
                              for dec_value in dec_values*core.G3Units.deg]
                y_indices  = [pixel_id//cmb_map_mK.shape[1]-1 \
                              for pixel_id in pixel_ids]
                
                make_histogram_and_cross_section_plot_of_a_field_map(
                    weight_map,
                    hist_xlabel="Normalized TT Weight [unitless] ",
                    hist_title=\
                        "Histogram of all non-zero values in the weight map",
                    cr_sec_xticks=y_indices,
                    cr_sec_xticklabels=dec_labels,
                    cr_sec_xlabel="Declination [degree]",
                    cr_sec_title=\
                        map_id+" "+map_type_str+\
                        "   Cross section view of the weight map "+\
                        "along the RA=0h contour"+"\n",
                    fig_title="More plots for the weight map",
                    save_fig=(not arguments.only_show_figures),
                    filename=fina_wts[0:-4]+"_cross_sectional_view.png")
            
            
            # - Also make a figure for noise levels
            
            if arguments.make_figures_for_noise_levels:
                if arguments.integrated_noise:
                    noise_data = map_frame["NoiseFromCoaddedMaps"]
                else:
                    noise_data = map_frame["NoiseFromIndividualMaps"]
                
                file_name  = map_id+"-"+map_type_str.replace(" ", "_")
                if arguments.integrated_noise:
                    file_name += "_noise_levels_from_running_coadds.png"
                else:
                    file_name += "_noise_levels_from_individual_maps.png"
                
                el_dict = {"ra0hdec-44.75": "El 0", "ra0hdec-52.25": "El 1",
                           "ra0hdec-59.75": "El 2", "ra0hdec-67.25": "El 3"}
                n_obss  = ["{"+str(len(obss))+" "+el_dict[source]+"}" \
                           for source, obss \
                           in  map_frame["CoaddedObservationIDs"].items()]
                n_obss  = " ".join(n_obss)
                
                fig_title_prefix = map_id + " " + map_type_str
                
                if arguments.integrated_noise:
                    fig_title_mid = "Noise levels of running coadded maps"
                else:
                    fig_title_mid = "Noise levels from individual observations"
                
                fig_title = \
                    fig_title_prefix + "   " + fig_title_mid + "\n"
                
                make_figure_for_noise_levels(
                    noise_data,
                    is_integrated_noise=arguments.integrated_noise,
                    xlim_left=arguments.left_xlimit_for_figures,
                    xlim_right=arguments.right_xlimit_for_figures,
                    fig_title=fig_title,
                    save_fig=(not arguments.only_show_figures),
                    file_name=arguments.directory_to_save_figures+"/"+file_name)
            
            
            # - Also make a figure for pointing discrepancies
            
            if arguments.make_figures_for_pointing_discrepancies:
                reorganized_data = {}
                for sub_field in map_frame["DeltaRasFromSources1"].keys():
                    reorganized_data[sub_field] = \
                        {"1": {"obs_ids": [], "delta_ras": [], "delta_decs": []},
                         "2": {"obs_ids": [], "delta_ras": [], "delta_decs": []},
                         "3": {"obs_ids": [], "delta_ras": [], "delta_decs": []}}
                    for coordinate in ["Ra", "Dec"]:
                        for ranking in ["1", "2", "3"]:
                            original_key = "Delta" + coordinate + \
                                           "sFromSources" + ranking
                            new_key      = "delta_" + coordinate.lower() + "s"
                            for obs_id, delta \
                            in  map_frame[original_key][sub_field].items():
                                if coordinate == "Ra":
                                    reorganized_data[sub_field][ranking]\
                                                    ["obs_ids"].append(obs_id)
                                delta /= core.G3Units.arcsec
                                reorganized_data[sub_field][ranking]\
                                                [new_key].append(delta)
                    
                    for coordinate in ["Ra", "Dec"]:
                        file_name  = map_id + "-" + \
                                     map_type_str.replace(" ", "_")
                        file_name += "_delta_" + coordinate + \
                                     "s_from_point_sources.png"
                        fig_title  = map_id+" "+map_type_str + " " + \
                                     "Difference between " + \
                                     "measured and true "+ coordinate + \
                                     " for several point sources" + "\n"
                    
                        make_figure_for_pointing_discrepancies(
                            reorganized_data,
                            ra_or_dec=coordinate.lower(),
                            xlim_left=arguments.left_xlimit_for_figures,
                            xlim_right=arguments.right_xlimit_for_figures,
                            fig_title=fig_title,
                            save_fig=(not arguments.only_show_figures),
                            file_name=arguments.directory_to_save_figures+\
                                      "/"+file_name)
            
            
            # - Also make a figure for flagging statistics
            
            if arguments.make_figures_for_flagging_statistics:
                file_name  = map_id + "-" + \
                             map_type_str.replace(" ", "_")
                file_name += "_average_numbers_of_flagged_detectors.png"
                fig_title  = map_id + " " + map_type_str + "   " + \
                             "Average number of flagged detectors " +\
                             "for each reason during each observation" + "\n"
                
                make_figure_for_flagging_statistics(
                    map_frame,
                    xlim_left=arguments.left_xlimit_for_figures,
                    xlim_right=arguments.right_xlimit_for_figures,
                    fig_title=fig_title,
                    save_fig=(not arguments.only_show_figures),
                    file_name=arguments.directory_to_save_figures+\
                              "/"+file_name)
            
            
            del sky_map
            del weight_map
            del map_frame
            gc.collect()
        
            if arguments.make_figures_for_one_frame_only:
                break

    del frames
    gc.collect()
        
    print()


# ==============================================================================

