# ============================================================================ #
#  This script is intended to be called by                                     #
#  spt_dq_system/summaryplot/cache_field_maps.py.                              #
#  However, it can also be used as an independent script (that can be called   #
#  from the command line or imported by another python script) that co-add     #
#  maps and do some analysis.                                                  #
#                                                                              #
#  When cache_field_maps.py calls this script, the former specifies            #
#  which input g3 files contain maps to be added together,                     #
#  which analyses should be done on individual or coadded maps, and so on.     #
#                                                                              #
#  Then, this script constructs a pipeline that co-adds maps and perfomes      #
#  the requested analyses. The types of analyses that are available are the    #
#  following:                                                                  #
#                                                                              #
#    Analyses that do not require actual map data:                             #
#                                                                              #
#    - Collecting flagging statistics                                          #
#        The script can calculate the average number of detectors              #
#        (average over all scans of an observation) that were flagged and      #
#        not used in the map-making, the average number of detectors that      #
#        were flagged due to a specific reason, and the average number of      #
#        detectors that were not flagged.                                      #
#                                                                              #
#    - Collecting temperature calibration factors                              #
#        The script can calculate the average value of the pW/K conversion     #
#        factors that were used to convert timestreams' units from power to    #
#        temperature during the map-making. The average is calculated for      #
#        each wafer and the whole detector array.                              #
#                                                                              #
#    - Calculating changes in responsivity                                     #
#        The script can collect the results from calibrator observations that  #
#        were taken at the bottom, middle, and top of each sub-field as part   #
#        part of the schedule for observing that field. Specifically, it       #
#        calculates the median response to the calibrator for each wafer and   #
#        and the whole detector array. In addtion, the median of the           #
#        fractional change in detectors' response to the calibrator at the     #
#        top with respect to the bottom of the field is calculated.            #
#                                                                              #
#    Analyses that do require map data:                                        #
#                                                                              #
#    - Calculating basic statistics of temperature and weight maps             #
#        The script can calculate standard deviations of the values in         #
#        temperature maps, mean values of weight maps, and numbers of          #
#        pixels in weight maps whose weights are higher than some normal       #
#        values. When these statistics are calculated, the pixels that are     #
#        close to the periphery of a map are ignored. Specifically, these      #
#        pixels are the ones whose RAs and DECs meet these conditions:         #
#        RA > 48 deg or < -48 deg, DEC > field_center_dec + 3 or < -3.         #
#        In addition, point sources are ignored when the standard deviations   #
#        are calculated.                                                       #
#                                                                              #
#    - Calculating pointing discrepancies                                      #
#        The script can calculate the differences between the true positions   #
#        and the measured ones of three bright point sources in each of the    #
#        four sub-fields. The measured positions are results of fitting a      #
#        2D gaussian to a small (4 arcmin. by 4 arcmin.) region of a map       #
#        that is roughly centered on each point source.                        #
#                                                                              #
#    - Comparing SPT maps with Planck maps                                     #
#        The script can compare an individual SPT map with a Planck map by     #
#        calculating the cross spectrum between the two, calculating the       #
#        autospectrum of the Planck map, dividing the first spectrum by the    #
#        second one, and calculating the average of the ratios in the ell      #
#        range [750, 1250], i.e., average of                                   #
#        Cl_(SPT x Planck) / Cl_(Planck_x_Planck).                             #
#        (the Planck in the numerator refers to the Planck full mission map,   #
#         and the Plancks in the denominator refer to two half mission maps;   #
#         this way, the denominator is less biased by Planck noise)            #
#        In order to save time and memory, only a 5 degree by 5 degree patch   #
#        centered around the corresponding sub-field center is used instead    #
#        of the entire sub-field.                                              #
#                                                                              #
#    - Calculating noise levels                                                #
#        The script can estimate noise levels of maps. For maps of             #
#        individual observations, if both a map made from only left-going      #
#        scans of an observation and one from right-going scans are            #
#        available, then the noise level is defined as the square root of      #
#        the Cl of the autospectrum of the difference map in the ell range     #
#        [3000, 5000]. If the maps are not split by the direction, then the    #
#        spectra of individual maps are just used. For coadded maps, if maps   #
#        are split by direction, then spectra of coadded 'left-going' maps     #
#        minus coadded 'right-going' maps are taken. If maps are not split     #
#        by direction, still two coadded maps will be maintained by            #
#        alternately adding individual maps to one or the other pile, and      #
#        the difference between them will be used for calculating the          #
#        noise. Once again, in order to save time and memory, only             #
#        the 5 degree by 5 degree patches mentioned above are used.            #
#                                                                              #
#                                                                              #
#  The script mainly has two parts. The first part consists of definitions of  #
#  functions that are relevent to performing these analyses, and the second    #
#  part defines a pipeline module (AnalyzeAndCoaddMaps) that adds maps         #
#  together and calls those functions to perform the analyses.                 #
#                                                                              #
#  The script is currently capable of dealing with temperature-only            #
#  map frames only. More code would need to be written if we were to use it    #
#  to co-add and analyze T, Q, and U maps as well.                             #
#                                                                              #
# ============================================================================ #


import argparse
import datetime
import logging
import os
import sys
import glob
import numpy
import pickle

from spt3g import core
from spt3g import mapmaker
from spt3g import util
from spt3g import maps
from spt3g import std_processing
from spt3g import mapspectra
from spt3g import pointing
from scipy import signal



# ==============================================================================
# Define functions needed to perform the analyses
# ------------------------------------------------------------------------------


def get_season_based_on_fields(some_fields):
    
    winter_fields = ["ra0hdec-44.75", "ra0hdec-52.25",
                     "ra0hdec-59.75", "ra0hdec-67.25"]
    summer_fields = ["ra5hdec-24.5" , "ra5hdec-31.5",
                     "ra5hdec-38.5" , "ra5hdec-45.5",
                     "ra5hdec-52.5" , "ra5hdec-59.5"]
    summer_fields_p = ["ra5hdec-29.75", "ra5hdec-33.25",
                       "ra5hdec-36.75", "ra5hdec-40.25"]
    summer_fields_b = ["ra1h40dec-29.75", "ra1h40dec-33.25",
                       "ra1h40dec-36.75", "ra1h40dec-40.25"]
    summer_fields_c = ["ra12h30dec-29.75", "ra12h30dec-33.25",
                       "ra12h30dec-36.75", "ra12h30dec-40.25"]
    
    if set(some_fields) <= set(winter_fields):
        return "winter"
    elif set(some_fields) <= set(summer_fields):
        return "summer"
    elif set(some_fields) <= set(summer_fields_p):
        return "summerp"
    elif set(some_fields) <= set(summer_fields_b):
        return "summerb"
    elif set(some_fields) <= set(summer_fields_c):
        return "summerc"
    else:
        raise Exception("Not clear what season to use!")




def add_two_map_frames(
        frame_one, frame_two, t_only=True,
        remove_weights_beforehand=False,
        multiply_2nd_maps_by=1.0,
        multiply_2nd_weights_by=1.0,
        remove_weights_afterward=False,
        divide_result_by=1.0,
        record_weights=True,
        logfun=print,
        iy=6000, ix=9000):
    
    logfun("$$$")
    logfun("$$$ adding two map frames ...")
    logfun("$$$ PARAMETERS:")
    for k, v in locals().items():
        if k in ["frame_one", "frame_two"]:
            continue
        logfun("    %-25s : %s" %(k, v))
    tu = core.G3Units.mK
    
    if t_only:
        
        if not frame_one["T"].weighted:
            maps.ApplyWeights(frame_one)
        if not frame_two["T"].weighted:
            maps.ApplyWeights(frame_two)
        
        logfun("$$$ use x = %d and y = %d "
               "for the indices of a particular pixel" %(ix, iy))
        logfun("$$$ t1 initial : %+7.3e"
               %(numpy.asarray(frame_one["T"])[iy][ix]/tu))
        logfun("$$$ t2 initial : %+7.3e"
               %(numpy.asarray(frame_two["T"])[iy][ix]/tu))
        logfun("$$$ w1 initial : %+7.3e [1/mK^2]"
               %(numpy.asarray(frame_one["Wunpol"].TT)[iy][ix]*tu*tu))
        logfun("$$$ w2 initial : %+7.3e [1/mK^2]"
               %(numpy.asarray(frame_two["Wunpol"].TT)[iy][ix]*tu*tu))
        
        # - First, remove weights if necessary
        
        if remove_weights_beforehand:
            maps.RemoveWeights(frame_one, zero_nans=True)
            maps.RemoveWeights(frame_two, zero_nans=True)
        
        logfun("$$$ t1 after potential weights removal : %+7.3e"
               %(numpy.asarray(frame_one["T"])[iy][ix]/tu))
        logfun("$$$ t2 after potential weights removal : %+7.3e"
               %(numpy.asarray(frame_two["T"])[iy][ix]/tu))
        
        maps.CompactMaps(frame_one, zero_nans=True)
        maps.CompactMaps(frame_two, zero_nans=True)
        
        # - Then, add field maps (and weight maps if applicable) together
        
        new_frame = core.G3Frame(core.G3FrameType.Map)

        new_frame["T"] = frame_one["T"] + \
                         frame_two["T"] * multiply_2nd_maps_by        
        new_frame["Wunpol"] = maps.G3SkyMapWeights()
        new_frame["Wunpol"].TT = frame_one["Wunpol"].TT + \
                                 frame_two["Wunpol"].TT * multiply_2nd_weights_by
        
        logfun("$$$ t resultant : %+7.3e"
               %(numpy.asarray(new_frame["T"])[iy][ix]/tu))
        logfun("$$$ w resultant : %+7.3e [1/mK^2]"
               %(numpy.asarray(new_frame["Wunpol"].TT)[iy][ix]*tu*tu))
        
        # - After that, remove weights if necessary
        
        if remove_weights_afterward:
            maps.RemoveWeights(new_frame, zero_nans=True)
        
        logfun("$$$ t after potential weights removal : %+7.3e [mK]"
               %(numpy.asarray(new_frame["T"])[iy][ix]))
        
        # - Furthur divide the map by some number if requested
        
        if divide_result_by != 1.0:
            new_frame["T"] = new_frame.pop("T") / divide_result_by
        
        logfun("$$$ t after division : %+7.3e [mK]"
               %(numpy.asarray(new_frame["T"])[iy][ix]/tu))
        
        # - Finally, check some numbers again
                
        logfun("$$$ t final : %+7.3e"
               %(numpy.asarray(new_frame["T"])[iy][ix]/tu))
        if record_weights:
            logfun("$$$ w final : %+7.3e [1/mK^2]"
                   %(numpy.asarray(new_frame["Wunpol"].TT)[iy][ix]*tu*tu))
        else:
            logfun("$$$ w final : N/A")
        maps.CompactMaps(new_frame, zero_nans=True)
    
    else:
        raise Exception("The case where map frames contain temperature "
                        "as well as polarization maps still needs to be "
                        "handled!")
    
    logfun("$$$ FRAME TO BE RETURNED:")
    for k, v in new_frame.iteritems():
        logfun("$$$   %-10s : %s" %(k, v))
    logfun("$$$ ... adding finished.")
    logfun("$$$")
    
    return new_frame




def create_new_map_frame_with_smaller_region(
        map_frame, center_ra, center_dec, length_x, length_y, t_only=True):
    
    if map_frame is None:
        return None
    
    new_map_frame = core.G3Frame(core.G3FrameType.Map)

    if t_only:
        
        map_parameters = {"x_len": int(length_x/map_frame["T"].res),
                          "y_len": int(length_y/map_frame["T"].res),
                          "res"  : map_frame["T"].res,
                          "proj" : map_frame["T"].proj,
                          "alpha_center": center_ra,
                          "delta_center": center_dec,
                          "coord_ref"   : map_frame["T"].coord_ref,
                          "pol_type"    : map_frame["T"].pol_type}
        
        smaller_field_map = maps.FlatSkyMap(**map_parameters)
        maps.reproj_map(map_frame["T"], smaller_field_map, 1)
        new_map_frame["T"] = smaller_field_map
        
        if "Wunpol" in map_frame.keys():
            smaller_weight_map = maps.FlatSkyMap(**map_parameters)
            maps.reproj_map(
                map_frame["Wunpol"].TT, smaller_weight_map, 1)
            new_map_frame["Wunpol"] = maps.G3SkyMapWeights()
            new_map_frame["Wunpol"].TT = smaller_weight_map
    
    return new_map_frame




def get_wafer_average_over_scans(
        bolo_names_from_each_scan, bolo_props_map,
        desired_band, desired_wafers):
    
    # * This one is called by the function below
    
    n_bolos_from_each_scan = {wafer: [] for wafer in desired_wafers}
    
    for bolo_names_from_one_scan in bolo_names_from_each_scan.values():
        n_bolos_from_one_scan = {wafer: 0 for wafer in desired_wafers}
        available_bolos = set(bolo_props_map.keys()) & \
                          set(bolo_names_from_one_scan)
        for bn in available_bolos:
            band = bolo_props_map[bn].band / core.G3Units.GHz
            if not numpy.isfinite(band):
                continue
            band = str(int(band))
            if band == desired_band:
                wafer = bolo_props_map[bn].wafer_id.upper()
                if wafer in desired_wafers:
                    n_bolos_from_one_scan[wafer] += 1
                    n_bolos_from_one_scan["AllBolos"] += 1
        
        for wafer, n_bolos in n_bolos_from_one_scan.items():
            n_bolos_from_each_scan[wafer].append(n_bolos)
    
    return {wafer: numpy.mean(ns_bolos) \
            for wafer, ns_bolos in n_bolos_from_each_scan.items()}




def collect_averages_from_flagging_info(
        pipe_info_frame, recognized_reasons,
        desired_wafers, desired_band, bolo_props_map):
    
    averages = {wafer: {} for wafer in desired_wafers}
    
    # - Collect average numbers of bolos not flagged
    
    bolos_not_flagged_each_scan = pipe_info_frame["SurvivingBolos"]
    # * This is a core.G3MapVectorString instance
    
    total_n_scans = numpy.max([len(scans) for reason, scans \
                               in pipe_info_frame["DroppedBolos"].items()])
    # * In case more scans are recorded in the DroppedBolos dictionary
    
    n_from_survival = len(bolos_not_flagged_each_scan.keys())
    diff_n_scans = total_n_scans - n_from_survival
    for i in range(diff_n_scans):
        bolos_not_flagged_each_scan["Dummy"+str(i)] = core.G3VectorString()
    
    avg_n_okay_bolos_each_wafer = \
        get_wafer_average_over_scans(
            bolos_not_flagged_each_scan,
            bolo_props_map, desired_band, desired_wafers)
    for wafer, average in avg_n_okay_bolos_each_wafer.items():
        averages[wafer]["TotalNotFlagged"] = average
    
    # - Collect average numbers of bolos flagged for each reason
    
    # -- Gather the numbers from each reason
    
    all_bolos_flagged_each_scan = pipe_info_frame["DroppedBolos"]
    # * This is a core.G3MapVectorVectorString instance
    
    for flagging_reason, bolos_flagged_each_scan \
    in  all_bolos_flagged_each_scan.items():
        bl_flgd_ea_sc_reorganized = core.G3MapVectorString()
        for scan_number, list_of_bolos in enumerate(bolos_flagged_each_scan):
            bl_flgd_ea_sc_reorganized[str(scan_number)] = list_of_bolos
        avg_bad_bolos_each_wafer = \
            get_wafer_average_over_scans(
                bl_flgd_ea_sc_reorganized,
                bolo_props_map, desired_band, desired_wafers)
        for wafer, average in avg_bad_bolos_each_wafer.items():
            averages[wafer][flagging_reason] = average
    
    # -- Collectively call unrecognized reasons "Others"
    #    and assign 0 to any recognized reason that is absent
    
    for wafer, reasons_and_numbers in averages.items():
        averages[wafer]["Others"] = 0.0
        unrecognized_reasons = []
        for reason, number in reasons_and_numbers.items():
            if reason not in recognized_reasons:
                averages[wafer]["Others"] += number
                unrecognized_reasons.append(reason)
        for reason in unrecognized_reasons:
            averages[wafer].pop(reason)
        for reason in recognized_reasons:
            if reason not in reasons_and_numbers.keys():
                averages[wafer][reason] = 0.0
    
    # - Collect average numbers of bolos removed
    
    all_bolos_flagged_each_scan = pipe_info_frame["DroppedBolos"]
    n_scans = numpy.max([len(bl_flgd_ea_sc) for flg_rsn, bl_flgd_ea_sc \
                         in all_bolos_flagged_each_scan.items()])
    bolos_removed_each_scan = core.G3MapVectorString()
    for i in range(n_scans):
        bolos_removed_one_scan = set([])
        for flagging_reason, bolos_flagged_each_scan \
        in  all_bolos_flagged_each_scan.items():
            try:
                bolos_removed_one_scan = \
                    bolos_removed_one_scan | set(bolos_flagged_each_scan[i])
            except IndexError:
                pass
        bolos_removed_each_scan[str(i)] = \
            core.G3VectorString(bolos_removed_one_scan)
    avg_n_removed_bolos_each_wafer = \
        get_wafer_average_over_scans(
            bolos_removed_each_scan, bolo_props_map,
            desired_band, desired_wafers)
    for wafer, average in avg_n_removed_bolos_each_wafer.items():
        averages[wafer]["TotalRemoved"] = average
    
    # * "TotalNotFlagged", "TotalRemoved", and "Others" keys
    #   need to be changed accordingly if the corresponding names
    #   in the initialization of AnalyzeAndCoaddMaps are changed.
    
    return averages




def group_suviving_detectors_by_wafer(
        desired_band, desired_wafers, cal_frame, pipe_info_frame, logfun=print):
    
    bolos_used_each_scan = pipe_info_frame["SurvivingBolos"]
    bolos_and_times_used = {}
    
    for scan_number, bolos_that_scan in bolos_used_each_scan.iteritems():
        for bolo in bolos_that_scan:
            bolos_and_times_used[bolo] = bolos_and_times_used.get(bolo, 0) + 1
    
    n_scans = len(bolos_used_each_scan.keys())
    mostly_surviving_bolos = {b: t for b, t in bolos_and_times_used.items() \
                              if t > 0.5 * n_scans}
    
    logfun("$$$ Number of bolos that survived at least once: %d"
           %(len(bolos_and_times_used.keys())))
    logfun("$$$ Number of bolos that survived more than 50%% of time: %d"
           %(len(mostly_surviving_bolos.keys())))
    logfun("$$$ (i.e. more than %d scans)" %(n_scans//2))
    
    bolo_props     = cal_frame["BolometerProperties"]
    bolos_by_wafer = {wafer: [] for wafer in desired_wafers}
    
    for bolo in mostly_surviving_bolos:
        band = bolo_props[bolo].band/core.G3Units.GHz
        if not numpy.isfinite(band):
            continue
        if str(int(band)) != desired_band:
            continue
        wafer = bolo_props[bolo].wafer_id.upper()
        if wafer not in desired_wafers:
            continue
        
        bolos_by_wafer[wafer].append(bolo)
    
    return bolos_by_wafer




def collect_medians_of_pW_per_K_factors(
        calframe, bolo_dict, desired_band, logfun=print):
    
    bpm = calframe["BolometerProperties"]
    flu = core.G3Units.K * core.G3Units.rad * core.G3Units.rad
    flx = {"RCW38": { "90": 4.0549662e-07*flu,
                     "150": 2.5601153e-07*flu,
                     "220": 2.8025804e-07*flu},
           "MAT5A": { "90": 2.5738063e-07*flu,
                     "150": 1.7319235e-07*flu,
                     "220": 2.1451640e-07*flu},
           "W28A2": { "90": 4.8580000e-08*flu,
                     "150": 3.5360000e-08*flu,
                     "220": 6.5600000e-08*flu},
           "IRAS17258": { "90": 7.561e-8*flu,
                         "150": 6.161e-8*flu,
                         "220": 1.134e-7*flu},
           "RCW122A": { "90": 8.852e-8*flu,
                       "150": 7.543e-8*flu,
                       "220": 1.341e-7*flu}}
           # * Copied from spt3g_software/calibration/python/apply_t_cal.py
    
    values_by_group = \
         {group: {"PicowattsPerKelvin": [],
                  "CalibratorResponse": [],
                  "HIIFluxCalibration": [],
                  "HIIIntegratedFluux": [],
                  "HIISkyTransmission": []} \
          for group in bolo_dict.keys()}
    
    for key in calframe.keys():
        if key.endswith("FluxCalibration"):
            flux_calibration_key = key
            source = key.replace("FluxCalibration", "")
        if key.endswith("IntegralFlux"):
            integral_flux_key = key
        if key.endswith("SkyTransmission"):
            sky_transmission_key = key
    
    available_bolos = set(calframe["BolometerProperties"].keys()) & \
                      set(calframe["CalibratorResponse"].keys())  & \
                      set(calframe[flux_calibration_key].keys())  & \
                      set(calframe[integral_flux_key].keys())
    
    logfun("$$$ Recording calibration chain-related quantities ...")
    
    for group, generally_surviving_bolos in bolo_dict.items():
        if group == "AllBolos":
            continue
        good_bolos_with_data = set(generally_surviving_bolos) & \
                               set(available_bolos)
        logfun("$$$ %s : %3d detectors" %(group, len(good_bolos_with_data)))
        
        for bolo in good_bolos_with_data:
            cr = calframe["CalibratorResponse"][bolo]
            fc = calframe[flux_calibration_key][bolo]
            fi = calframe[integral_flux_key][bolo]
            if desired_band in calframe[sky_transmission_key].keys():
                st = calframe[sky_transmission_key][desired_band]
                pk = cr * fc * fi * st / flx[source][desired_band]
            else:
                st = numpy.nan
                pk = numpy.nan
            
            values_by_group[group]["PicowattsPerKelvin"].append(pk)
            values_by_group[group]["CalibratorResponse"].append(cr)
            values_by_group[group]["HIIFluxCalibration"].append(fc)
            values_by_group[group]["HIIIntegratedFluux"].append(fi)
            values_by_group[group]["HIISkyTransmission"].append(st)
            values_by_group["AllBolos"]["PicowattsPerKelvin"].append(pk)
            values_by_group["AllBolos"]["CalibratorResponse"].append(cr)
            values_by_group["AllBolos"]["HIIFluxCalibration"].append(fc)
            values_by_group["AllBolos"]["HIIIntegratedFluux"].append(fi)
            values_by_group["AllBolos"]["HIISkyTransmission"].append(st)
    
    for group in values_by_group.keys():
        for quantity in values_by_group[group].keys():
            if len(values_by_group[group][quantity]) == 0:
                values_by_group[group][quantity] = numpy.nan
            else:
                values_by_group[group][quantity] = \
                    numpy.nanmedian(values_by_group[group][quantity])
    
    logfun("$$$ Median pW/K conversion factors:")
    for group in values_by_group.keys():
        u = core.G3Units.pW / core.G3Units.K
        logfun("$$$    %s : %7.3e"
               %(group, values_by_group[group]["PicowattsPerKelvin"]/u))
    logfun("$$$ ... that's all I wanted to check.")
    
    return values_by_group




def maybe_get_gfal_copied_file(path, outdir=None):
    
    if path.startswith("/sptgrid"):
        origin = "gsiftp://osg-gridftp.grid.uchicago.edu:2811" + path
        if outdir is None:
            outdir = os.getcwd()
        gfldir = os.path.join(outdir, "temp_dir_for_gfalcp")
        if not os.path.isdir(gfldir):
            os.system("mkdir {}".format(gfldir))
        destin = os.path.join(gfldir, os.path.basename(path))
        # destin = "/sptlocal/analysis/dqtemporary/gfal_stage"
        os.system("gfal-copy {} file://{} -f -t 86400".format(origin, destin))
        path = destin
    
    return path




def maybe_del_gfal_copied_file(path):
    
    if "temp_dir_for_gfalcp" in path and os.path.isfile(path):
        os.system("rm {}".format(path))




def calculate_change_in_cal_response_vs_elevation(
        field_obs_id, dec_center, dec_height,
        cal_obs_ids, cal_ts_dir, cal_autoproc_dir,
        bolo_dict, logfun=print):
    
    response_dict = {group: {"AtBottomOfFields": [],
                             "AtMiddleOfFields": [],
                             "AtTopOfFields"   : [],
                             "FractionalChangesTopToBottom": []} \
                     for group in list(bolo_dict.keys())}
    
    def fill_with_nans():
        for wafer in response_dict.keys():
            for location in response_dict[wafer].keys():
                response_dict[wafer][location] = numpy.nan
    
    
    # - Find candidate observation IDs of the three calibrator observations
    
    cal_obs_ids = numpy.unique(cal_obs_ids)
    idx_closest_cal_obs_id = numpy.argmin(numpy.abs(cal_obs_ids - field_obs_id))
    
    try:
        id_at_bot = cal_obs_ids[idx_closest_cal_obs_id]
        id_at_mid = cal_obs_ids[idx_closest_cal_obs_id-1]
        id_at_top = cal_obs_ids[idx_closest_cal_obs_id+1]
    except:
        fill_with_nans()
        logfun("$$$ One (or more) candidate calibrator observation does not "
               "seem to have timestream files!")
        return response_dict
    
    if (id_at_bot > field_obs_id)        or \
       (id_at_mid < field_obs_id - 1200) or \
       (id_at_top > field_obs_id + 9000):
        fill_with_nans()
        logfun("$$$ One (or more) candidate observation ID does not "
               "seem to have a roughly correct number!")
        return response_dict
    
    cal_obs_ids = {"AtBottomOfFields": id_at_bot,
                   "AtMiddleOfFields": id_at_mid,
                   "AtTopOfFields"   : id_at_top}
    
    
    # - Check the elevation at which each observation was taken
    
    elevations  = {"AtBottomOfFields": 0.0,
                   "AtMiddleOfFields": 0.0,
                   "AtTopOfFields"   : 0.0}
    
    for location, cal_id in cal_obs_ids.items():
        cal_ts_file = os.path.join(cal_ts_dir, str(cal_id), "0000.g3")
        try:
            cal_ts_file = maybe_get_gfal_copied_file(cal_ts_file)
            iterator = core.G3File(cal_ts_file)
            while True:
                fr = iterator.next()
                if fr.type == core.G3FrameType.Scan:
                    mean_el = numpy.mean(fr["RawBoresightEl"])/core.G3Units.deg
                    elevations[location] = mean_el
                    break
            maybe_del_gfal_copied_file(cal_ts_file)
        except:
            pass
    
    el_center = -1 * dec_center / core.G3Units.deg
    dec_height /= core.G3Units.deg
    margin = 1.0
    half_height_ub = (dec_height / 2.0) + margin
    half_height_lb = (dec_height / 2.0) - margin
    for location, recorded_el in elevations.items():
        if location == "AtBottomOfFields":
            if (recorded_el < el_center - half_height_ub) or \
               (recorded_el > el_center - half_height_lb):
                fill_with_nans()
                logfun("$$$ The bottom elevation seems wrong!")
                return response_dict
        if location == "AtMiddleOfFields":
            if (recorded_el < el_center - margin) or \
               (recorded_el > el_center + margin):
                fill_with_nans()
                logfun("$$$ The middle elevation seems wrong!")
                return response_dict
        if location == "AtTopOfFields":
            if (recorded_el > el_center + half_height_ub) or \
               (recorded_el < el_center + half_height_lb):
                fill_with_nans()
                logfun("$$$ The top elevation seems wrong!")
                return response_dict
    
    
    # - Gather auto-processed CalibratorResponse 
    
    cal_autoproc_file_bot = \
        os.path.join(cal_autoproc_dir,
                     str(cal_obs_ids["AtBottomOfFields"])+".g3")
    cal_autoproc_file_mid = \
        os.path.join(cal_autoproc_dir,
                     str(cal_obs_ids["AtMiddleOfFields"])+".g3")
    cal_autoproc_file_top = \
        os.path.join(cal_autoproc_dir,
                     str(cal_obs_ids["AtTopOfFields"])+".g3")
    
    if (not os.path.isfile(cal_autoproc_file_bot)) or \
       (not os.path.isfile(cal_autoproc_file_mid)) or \
       (not os.path.isfile(cal_autoproc_file_top)):
        fill_with_nans()
        logfun("$$$ One or more calibration file not found!")
        return response_dict
    
    cal_autoproc_file_bot,  \
    cal_autoproc_file_mid,  \
    cal_autoproc_file_top = \
        [maybe_get_gfal_copied_file(f) for f in
         (cal_autoproc_file_bot,
          cal_autoproc_file_mid,
          cal_autoproc_file_top)]
    cal_data_bot = list(core.G3File(cal_autoproc_file_bot))[0]\
                   ["CalibratorResponse"]
    cal_data_mid = list(core.G3File(cal_autoproc_file_mid))[0]\
                   ["CalibratorResponse"]
    cal_data_top = list(core.G3File(cal_autoproc_file_top))[0]\
                   ["CalibratorResponse"]
    for f in (cal_autoproc_file_bot,
              cal_autoproc_file_mid,
              cal_autoproc_file_top):
        maybe_del_gfal_copied_file(f)
    
    
    # - Organize the data and calculate the fractional changes
    
    bolos_with_three_data = set(cal_data_bot.keys()) & \
                            set(cal_data_mid.keys()) & \
                            set(cal_data_top.keys())
    
    for group, generally_surviving_bolos in bolo_dict.items():
        if group == "AllBolos":
            continue
        good_bolos_with_data = set(generally_surviving_bolos) & \
                               set(bolos_with_three_data)
        frac_key = "FractionalChangesTopToBottom"
        for bolo in good_bolos_with_data:
            cal_resp_bot = cal_data_bot[bolo]
            cal_resp_mid = cal_data_mid[bolo]
            cal_resp_top = cal_data_top[bolo]
            try:
                frac_dif = (cal_resp_top - cal_resp_bot) / cal_resp_bot
            except:
                frac_dif = numpy.nan
            
            response_dict[group]["AtBottomOfFields"].append(cal_resp_bot)
            response_dict[group]["AtMiddleOfFields"].append(cal_resp_mid)
            response_dict[group]["AtTopOfFields"].append(cal_resp_top)
            response_dict[group][frac_key].append(frac_dif)
            response_dict["AllBolos"]["AtBottomOfFields"].append(cal_resp_bot)
            response_dict["AllBolos"]["AtMiddleOfFields"].append(cal_resp_mid)
            response_dict["AllBolos"]["AtTopOfFields"].append(cal_resp_top)
            response_dict["AllBolos"][frac_key].append(frac_dif)
    
    for group in response_dict.keys():
        for location in response_dict[group].keys():
            if len(response_dict[group][location]) == 0:
                response_dict[group][location] = numpy.nan
            else:
                response_dict[group][location] = \
                    numpy.nanmedian(response_dict[group][location])
    
    logfun("$$$ Observation IDs and the elvations:")
    for location in cal_obs_ids.keys():
        logfun("$$$    ID : %9d, Elevation : %5.2f [degree]"
               %(cal_obs_ids[location], elevations[location]))
    logfun("$$$ Median responses:")
    for key, value in response_dict["AllBolos"].items():
        if key != frac_key:
            logfun("$$$    %-16s : %4.2f [fW]" 
                   %(key, value/(1e-3*core.G3Units.pW)))
        else:
            logfun("$$$    %s : %+5.1f %%"
                   %(key, value*1e2))
    
    return response_dict




def identify_pixels_of_non_atypical_region(
        field_name, map_frame, point_source_list_file):
    
    ras, decs = maps.get_ra_dec_map(map_frame["T"])
    ras  = numpy.asarray(ras/core.G3Units.deg)
    decs = numpy.asarray(decs/core.G3Units.deg)
    
    ptsrc_msk = numpy.asarray(mapspectra.apodmask.make_apodized_ptsrc_mask(
                              map_frame, point_source_list_file,
                              apod_type="none", zero_border_arcmin=0.0))
    
    season = get_season_based_on_fields([field_name])
    if season == "winter":
        ra_max =  40.0
        ra_min = -40.0
        dccntr = float(field_name[-6:])
        de_max = dccntr + 3.00
        de_min = dccntr - 3.00
    elif season == "summer":
        ra_max = 95.0
        ra_min = 55.0
        dccntr = float(field_name[-5:])
        de_max = dccntr + 2.75
        de_min = dccntr - 2.75
    elif season == "summerp":
        ra_max = 95.0
        ra_min = 55.0
        dccntr = float(field_name[-6:])
        de_max = dccntr + 1.00
        de_min = dccntr - 1.00
    elif season == "summerb":
        ra_max = 45.0
        ra_min =  5.0
        dccntr = float(field_name[-6:])
        de_max = dccntr + 1.00
        de_min = dccntr - 1.00
    elif season == "summerc":
        ra_max = 220.0
        ra_min = 155.0
        dccntr = float(field_name[-6:])
        de_max = dccntr + 1.00
        de_min = dccntr - 1.00
    else:
        raise Exception("Not clear what season to use!")

    typical_pixels = numpy.where((ras  > ra_min) & (ras  < ra_max) &
                                 (decs > de_min) & (decs < de_max) &
                                 (ptsrc_msk == 1.0))
    
    return typical_pixels




def calculate_map_fluctuation_metrics(
        map_frame, band, sub_field, pixs=None, t_only=True, logfun=print):
    
    if t_only:
        field_map_keys = ["T"]
        weight_types   = ["TT"]
    
    k_suffix_stddev = "MapStandardDeviations"
    k_template_avgw = "MeansOf{}Weights"
    k_template_npgw = "NumbersOfPixelsWithGood{}Weights"
    
    flc_keys = []
    for fmk in field_map_keys:
        flc_keys.append(fmk+k_suffix_stddev)
    for wtp in weight_types:
        flc_keys.append(k_template_avgw.replace("{}", wtp))
        flc_keys.append(k_template_npgw.replace("{}", wtp))
    
    logfun("$$$ print some debugging info here ...")
    logfun("$$$ types of metrics to calculate:")
    for k in flc_keys:
        logfun("$$$    %s" %(k))
    
    fluctuation_metrics = {key: numpy.nan for key in flc_keys}
    
    for fmk in field_map_keys:
        assert (not map_frame[fmk].weighted), fmk + " map is still weighted!"
    
    if t_only:
        wu = 1 / (core.G3Units.mK * core.G3Units.mK)
        nominal_weights_dict = \
            { "90": {"ra0hdec-44.75": 12.0 * wu,
                     "ra0hdec-52.25": 14.0 * wu,
                     "ra0hdec-59.75": 18.0 * wu,
                     "ra0hdec-67.25": 24.0 * wu,
                     "ra5hdec-24.5": 17.0 * wu,
                     "ra5hdec-31.5": 19.0 * wu,
                     "ra5hdec-38.5": 17.0 * wu,
                     "ra5hdec-45.5": 17.0 * wu,
                     "ra5hdec-52.5": 22.0 * wu,
                     "ra5hdec-59.5": 24.0 * wu,
                     "ra5hdec-29.75": 33.0 * wu,
                     "ra5hdec-33.25": 33.0 * wu,
                     "ra5hdec-36.75": 33.0 * wu,
                     "ra5hdec-40.25": 33.0 * wu,
                     "ra1h40dec-29.75": 33.0 * wu,
                     "ra1h40dec-33.25": 33.0 * wu,
                     "ra1h40dec-36.75": 33.0 * wu,
                     "ra1h40dec-40.25": 33.0 * wu,
                     "ra12h30dec-29.75": 22.0 * wu,
                     "ra12h30dec-33.25": 22.0 * wu,
                     "ra12h30dec-36.75": 22.0 * wu,
                     "ra12h30dec-40.25": 22.0 * wu},
             "150": {"ra0hdec-44.75": 16.0 * wu,
                     "ra0hdec-52.25": 21.0 * wu,
                     "ra0hdec-59.75": 25.0 * wu,
                     "ra0hdec-67.25": 35.0 * wu,
                     "ra5hdec-24.5": 13.0 * wu,
                     "ra5hdec-31.5": 16.0 * wu,
                     "ra5hdec-38.5": 16.0 * wu,
                     "ra5hdec-45.5": 19.0 * wu,
                     "ra5hdec-52.5": 23.0 * wu,
                     "ra5hdec-59.5": 33.0 * wu,
                     "ra5hdec-29.75": 33.0 * wu,
                     "ra5hdec-33.25": 33.0 * wu,
                     "ra5hdec-36.75": 33.0 * wu,
                     "ra5hdec-40.25": 33.0 * wu,
                     "ra1h40dec-29.75": 33.0 * wu,
                     "ra1h40dec-33.25": 33.0 * wu,
                     "ra1h40dec-36.75": 33.0 * wu,
                     "ra1h40dec-40.25": 33.0 * wu,
                     "ra12h30dec-29.75": 22.0 * wu,
                     "ra12h30dec-33.25": 22.0 * wu,
                     "ra12h30dec-36.75": 22.0 * wu,
                     "ra12h30dec-40.25": 22.0 * wu},
             "220": {"ra0hdec-44.75": 1.2 * wu,
                     "ra0hdec-52.25": 1.6 * wu,
                     "ra0hdec-59.75": 2.1 * wu,
                     "ra0hdec-67.25": 2.8 * wu,
                     "ra5hdec-24.5": 1.1 * wu,
                     "ra5hdec-31.5": 1.3 * wu,
                     "ra5hdec-38.5": 1.4 * wu,
                     "ra5hdec-45.5": 1.5 * wu,
                     "ra5hdec-52.5": 1.8 * wu,
                     "ra5hdec-59.5": 2.6 * wu,
                     "ra5hdec-29.75": 2.0 * wu,
                     "ra5hdec-33.25": 2.0 * wu,
                     "ra5hdec-36.75": 2.0 * wu,
                     "ra5hdec-40.25": 2.0 * wu,
                     "ra1h40dec-29.75": 2.0 * wu,
                     "ra1h40dec-33.25": 2.0 * wu,
                     "ra1h40dec-36.75": 2.0 * wu,
                     "ra1h40dec-40.25": 2.0 * wu,
                     "ra12h30dec-29.75": 1.3 * wu,
                     "ra12h30dec-33.25": 1.3 * wu,
                     "ra12h30dec-36.75": 1.3 * wu,
                     "ra12h30dec-40.25": 1.3 * wu}}
        # * Thresholds were decided in the following ways:
        #       Winter subfields: around 0.2 x median weights from 2020 data
        #       Summer subfields: around 0.2 x median weights from 2019 data
        #       Summer-b subfields: around 0.2 x median weights from 2020 data (1st week)
        #       Summer-p subfields: same as the numbers from Summer-b subfields
        #         (Summer-p refers to the subdivided Summer el1 & el2 subfield)
        #       Summer-c subfields: 2/3 of the numbers from Summer-b subfields
        #         (Each Summer-c subfield is larger than
        #          the corresponding Summer-b subfield by 50%,
        #          but both an observation of the former and an observation of the latter
        #          last for the same amount of time.)
        
        t_vals = numpy.asarray(map_frame["T"])
        
        if "Wunpol" in map_frame.keys():
            w_vals = numpy.asarray(map_frame["Wunpol"].TT)
            pnzw = numpy.where((w_vals>0.0) & numpy.isfinite(w_vals))
            logfun("$$$ number of pixels with non-zero weights : %d"
                   %(len(pnzw[0])))
            if len(pnzw[0]) == 0:
                return fluctuation_metrics
        else:
            w_vals = None
        
        if pixs is None:
            w_cut = numpy.percentile(w_vals[pnzw], 25)
            igw = numpy.where(w_vals > w_cut)            
            good_w_vals = w_vals[igw]
            good_t_vals = t_vals[igw]
            pctl_hi = numpy.nanpercentile(good_t_vals, 99.5)
            pctl_lo = numpy.nanpercentile(good_t_vals,  0.5)
            inx = numpy.where(numpy.isfinite(good_t_vals) &
                              (good_t_vals < pctl_hi)     &
                              (good_t_vals > pctl_lo))
            better_t_vals = good_t_vals[inx]
            better_w_vals = good_w_vals[inx]
        else:
            good_t_vals = t_vals[pixs]
            if w_vals is None:
                better_t_vals = good_t_vals
                better_w_vals = None
            else:
                good_w_vals = w_vals[pixs]
                igw = numpy.where((good_w_vals>0.0) &
                                   numpy.isfinite(good_w_vals))
                better_t_vals = good_t_vals[igw]
                better_w_vals = good_w_vals[igw]
        
        map_stdddev = numpy.std(better_t_vals)
        if better_w_vals is None:
            mean_weight = numpy.nan
            n_pix_g_wgt = numpy.nan
        else:
            mean_weight = numpy.mean(better_w_vals)
            nominal_wgt = nominal_weights_dict[band][sub_field]
            n_pix_g_wgt = len(numpy.where(better_w_vals/wu > nominal_wgt)[0])
        
        logfun("$$$ number of pixels looked at : %s" %(len(pixs[0])))
        logfun("$$$ number of pixels with good weights : %s" %(n_pix_g_wgt))
        
        fluctuation_metrics["TMapStandardDeviations"] = map_stdddev
        fluctuation_metrics["MeansOfTTWeights"]       = mean_weight
        fluctuation_metrics["NumbersOfPixelsWithGoodTTWeights"] = n_pix_g_wgt
    
    logfun("$$$ ... that's all I wanted to know.")
    maps.CompactMaps(map_frame, zero_nans=True)
    
    return fluctuation_metrics




def calculate_pointing_discrepancies(
        map_frame, sub_field, map_stddev=None,
        thumbnails_file=None, thumbnails_band=None):
    
    # * For northern testing purpose only
    #   thumbnails_file = thumbnails_file.replace("bright_sources", "thumbnails")
    
    discrep_dict = {"1st": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan},
                    "2nd": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan},
                    "3rd": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan}}
    flux_dict = {"1st": numpy.nan, "2nd": numpy.nan, "3rd": numpy.nan}
    snr_dict  = {"1st": numpy.nan, "2nd": numpy.nan, "3rd": numpy.nan}

    
    # * If there are thumbnail maps of the sources
    #   made by autoprocessing, then those maps will be used.
    #   (the sources were probably interpolated over
    #    when the field map was made)
    
    if (thumbnails_file is not None) and \
       (thumbnails_band is not None) and \
       (os.path.isfile(thumbnails_file)):
        
        # * The dictionary below contains the map IDs that should be looked for
        #   when thumbnail maps are extracted from the file
        #   The IDs are the ones returned by calling
        #   sources.source_utils.get_brightest_sources(
        #       "1500d_ptsrc_150GHz_50mJy_2021.txt", field, 3)
        
        some_brightest_sources = \
            {"ra0hdec-44.75": {"1st":  "95", "2nd": "161", "3rd":  "67"},
             # * For northern testing purpose only
             #   "ra0hdec-44.75": {"1st": "829", "2nd": "401", "3rd": "135"},
             "ra0hdec-52.25": {"1st":  "53", "2nd":  "72", "3rd": "179"},
             "ra0hdec-59.75": {"1st":  "78", "2nd":  "23", "3rd":  "18"},
             "ra0hdec-67.25": {"1st": "118", "2nd":  "46", "3rd": "130"}}

        relevant_point_sources = some_brightest_sources[sub_field]
        
        for point_source_rank, point_source_id in relevant_point_sources.items():
            thumbnail_id = "{}-{}".format(point_source_id, thumbnails_band)
            iterator = core.G3File(thumbnails_file)
            while True:
                new_frame = iterator.next()
                if ("Id" in new_frame.keys()) and (new_frame["Id"] == thumbnail_id):
                    thumbnail_frame = new_frame
                    break
            deltas = pointing.astrometry.check_astrometry_at20g_templated(
                         [thumbnail_frame],
                         dist_thresh=1e9*core.G3Units.arcsec, nsigma=0.0)
            discrep_dict[point_source_rank]["DeltaRa"] = deltas["dxdec"][0]
            discrep_dict[point_source_rank]["DeltaDec"] = deltas["ddec"][0]
            snr_dict[point_source_rank] = deltas["sn"][0]
            flux_dict[point_source_rank] = numpy.nan
    
    else:
        assert (not map_frame["T"].weighted), \
               "The temperature map is still weighted!"
    
        usys = core.G3Units
    
        # * The dictionary below contains ra, dec, and flux of
        #   three brightest point sources recorded in
        #   spt3g_software/sources/1500d_ptsrc_3band_50mJy.txt as of July, 2019.
        #   The ra and dec of each source are from the ATCA20GHz catalog,
        #   but the flux is from 2018 SPT3G maps (actually, flux is not used
        #   in any calculation)."
    
        some_brightest_sources = \
            {"ra0hdec-44.75": \
                 {"1st": numpy.array([352.32358, -47.50531, 1408]),
                  "2nd": numpy.array([314.06833, -47.24664, 1375]),
                  "3rd": numpy.array([ 41.50037, -46.85467,  715])},
             "ra0hdec-52.25": \
                 {"1st": numpy.array([ 32.69288, -51.01703, 3820]),
                  "2nd": numpy.array([ 23.27404, -52.00094, 1027]),
                  "3rd": numpy.array([359.47312, -53.18686,  864])},
             "ra0hdec-59.75": \
                 {"1st": numpy.array([ 47.48363, -60.97761,  870]),
                  "2nd": numpy.array([ 45.96104, -62.19042,  833]),
                  "3rd": numpy.array([ 14.69433, -56.98650,  786])},
             "ra0hdec-67.25": \
                 {"1st": numpy.array([329.27542, -69.68981, 1115]),
                  "2nd": numpy.array([337.25092, -69.17492,  445]),
                  "3rd": numpy.array([325.44375, -64.18742,  388])},
             "ra5hdec-24.5": \
                 {"1st": numpy.array([ 74.26346, -23.41439, 3841]),
                  "2nd": numpy.array([ 52.47542, -23.95239, 1458]),
                  "3rd": numpy.array([ 92.24900, -22.33922, 1003])},
             "ra5hdec-31.5": \
                 {"1st": numpy.array([ 69.40233, -29.90108,  564]),
                  "2nd": numpy.array([ 61.89133, -33.06258,  563]),
                  "3rd": numpy.array([ 60.58879, -31.79044,  467])},
             "ra5hdec-38.5": \
                 {"1st": numpy.array([ 60.97404, -36.08358, 4013]),
                  "2nd": numpy.array([ 80.74142, -36.45844, 3909]),
                  "3rd": numpy.array([ 67.16821, -37.93867, 1850])},
             "ra5hdec-45.5": \
                 {"1st": numpy.array([ 79.95708, -45.77881, 6320]),
                  "2nd": numpy.array([ 84.70979, -44.08572, 5294]),
                  "3rd": numpy.array([ 73.96163, -46.26628, 4157])},
             "ra5hdec-52.5": \
                 {"1st": numpy.array([ 85.19104, -54.30606, 1127]),
                  "2nd": numpy.array([ 57.86771, -51.71525,  401]),
                  "3rd": numpy.array([ 92.20471, -54.94525,  386])},
             "ra5hdec-59.5": \
                 {"1st": numpy.array([ 76.68342, -61.16150, 1655]),
                  "2nd": numpy.array([ 87.53988, -57.54017, 1005]),
                  "3rd": numpy.array([ 80.64333, -61.13250,  573])},
             "ra5hdec-29.75": \
                 {"1st": numpy.array([ 84.97571, -28.66561,  677]),
                  "2nd": numpy.array([ 69.40233, -29.90108,  564]),
                  "3rd": numpy.array([ 95.12212, -28.46006,  503])},
             "ra5hdec-33.25": \
                 {"1st": numpy.array([ 84.11850, -34.01967,  842]),
                  "2nd": numpy.array([ 61.89133, -33.06258,  563]),
                  "3rd": numpy.array([ 60.58879, -31.79044,  467])},
             "ra5hdec-36.75": \
                 {"1st": numpy.array([ 60.97404, -36.08358, 4013]),
                  "2nd": numpy.array([ 80.74142, -36.45844, 3909]),
                  "3rd": numpy.array([ 67.16821, -37.93867, 1850])},
             "ra5hdec-40.25": \
                 {"1st": numpy.array([ 53.55675, -40.14036, 1274]),
                  "2nd": numpy.array([ 90.13046, -39.61711,  586]),
                  "3rd": numpy.array([ 72.42625, -39.18600,  401])},
             "ra1h40dec-29.75": \
                 {"1st": numpy.array([  2.64967, -30.46339,  741]),
                  "2nd": numpy.array([ 39.12962, -29.89864,  612]),
                  "3rd": numpy.array([ 18.94342, -30.82197,  472])},
             "ra1h40dec-33.25": \
                 {"1st": numpy.array([ 35.73500, -34.69103, 1019]),
                  "2nd": numpy.array([ 28.29246, -33.17408,  542]),
                  "3rd": numpy.array([ 34.20079, -32.79458,  540])},
             "ra1h40dec-36.75": \
                 {"1st": numpy.array([  6.56833, -35.21369, 1123]),
                  # "2nd": numpy.array([ 23.49125, -36.49314,  440]),
                  # * The measured offsets are too large for this source for some reason...
                  "2nd": numpy.array([ 42.73092, -36.27647,  340]),
                  "3rd": numpy.array([ 37.36842, -36.73242,  263])},
             "ra1h40dec-40.25": \
                 {"1st": numpy.array([ 16.68796, -40.57208, 2146]),
                  "2nd": numpy.array([  3.24954, -39.90733, 1609]),
                  "3rd": numpy.array([ 23.63392, -38.72603,  685])},
             "ra12h30dec-29.75": \
                 {"1st": numpy.array([159.31675, -29.56744, 2681]),
                  # "2nd": numpy.array([205.57287, -29.01333,  454]),
                  "2nd": numpy.array([154.61983, -31.39811,  562]),
                  "3rd": numpy.array([176.60854, -28.98883,  425])},
             "ra12h30dec-33.25": \
                 {# "1st": numpy.array([199.03371, -33.64967, 1497]),
                  # "2nd": numpy.array([194.49667, -31.92089, 1460]),
                  # "3rd": numpy.array([216.92212, -33.09219,  795])},
                  "1st": numpy.array([151.88067, -33.55183,  652]),
                  "2nd": numpy.array([156.00204, -32.57092,  491]),
                  "3rd": numpy.array([162.77021, -31.63742,  439])},
             "ra12h30dec-36.75": \
                 {# "1st": numpy.array([223.61446, -37.79250, 1886]),
                  "1st": numpy.array([176.75608, -38.20297, 1376]),
                  "2nd": numpy.array([178.59079, -35.09142,  889]),
                  "3rd": numpy.array([159.22262, -37.73753,  812])},
                  # "3rd": numpy.array([196.30371, -36.92919,  239])},
                  # * Too large of offsets (RA)...
             "ra12h30dec-40.25": \
                 {# "1st": numpy.array([197.45250, -39.80878,  416]),
                  # * The measured offsets are too large for this source...
                  #   Actually, for a few sources, the RA offsets are
                  #   on the order of 1e5 arcseconds...
                  "1st": numpy.array([162.15954, -41.23331,  369]),
                  "2nd": numpy.array([166.29617, -39.47856,  285]),
                  "3rd": numpy.array([160.68488, -41.72386,  233])}}
            
        
        relevant_point_sources = some_brightest_sources[sub_field]
        
        for point_source_rank, info in relevant_point_sources.items():
            
            true_right_ascension = info[0]
            if true_right_ascension > 180:
                true_right_ascension -= 360
            true_right_ascension *= usys.deg
            true_declination  = info[1]
            true_declination *= usys.deg
            true_right_ascension = float(true_right_ascension)
            true_declination = float(true_declination)
            
            true_x, true_y = \
                map_frame["T"].angle_to_xy(true_right_ascension,
                                           true_declination)
            nominal_pixel_number = \
                map_frame["T"].angle_to_pixel(true_right_ascension,
                                              true_declination)
            nominal_y_index, nominal_x_index = \
                numpy.unravel_index(
                    nominal_pixel_number, map_frame["T"].shape)
            
            nominal_x_index = int(nominal_x_index)
            nominal_y_index = int(nominal_y_index)
            
            mini_map_width    = 4 * usys.arcmin
            mini_map_height   = 4 * usys.arcmin
            mini_map_x_length = int(mini_map_width  / map_frame["T"].x_res)
            mini_map_y_length = int(mini_map_height / map_frame["T"].y_res)
            half_x_length     = mini_map_x_length // 2
            half_y_length     = mini_map_y_length // 2
            mini_map_x_left   = nominal_x_index - half_x_length
            mini_map_x_right  = nominal_x_index + half_x_length
            mini_map_y_bottom = nominal_y_index - half_y_length
            mini_map_y_top    = nominal_y_index + half_y_length
            mini_map = numpy.asarray(map_frame["T"])\
                           [mini_map_y_bottom:mini_map_y_top+1,
                            mini_map_x_left:mini_map_x_right+1]
            
            try:
                fit_parameters = util.fitting.fit_gaussian2d(mini_map)
            except ValueError:
                discrep_dict[point_source_rank]["DeltaRa"]  = numpy.nan
                discrep_dict[point_source_rank]["DeltaDec"] = numpy.nan
                continue
            
            gaussian_center_x = fit_parameters[1]
            gaussian_center_y = fit_parameters[2]
            
            measured_x = gaussian_center_x + mini_map_x_left
            measured_y = gaussian_center_y + mini_map_y_bottom
            measured_x = float(measured_x)
            measured_y = float(measured_y)
            
            measured_right_ascension, measured_declination = \
                map_frame["T"].xy_to_angle(measured_x, measured_y)
            
            delta_ra  = measured_right_ascension - true_right_ascension
            delta_dec = measured_declination     - true_declination
            effective_delta_ra = delta_ra * numpy.cos(true_declination/usys.rad)
            
            discrep_dict[point_source_rank]["DeltaRa"]  = effective_delta_ra
            discrep_dict[point_source_rank]["DeltaDec"] = delta_dec
            
            """
            nearest_x = int(numpy.round(measured_x))
            nearest_y = int(numpy.round(measured_y))
            source_centered_mini_map = \
                numpy.asarray(map_frame["T"])\
                    [nearest_y-half_y_length:nearest_y+half_y_length+1,
                     nearest_x-half_x_length:nearest_x+half_x_length+1]
            """
            source_centered_mini_map = mini_map
            # * If the measured positions are crazy,
            #   the source-centered mini map can be empty...
            source_flux = numpy.sum(source_centered_mini_map) * \
                          map_frame["T"].x_res * map_frame["T"].y_res
            flux_dict[point_source_rank] = source_flux
            
            if map_stddev is None:
                source_snr = numpy.nan
            else:
                source_snr = numpy.max(source_centered_mini_map) / map_stddev
            snr_dict[point_source_rank] = source_snr
        
            maps.CompactMaps(map_frame, zero_nans=True)
    
    return discrep_dict, flux_dict, snr_dict




def create_mask_for_powspec_calc_of_small_region(
        map_frame, point_source_file,
        center_ra, center_dec, height, margin=0.25*core.G3Units.deg):
    
    """
    ptsrc_mask = mapspectra.apodmask.make_apodized_ptsrc_mask(
                     map_frame, point_source_file)
    
    edge_mask = maps.FlatSkyMap(
                    x_len=map_frame["T"].shape[1],
                    y_len=map_frame["T"].shape[0],
                    res=map_frame["T"].res, proj=map_frame["T"].proj,
                    coord_ref=map_frame["T"].coord_ref,
                    alpha_center=map_frame["T"].alpha_center,
                    delta_center=map_frame["T"].delta_center,
                    pol_type=map_frame["T"].pol_type, units=None)
    
    twod_window  = numpy.ones(edge_mask.shape)
    twod_window *= signal.tukey(edge_mask.shape[1], alpha=0.2)
    twod_window  = (twod_window.transpose() *\
                    signal.tukey(edge_mask.shape[0], alpha=0.2)).\
                   transpose()
    numpy.asarray(edge_mask)[:,:] = twod_window
    
    full_mask = ptsrc_mask * edge_mask
    
    return full_mask
    """
    
    center_ra /= core.G3Units.deg
    center_dec /= core.G3Units.deg
    height /= 2.0*core.G3Units.deg
    
    source_mask = mapspectra.apodmask.make_apodized_ptsrc_mask(
                      map_frame, point_source_file)
    
    edge_mask = source_mask.clone(False)
    ras, decs = maps.get_ra_dec_map(source_mask)
    ras  = numpy.asarray(ras/core.G3Units.deg)
    decs = numpy.asarray(decs/core.G3Units.deg)
    numpy.asarray(edge_mask)[(decs>center_dec-height) &
                             (decs<center_dec+height)] = 1.0    
    del ras, decs
    
    margin_npix = int(margin/edge_mask.res)
    numpy.asarray(edge_mask)[:margin_npix,:] = 0.0
    numpy.asarray(edge_mask)[-margin_npix:,:] = 0.0
    numpy.asarray(edge_mask)[:,:margin_npix] = 0.0
    numpy.asarray(edge_mask)[:,-margin_npix:] = 0.0
    
    edge_mask = mapspectra.apodmask.make_border_apodization(
                edge_mask, smooth_weights_arcmin=0.0,
                apod_type="cosine", radius_arcmin=60)
    
    full_mask = source_mask * edge_mask
    
    return full_mask




def create_mini_planck_map_frame(
        fits_file, spt_map_frame,
        center_ra, center_dec, length_x, length_y,
        t_only=True):
    
    planck_healpix_maps = maps.fitsio.load_skymap_fits(fits_file)
    
    if t_only:
        planck_flatsky_map = maps.maputils.healpix_to_flatsky(
                                 planck_healpix_maps["T"],
                                 res=spt_map_frame["T"].res,
                                 x_len=spt_map_frame["T"].shape[1],
                                 y_len=spt_map_frame["T"].shape[0],
                                 alpha_center=spt_map_frame["T"].alpha_center,
                                 delta_center=spt_map_frame["T"].delta_center,
                                 proj=spt_map_frame["T"].proj,
                                 coord_ref=spt_map_frame["T"].coord_ref,
                                 pol_type=spt_map_frame["T"].pol_type)
        new_map_frame = core.G3Frame(core.G3FrameType.Map)
        new_map_frame["T"] = planck_flatsky_map
    
        map_frame_smaller = create_new_map_frame_with_smaller_region(
                                new_map_frame,
                                center_ra, center_dec, length_x, length_y)

    return map_frame_smaller




def calculate_average_ratio_of_spt_planck_xspectra(
        spt_map_frame, planck_map_frames, mask, t_only=True):
    
    assert (not spt_map_frame["T"].weighted), \
           "The SPT T map is still weighted!"
    
    for mission, map_frame in planck_map_frames.items():
        if mission == "fullmission":
            fullmission_planck_map_frame  = map_frame
        if mission == "halfmission-1":
            halfmission1_planck_map_frame = map_frame
        if mission == "halfmission-2":
            halfmission2_planck_map_frame = map_frame
        
    average_ratios = {}
    
    ell_lo = 790
    ell_hi = 1410
    
    if t_only:
        if not numpy.all(numpy.isfinite(numpy.asarray(halfmission1_planck_map_frame["T"]))):
            halfmission1_planck_map_frame = fullmission_planck_map_frame
            # * The 217/220 GHz mini Planck maps for the Summer-c field
            #   has NaN values for some reason...
        
        planck_x_planck = mapspectra.map_analysis.calculate_powerspectra(
                              halfmission1_planck_map_frame,
                              input2=halfmission2_planck_map_frame,
                              delta_l=200, lmin=700, lmax=1500,
                              apod_mask=mask, realimag="real", flatten=False)
        spt_x_planck = mapspectra.map_analysis.calculate_powerspectra(
                           spt_map_frame, input2=fullmission_planck_map_frame,
                           delta_l=200, lmin=700, lmax=1500,
                           apod_mask=mask, realimag="real", flatten=False)
        
        cl_ratios   = spt_x_planck["TT"] / planck_x_planck["TT"]
        bin_centers = planck_x_planck["TT"].bin_centers
        idx = numpy.where((bin_centers > ell_lo) & (bin_centers < ell_hi))[0]
        average_ratio_tt = numpy.mean(cl_ratios[idx])
        
        average_ratios["TT"] = average_ratio_tt
    
    return average_ratios




def calculate_noise_levels(map_frame, mask, t_only=True):
    
    assert (not map_frame["T"].weighted), \
            "The T maps is still weighted!"
    
    noises = {}
    
    ell_lo = 2990
    ell_hi = 5010
    
    if t_only:
        cls = mapspectra.map_analysis.calculate_powerspectra(
                  map_frame, input2=None,
                  delta_l=200, lmin=2500, lmax=5500,
                  apod_mask=mask, realimag="real", flatten=False)
        
        bin_centers = cls["TT"].bin_centers
        idx = numpy.where((bin_centers > ell_lo) & (bin_centers < ell_hi))[0]
        noise_t = numpy.sqrt(numpy.mean(cls["TT"][idx]))
        
        noises["T"] = noise_t
    
    return noises


# ==============================================================================




# ==============================================================================
# Define the pipeline module that calls the functions defined above to
# add maps together and perform the anlyses
# ------------------------------------------------------------------------------


class AnalyzeAndCoaddMaps(object):
    
    def __init__(self,
                 map_ids=["90GHz"], map_sources=["ra0hdec-44.75"], t_only=True,
                 maps_split_by_scan_direction=False,
                 combine_left_right=False,
                 combine_different_wafers=False,
                 allow_subtraction=False,
                 collect_averages_from_flagging_stats=False,
                 calculate_pW_to_K_conversion_factors=False,
                 calculate_calibrator_response_vs_el=False,
                 calibration_data_dir=None,
                 bolo_timestreams_dir=None,
                 min_field_obs_id=None,
                 max_field_obs_id=None,
                 calculate_map_rmss_and_weight_stats=False,
                 rmss_and_wgts_from_coadds_or_individuals=None,
                 rmss_and_wgts_from_signals_or_noises=None,
                 calculate_pointing_discrepancies=False,
                 calculate_noise_from_individual_maps=False,
                 calculate_noise_from_coadded_maps=False,
                 point_source_list_file=None,
                 calculate_cross_spectra_with_planck_map=False,
                 planck_map_fits_file=None,
                 auxiliary_files_directory=None,
                 logging_function=logging.info,
                 less_verbose=False):
        
        self.log = logging_function
        if less_verbose:
            def dummy(x):
                return
            self.detail = dummy
        else:
            self.detail = self.log
        
        # - Initialize variables related to map IDs

        self.combine_left_right = combine_left_right
        if self.combine_left_right:
            directions = ["Left", "Right"]
        else:
            directions = [""]
        
        self.combine_different_wafers = combine_different_wafers
        if self.combine_different_wafers:
            wafers = ["W172", "W174", "W176", "W177", "W180",
                      "W181", "W188", "W203", "W204", "W206"]
        else:
            wafers = [""]
        
        self.map_ids = map_ids
        self.id_mapping = {}
        for map_id in self.map_ids:
            for direction in directions:
                for wafer in wafers:
                    self.id_mapping[direction+map_id+wafer] = map_id
        
        if self.combine_left_right or self.combine_different_wafers:
            self.log("- Since maps from different directions/wafers")
            self.log("  are to be combined:")
            self.log("  an ID stored in a map frame will be regarded")
            self.log("  as a different ID according to the relations below:")
            for id_from in sorted(self.id_mapping.keys()):
                id_to = self.id_mapping[id_from]
                self.log("      %16s ==> %16s", id_from, id_to)
        
        
        # - Initialize variables related to storing coadded maps and
        #   their basic information
        
        self.map_sources = map_sources
        self.t_only      = t_only
        
        self.coadded_map_frames = \
            {mid: core.G3Frame(core.G3FrameType.Map) for mid in self.map_ids}
        self.coadded_map_ids = \
            {mid: core.G3MapVectorString() for mid in self.map_ids}
        self.coadded_obs_ids = \
            {mid: core.G3MapVectorInt() for mid in self.map_ids}
        self.ignored_obs_ids = \
            {mid: core.G3MapVectorInt() for mid in self.map_ids}
        self.obser_durations = \
            {mid: core.G3MapMapDouble() for mid in self.map_ids}
        self.allow_subtraction = allow_subtraction
        self.obs_info = None
        
        
        # - Initialize variables related to calculating flagging stats,
        #   calculating pW/K numbers, and calibrator response changes
        
        self.calc_median_pW_to_K = calculate_pW_to_K_conversion_factors
        self.calc_cal_vs_el      = calculate_calibrator_response_vs_el
        self.get_avgs_flagging_stats = collect_averages_from_flagging_stats
        
        self.collect_by_band_wafer_quantities = self.calc_median_pW_to_K or \
                                                self.calc_cal_vs_el      or \
                                                self.get_avgs_flagging_stats
        
        if self.collect_by_band_wafer_quantities:
            self.wafers = ["W172", "W174", "W176", "W177", "W180",
                           "W181", "W188", "W203", "W204", "W206",
                           "AllBolos"]
            self.pipe_info_frame  = None
            self.cal_frame        = None
            self.wafers_and_bolos = None
        
        
        # - Initialize variables related to collecting average numbers
        #   of various bolometer flagging statistics
        
        if self.get_avgs_flagging_stats:
            self.flagging_reasons = \
                ["BadCalSn", "BadWeight",  "Glitchy", "Oscillating",
                 "Latched",  "Overbiased", "BadHk",
                 "PostCalibrationNaNs",    "UnphysicalLowVariance",
                 "MissingFluxCalibration",
                 "Others",  "TotalNotFlagged", "TotalRemoved"]
            self.key_prefix_flg = "FlaggingStatisticsAverageNumbersOf"
            
            self.avgs_flagging_stats =  \
                {wafer: {flag_rea: {map_id: core.G3MapMapDouble() \
                                    for map_id in self.map_ids}   \
                         for flag_rea in self.flagging_reasons}   \
                 for wafer in self.wafers}
        
        
        # - Initialize variables related to calculating the median value of
        #   the pW/K conversion factors for each wafer
        
        if self.calc_median_pW_to_K:
            self.calfactors = ["PicowattsPerKelvin", "CalibratorResponse",
                               "HIIFluxCalibration", "HIIIntegratedFluux",
                               "HIISkyTransmission"]
            usys = core.G3Units
            self.cal_fat_unis = \
                {"PicowattsPerKelvin": usys.pW / usys.K,
                 "CalibratorResponse": usys.pW / 1000,
                 "HIIFluxCalibration": 1.0,
                 "HIIIntegratedFluux": usys.arcmin * usys.arcmin,
                 "HIISkyTransmission": 1.0}
            self.cal_fat_ustr = \
                {"PicowattsPerKelvin": "pW/K",
                 "CalibratorResponse": "fW",
                 "HIIFluxCalibration": "unitless",
                 "HIIIntegratedFluux": "arcmin^2",
                 "HIISkyTransmission": "unitless"}
            self.key_prefix_pwk = "MediansOfTemperatureCalRelatedFactors"
            
            self.medians_of_temp_cal_factors = \
                {wafer: {factor: {map_id: core.G3MapMapDouble() \
                                  for map_id in self.map_ids}   \
                         for factor in self.calfactors}         \
                 for wafer in self.wafers}
        
        
        # - Initialize variables related to calculating changeds in detectors'
        #   response to the calibrator when the telescope elevation was at
        #   the bottom and the top of the field.
        
        if self.calc_cal_vs_el:
            self.cal_ts_dir  = \
                os.path.join(bolo_timestreams_dir, "calibrator")
            self.cal_autoproc_dir = \
                os.path.join(calibration_data_dir, "calibrator")
            self.cal_obs_ids = \
                [int(obs_id) for obs_id in os.listdir(self.cal_ts_dir) \
                 if (int(obs_id) >= min_field_obs_id-1200) and         \
                    (int(obs_id) <= max_field_obs_id+9000)]
            
            self.cal_locations = ["AtBottomOfFields",
                                  "AtMiddleOfFields",
                                  "AtTopOfFields",
                                  "FractionalChangesTopToBottom"]
            self.key_prefix_calel = "MedianCalibratorResponse"
            
            self.cal_resp_diffs_els = \
                {wafer: {location: {map_id: core.G3MapMapDouble() \
                                    for map_id in self.map_ids}   \
                         for location in self.cal_locations}      \
                 for wafer in self.wafers}
        
        
        # - Initialize variables related to calculating a Std. dev. of
        #   a map, a mean TT weight, and number of pixels in the weight map
        #   whose weights are above some thresholds.
        
        self.calc_map_stddev_etc = calculate_map_rmss_and_weight_stats
        
        if self.calc_map_stddev_etc:
            if self.t_only:
                self.fluctuation_keys = ["TMapStandardDeviations",
                                         "MeansOfTTWeights",
                                         "NumbersOfPixelsWithGoodTTWeights"]
                tu = core.G3Units.mK
                self.flu_met_unis = {"TMapStandardDeviations": tu,
                                     "MeansOfTTWeights"      : 1/(tu*tu),
                                     "NumbersOfPixelsWithGoodTTWeights": 1}
                self.flu_met_ustr = {"TMapStandardDeviations": "mK",
                                     "MeansOfTTWeights"      : "1/mK^2",
                                     "NumbersOfPixelsWithGoodTTWeights": ""}
            
            self.coadds_or_individs_for_flc = []
            if "c" in rmss_and_wgts_from_coadds_or_individuals:
                self.coadds_or_individs_for_flc.append("Coadded")
            if "i" in rmss_and_wgts_from_coadds_or_individuals:
                self.coadds_or_individs_for_flc.append("Individual")
            
            self.signals_or_noises_for_flc = []
            if "s" in rmss_and_wgts_from_signals_or_noises:
                self.signals_or_noises_for_flc.append("Signal")
            if "n" in rmss_and_wgts_from_signals_or_noises:
                self.signals_or_noises_for_flc.append("Noise")
            
            self.map_types_for_flc = []
            for c_or_i in self.coadds_or_individs_for_flc:
                for s_or_n in self.signals_or_noises_for_flc:
                    self.map_types_for_flc.append(c_or_i+s_or_n+"Maps")
            
            self.key_prefix_flc = "FluctuationMetrics"
            
            self.map_fluctuation_metrics = \
                {map_type: {fluc_key: {map_id: core.G3MapMapDouble()  \
                                       for map_id in self.map_ids}    \
                            for fluc_key in self.fluctuation_keys}    \
                 for map_type in self.map_types_for_flc}
            
            self.pixels_to_use_for_flc_calc = \
                {sub_field: None for sub_field in self.map_sources}
        
        
        # - Initialize variables related to pointing dicrepancy calculations
        
        self.calc_pointing_discrepancies = calculate_pointing_discrepancies
        
        if self.calc_pointing_discrepancies:
            self.ptsrc_ranks        = ["1st", "2nd", "3rd"]
            self.ptsrc_offset_types = ["DeltaRa", "DeltaDec"]
            key_suffix_ptsrc = "OfNthBrightestSourceFromEachSubfield"
            self.key_template_ptsrc_d = "Deltas" + key_suffix_ptsrc
            self.key_template_ptsrc_f = "Fluxes" + key_suffix_ptsrc
            self.key_template_ptsrc_s = "SNRs"   + key_suffix_ptsrc
            
            self.ptsrc_delta_ras_and_decs = \
                {rank: {delta_type: {map_id: core.G3MapMapDouble()  \
                                     for map_id in self.map_ids}    \
                        for delta_type in self.ptsrc_offset_types}  \
                 for rank in self.ptsrc_ranks}
            self.ptsrc_fluxes = \
                {rank: {map_id: core.G3MapMapDouble()  \
                        for map_id in self.map_ids}    \
                 for rank in self.ptsrc_ranks}
            self.ptsrc_snrs = \
                {rank: {map_id: core.G3MapMapDouble()  \
                        for map_id in self.map_ids}    \
                 for rank in self.ptsrc_ranks}
            
            self.field_autoproc_dir = os.path.join(calibration_data_dir, "FIELD")
        
        
        # - Initialize variables related to calculations of cross spectra
        #   with Planck maps
        
        self.calc_xspec_with_plck_map = calculate_cross_spectra_with_planck_map
        
        if self.calc_xspec_with_plck_map:
            self.planck_map_fits_file = planck_map_fits_file
            self.spt_to_planck_bands  = \
                {"90": "100", "150": "143", "220": "217"}
            self.planck_missions = \
                ["fullmission", "halfmission-1", "halfmission-2"]
            self.mini_planck_map_frames = \
                {mission: {band: {sub_field: None                       \
                                  for sub_field in self.map_sources}    \
                           for band in self.spt_to_planck_bands.keys()} \
                 for mission in self.planck_missions}
            
            self.av_xspec_ratio_key = "AveragesOfRatiosOfSPTxPlancktoPlckxPlck"
            
            if self.t_only:
                self.xspec_types = ["TT"]
            
            self.avgs_xspec_ratios = \
                {xspec_type: {map_id: core.G3MapMapDouble() \
                              for map_id in self.map_ids}   \
                 for xspec_type in self.xspec_types}
        
        
        # - Initialize variables related to noise calculations
        
        self.calc_noise_from_individ_maps = calculate_noise_from_individual_maps
        self.calc_noise_from_coadded_maps = calculate_noise_from_coadded_maps
        
        if self.t_only:
            self.map_types_for_noise_calc = ["T"]
        
        if self.calc_noise_from_individ_maps:
            self.indvid_noise_key = "NoiseLevelsFromIndividual{}Maps"
            self.noise_from_individual_maps = \
                {map_type: {map_id: core.G3MapMapDouble() \
                            for map_id in self.map_ids}   \
                 for map_type in self.map_types_for_noise_calc}
        if self.calc_noise_from_coadded_maps:
            self.coadds_noise_key = "NoiseLevelsFromCoadded{}Maps"
            self.noise_from_coadded_maps = \
                {map_type: {map_id: core.G3MapMapDouble() \
                            for map_id in self.map_ids}   \
                 for map_type in self.map_types_for_noise_calc}
        
        
        # - Initialize other variables needed for some of the four analyses,
        #   i.e., calculating basic stats of fluctuations in maps,
        #   cross spectra with Planck maps, noise levels, and pointing offsets.        
        
        self.analyze_maps = self.calc_map_stddev_etc          or \
                            self.calc_pointing_discrepancies  or \
                            self.calc_xspec_with_plck_map     or \
                            self.calc_noise_from_individ_maps or \
                            self.calc_noise_from_coadded_maps
        
        if self.analyze_maps:    
            self.maps_split_by_scan_direction = maps_split_by_scan_direction
            self.point_source_list_file       = point_source_list_file
            self.masks_for_powspec_calculations = \
                {sub_field: None for sub_field in self.map_sources}
            self.aux_files_dir = auxiliary_files_directory
            
            if self.calc_map_stddev_etc:
                for sub_field in self.pixels_to_use_for_flc_calc.keys():
                    pf = "pixel_numbers_for_calculating_"+\
                         "fluctuation_metrics_of_{}_sub_field.pickle".\
                         format(sub_field)
                    pf = os.path.join(self.aux_files_dir, pf)
                    if os.path.isfile(pf):
                        with open(pf, "rb") as f_obj:
                            indices = pickle.load(f_obj)
                            self.pixels_to_use_for_flc_calc[sub_field] = indices
            
            if self.calc_xspec_with_plck_map:
                for band in self.spt_to_planck_bands.keys():
                    for mission in self.planck_missions:
                        for sub_field in self.map_sources:
                            mf = ("mini_{}_planck_map_for_spt_"+\
                                  "{}GHz_{}_sub_field.g3").\
                                 format(mission, band, sub_field)
                            mf = os.path.join(self.aux_files_dir, mf)
                            if os.path.isfile(mf):
                                fr = list(core.G3File(mf))[0]
                                self.mini_planck_map_frames\
                                    [mission][band][sub_field] = fr
            
            if self.calc_xspec_with_plck_map     or \
               self.calc_noise_from_individ_maps or \
               self.calc_noise_from_coadded_maps:
                for sub_field in self.map_sources:
                    mf = "mask_for_power_spectrum_calculations_"+\
                         "for_{}_sub_field.g3".format(sub_field)
                    mf = os.path.join(self.aux_files_dir, mf)
                    if os.path.isfile(mf):
                        mask = list(core.G3File(mf))[0]["Mask"]
                        self.masks_for_powspec_calculations[sub_field] = mask
            
            if self.maps_split_by_scan_direction:
                self.direction_independent_map_ids = []
                self.map_fr_arrival_counters       = {}
                self.map_cache_for_both_dirs       = {}
                for map_id in self.map_ids:
                    if "Left" in map_id:
                        self.direction_independent_map_ids.append(
                            map_id.replace("Left", "SomeDirection"))
                        self.map_fr_arrival_counters[map_id] = \
                            {sub_field: 0 for sub_field in self.map_sources}
                        self.map_cache_for_both_dirs[map_id] = \
                            {sub_field: None for sub_field in self.map_sources}
                    elif "Right" in map_id:
                        self.direction_independent_map_ids.append(
                            map_id.replace("Right", "SomeDirection"))
                        self.map_fr_arrival_counters[map_id] = \
                            {sub_field: 0 for sub_field in self.map_sources}
                        self.map_cache_for_both_dirs[map_id] = \
                            {sub_field: None for sub_field in self.map_sources}
                    elif self.combine_left_right:
                        self.direction_independent_map_ids.append(map_id)
                        self.map_fr_arrival_counters["Left"+map_id]  = \
                            {sub_field: 0 for sub_field in self.map_sources}
                        self.map_fr_arrival_counters["Right"+map_id] = \
                            {sub_field: 0 for sub_field in self.map_sources}
                        self.map_cache_for_both_dirs["Left"+map_id]  = \
                            {sub_field: None for sub_field in self.map_sources}
                        self.map_cache_for_both_dirs["Right"+map_id] = \
                            {sub_field: None for sub_field in self.map_sources}
                self.direction_independent_map_ids = \
                    set(self.direction_independent_map_ids)
    
    
    
    def combine_ids(self, old_dict, new_dict, report_commonness=False):
        
        common_ids_present = False
        for map_id in new_dict.keys():
            for sub_field, new_ids in new_dict[map_id].items():
                if sub_field not in old_dict[map_id].keys():
                    old_dict[map_id][sub_field] = new_ids
                else:
                    old_ids    = old_dict[map_id][sub_field]
                    common_ids = set(new_ids) & set(old_ids)
                    if len(common_ids) > 0:
                        common_ids_present = True
                    else:
                        for new_id in new_ids:
                            old_dict[map_id][sub_field].append(new_id)
        
        if report_commonness:
            return old_dict, common_ids_present
        else:
            return old_dict
    
    
    def remove_ids(self, old_dict, new_dict):
        
        for map_id in new_dict.keys():
            for sub_field, ids_to_remove in new_dict[map_id].items():
                updated_ids = []
                for old_id in old_dict[map_id][sub_field]:
                    if old_id not in ids_to_remove:
                        updated_ids.append(old_id)
                updated_ids = type(old_dict[map_id][sub_field])(updated_ids)
                old_dict[map_id].pop(sub_field)
                old_dict[map_id][sub_field] = updated_ids
        
        return old_dict
    
    
    def create_mmd_for_one_value(self, sub_field, obs_id, value):
        
        mmd = core.G3MapMapDouble()
        mmd[sub_field] = core.G3MapDouble()
        mmd[sub_field][str(obs_id)] = value
        
        return mmd
    
    
    def combine_mapmapdoubles(self, old_dict, new_dict):
        
        for map_id in new_dict.keys():
            for sub_field, some_map in new_dict[map_id].items():
                if sub_field not in old_dict[map_id].keys():
                    old_dict[map_id][sub_field] = some_map
                else:
                    for obs_id, some_value in some_map.items():
                        if obs_id not in old_dict[map_id][sub_field].keys():
                            old_dict[map_id][sub_field][obs_id] = some_value
        
        return old_dict
    
    
    def remove_partial_mapmapdouble(self, old_dict, new_dict):
        
        for map_id in new_dict.keys():
            for sub_field, obs_ids_to_remove in new_dict[map_id].items():
                for obs_id in obs_ids_to_remove:
                    old_dict[map_id][sub_field].pop(str(obs_id))
        
        return old_dict
    
    
    def print_mapmapdouble(self, mapmapdouble, u, initial_indent):
        
        indlev = initial_indent
        self.log("")
        for sub_field, data in mapmapdouble.items():
            self.log(" "*indlev + "- %s", sub_field)
            self.log("")
            oids = data.keys()
            for oid, datum in data.items():
                self.log(" "*(indlev+3) + "%12s : %+7.3e", oid, datum/u)
            self.log("")
    
    
    def map_seems_fine(self, map_frame):
        
        if "IgnoreThisMap" in map_frame.keys():
            self.log("")
            self.log("* This map was determined to be bad from ")
            self.log("* a previous processing step!")
            self.log("")
            return False
        
        if map_frame["T"].units != core.G3TimestreamUnits.Tcmb:
            self.log("")
            self.log("* The units of this map are not in Tcmb!")
            self.log("")
            return False
        
        mean = numpy.nanmean(numpy.asarray(map_frame["T"]))
        maps.CompactMaps(map_frame, zero_nans=True)
        if not numpy.isfinite(mean):
            self.log("")
            self.log("* There seem to be only NaNs in the map!")
            self.log("")
            return False
        
        return True
    
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Observation:
            current_time = datetime.datetime.utcnow()
            self.log("")
            self.log("-----------------------------------------------------")
            self.log("Found an observation frame!")
            self.log("(Probably a new g3 file has arrived to the pipeline.)")
            self.log("(Current UTC time is %s.)", current_time)
            self.log("-----------------------------------------------------")
            self.log("")
            self.obs_info = frame
            return []
        
        if frame.type == core.G3FrameType.Calibration:
            if self.collect_by_band_wafer_quantities:
                if "BolometerProperties" in frame.keys():
                    self.cal_frame = frame
            return []
        
        if frame.type == core.G3FrameType.PipelineInfo:
            if self.collect_by_band_wafer_quantities:
                if ["DroppedBolos", "SurvivingBolos"] <= frame.keys():
                    self.pipe_info_frame = frame
            return []
        
        
        if frame.type == core.G3FrameType.Map:
            
            # - Load a map frame and figure out
            #   whether its data should go into the output file
            #   based on some of its basic information
            #   The frame can contain either coadded maps and
            #   previous analysis results or an individual map
            
            id_for_coadds = None
            
            # * The above will be the same as frame["Id"]
            #   if we don't combine maps from different directions or wafers
            
            if "AreCoaddedMapsContained" in frame.keys():
                if frame["AreCoaddedMapsContained"] == True:
                    frame_has_old_coadds = True
                else:
                    frame_has_old_coadds = False
            else:
                frame_has_old_coadds = False
            
            if frame_has_old_coadds:
                if frame["Id"] in self.id_mapping.values():
                    id_for_coadds = frame["Id"]
                elif frame["Id"] in self.id_mapping.keys():
                    id_for_coadds = self.id_mapping[frame["Id"]]
                
                if id_for_coadds == None:
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because its map ID is ")
                    self.log("* not one of the desired IDs.")
                    self.log("")
                    return []
                else:
                    season = get_season_based_on_fields(
                                 frame["CoaddedObservationIDs"].keys())
                    center_ra = frame["T"].alpha_center
                    center_dec = frame["T"].delta_center
                    """
                    if season == "winter":
                        center_ra  = core.G3Units.deg *   0.0
                        center_dec = core.G3Units.deg * -56.0
                    elif season == "summer":
                        center_ra  = core.G3Units.deg *  75.0
                        center_dec = core.G3Units.deg * -42.0
                    elif season == "summerb":
                        center_ra  = core.G3Units.deg *  25.0
                        center_dec = core.G3Units.deg * -35.0
                    """
                    center_y, center_x = \
                        numpy.unravel_index(
                            frame["T"].angle_to_pixel(center_ra, center_dec),
                            frame["T"].shape)
                    obs_ids_from_this_frame = \
                        {id_for_coadds: frame["CoaddedObservationIDs"]}
                    map_ids_from_this_frame = \
                        {id_for_coadds: frame["CoaddedMapIDs"]}
                    obs_tms_from_this_frame = \
                        {id_for_coadds: frame["ObservationDurations"]}
                    bad_ids_from_this_frame = \
                        {id_for_coadds: frame["IgnoredObservationIDs"]}
                    self.ignored_obs_ids = \
                        self.combine_ids(self.ignored_obs_ids,
                                         bad_ids_from_this_frame)
            
            else:
                if frame["Id"] not in self.id_mapping.keys():
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because its map ID is")
                    self.log("* not one of the desired IDs.")
                    self.log("")
                    return []
                elif self.obs_info is None:
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because it is unclear what type of observation")
                    self.log("* the map is from.")
                    self.log("")
                    return []
                elif self.obs_info["SourceName"] not in self.map_sources:
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because of the wrong sub-field.")
                    self.log("")
                    return []
                else:
                    id_for_coadds = self.id_mapping[frame["Id"]]
                    oid  = self.obs_info["ObservationID"]
                    sbfd = self.obs_info["SourceName"]
                    dmid = str(oid) + "_" + frame["Id"]   # * more detailed id
                    t_i  = self.obs_info["ObservationStart"]
                    t_f  = self.obs_info["ObservationStop"]
                    d_i  = std_processing.time_to_obsid(t_i)
                    d_f  = std_processing.time_to_obsid(t_f)
                    dura = {str(oid): (d_f - d_i) * core.G3Units.s}
                    season = get_season_based_on_fields([sbfd])
                    if season == "winter":
                        center_dec = core.G3Units.deg * float(sbfd[-6:])
                        dec_height = 7.5 * core.G3Units.deg
                    elif season == "summer":
                        center_dec = core.G3Units.deg * float(sbfd[-5:])
                        dec_height = 7.0 * core.G3Units.deg
                    elif season == "summerp":
                        center_dec = core.G3Units.deg * float(sbfd[-6:])
                        dec_height = 3.5 * core.G3Units.deg
                    elif season == "summerb":
                        center_dec = core.G3Units.deg * float(sbfd[-6:])
                        dec_height = 3.5 * core.G3Units.deg
                    elif season == "summerc":
                        center_dec = core.G3Units.deg * float(sbfd[-6:])
                        dec_height = 3.5 * core.G3Units.deg
                    else:
                        raise Exception("Not clear what season to use!")
                    center_ra = frame["T"].alpha_center
                    minimap_length_x = 12.0 * core.G3Units.deg
                    minimap_length_y = 12.0 * core.G3Units.deg
                    center_y, center_x = \
                        numpy.unravel_index(
                            frame["T"].angle_to_pixel(center_ra, center_dec),
                            frame["T"].shape)
                    for b in ["90GHz", "150GHz", "220GHz"]:
                        if b in id_for_coadds:
                            band = b.replace("GHz", "")
                            break
                    obs_ids_from_this_frame = \
                        {id_for_coadds: {sbfd: core.G3VectorInt([oid])}}
                    map_ids_from_this_frame = \
                        {id_for_coadds: {sbfd: core.G3VectorString([dmid])}}
                    obs_tms_from_this_frame = \
                        {id_for_coadds: {sbfd: core.G3MapDouble(dura)}}
            
            
            # - Perform a few more checks on the frame
            #   to furthur figure out whether it can go into the coadds
            
            if not self.map_seems_fine(frame):
                self.log("")
                self.log("* Well, the map doesn't look good,")
                self.log("* so, it will be skipped.")
                self.log("")
                self.ignored_obs_ids = \
                    self.combine_ids(self.ignored_obs_ids,
                                     obs_ids_from_this_frame)
                return []
            
            self.coadded_obs_ids = \
                self.combine_ids(self.coadded_obs_ids,
                                 obs_ids_from_this_frame)
            self.coadded_map_ids, this_frame_already_added = \
                self.combine_ids(self.coadded_map_ids,
                                 map_ids_from_this_frame,
                                 report_commonness=True)
            self.obser_durations = \
                self.combine_mapmapdoubles(self.obser_durations,
                                           obs_tms_from_this_frame)
            
            if this_frame_already_added and \
               (not self.allow_subtraction):
                self.log("")
                self.log("* Well, it looks that this map frame")
                self.log("* has already been added (at least partially),")
                self.log("* so, it will be skipped.")
                self.log("")
                return []
            
            subtract_maps_in_this_frame = this_frame_already_added
            
            
            # - Now that the map is confirmed to be analyzed and added,
            #   various processing steps will act on maps
            
            # -- Add (including plus -1 times) the frame to the existing coadds
            
            if subtract_maps_in_this_frame:
                self.log("")
                self.log("* Removing a set of ")
                self.log("* sky map and weight map (ID: %s)", frame["Id"])
                self.log("* from the coadded maps (ID %s) ...", id_for_coadds)
                self.coadded_obs_ids = \
                    self.remove_ids(self.coadded_obs_ids,
                                    obs_ids_from_this_frame)
                self.coadded_map_ids = \
                    self.remove_ids(self.coadded_map_ids,
                                    map_ids_from_this_frame)
                self.obser_durations = \
                    self.remove_partial_mapmapdouble(self.obser_durations,
                                                     obs_ids_from_this_frame)
                self.coadded_map_frames[id_for_coadds] = \
                    add_two_map_frames(
                        self.coadded_map_frames[id_for_coadds],
                        frame,
                        t_only=self.t_only,
                        remove_weights_beforehand=False,
                        multiply_2nd_maps_by=-1.0,
                        multiply_2nd_weights_by=-1.0,
                        remove_weights_afterward=False,
                        divide_result_by=1.0,
                        record_weights=True,
                        logfun=self.detail,
                        iy=center_y, ix=center_x)
                self.log("* Done.")
                self.log("")
            else:
                self.log("")
                self.log("* Adding a set of ")
                self.log("* sky map and weight map (ID: %s)", frame["Id"])
                self.log("* to the coadded maps (ID: %s) ...", id_for_coadds)
                if len(self.coadded_map_frames[id_for_coadds].keys()) == 0:
                    if self.t_only:
                        self.coadded_map_frames[id_for_coadds]["T"] = frame["T"]
                        self.coadded_map_frames[id_for_coadds]["Wunpol"] = \
                            frame["Wunpol"]
                        tu = core.G3Units.mK
                        self.detail("$$$ t of original map @ "
                                    "(x, y) = (%d, %d): %7.3e"
                                    %(center_x, center_y,
                                      numpy.asarray(
                                      self.coadded_map_frames\
                                      [id_for_coadds]["T"])\
                                      [center_y][center_x]/tu))
                        self.detail("$$$ w of original map @ "
                                    "(x, y) = (%d, %d): %7.3e"
                                    %(center_x, center_y,
                                      numpy.asarray(
                                      self.coadded_map_frames\
                                      [id_for_coadds]["Wunpol"].TT)\
                                      [center_y][center_x]*tu*tu))
                    self.log("* ... field and weight maps were added "
                             "to the empty cache.")
                    maps.CompactMaps(self.coadded_map_frames[id_for_coadds],
                                     zero_nans=True)
                else:
                    self.coadded_map_frames[id_for_coadds] = \
                        add_two_map_frames(
                            self.coadded_map_frames[id_for_coadds],
                            frame,
                            t_only=self.t_only,
                            remove_weights_beforehand=False,
                            multiply_2nd_maps_by=1.0,
                            multiply_2nd_weights_by=1.0,
                            remove_weights_afterward=False,
                            divide_result_by=1.0,
                            record_weights=True,
                            logfun=self.detail,
                            iy=center_y, ix=center_x)
                self.log("* Done.")
                self.log("")
            
            
            # -- Collect flagging and calibration related quantities
            #    for each wafer
            
            # --- Collect averages related to detector flagging
            
            if self.get_avgs_flagging_stats:
                
                if subtract_maps_in_this_frame:
                    for waf in self.wafers:
                        for fr in self.flagging_reasons:
                            self.avgs_flagging_stats[waf][fr] = \
                                self.remove_partial_mapmapdouble(
                                    self.avgs_flagging_stats[waf][fr],
                                    obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering average numbers related to")
                        self.log("* detector flagging that were")
                        self.log("* calculated previously ...")
                        avgs_flg_stats_from_this_fr = {}
                        for key in frame.keys():
                            if self.key_prefix_flg in key:
                                for waf in self.wafers:
                                    avgs_flg_stats_from_this_fr[waf] = {}
                                    for fr in self.flagging_reasons:
                                        fk = self.key_prefix_flg + waf + fr
                                        avgs_flg_stats_from_this_fr[waf][fr] = \
                                            {id_for_coadds: frame[fk]}
                    else:
                        self.log("")
                        self.log("* Gathering average numbers related to")
                        self.log("* detector flagging ...")
                        avgs_flg_stats_from_this_fr = {}
                        avgs_from_each_wafer = \
                            collect_averages_from_flagging_info(
                                self.pipe_info_frame,
                                self.flagging_reasons,
                                self.wafers, band,
                                self.cal_frame["BolometerProperties"])
                        for wafer in avgs_from_each_wafer.keys():
                            avgs_flg_stats_from_this_fr[wafer] = {}
                            for reason in avgs_from_each_wafer[wafer].keys():
                                avg = avgs_from_each_wafer[wafer][reason]
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, avg)
                                avgs_flg_stats_from_this_fr[wafer][reason] = \
                                    {id_for_coadds: mmd}
                    for waf in avgs_flg_stats_from_this_fr.keys():
                        for rs in avgs_flg_stats_from_this_fr[waf].keys():
                            self.avgs_flagging_stats[waf][rs] = \
                                self.combine_mapmapdoubles(
                                    self.avgs_flagging_stats[waf][rs],
                                    avgs_flg_stats_from_this_fr[waf][rs])
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Make a dictionary specifying which detectors
            #     belong to which wafers so that it can be used
            #     when collecting calibration related quantities
            
            if (self.calc_median_pW_to_K or self.calc_cal_vs_el) and \
               (not frame_has_old_coadds)                        and \
               (not subtract_maps_in_this_frame):
                self.log("")
                self.log("* Grouping detectors that were not flagged")
                self.log("* during more than half of the scans by wafer ...")
                self.wafers_and_bolos = \
                     group_suviving_detectors_by_wafer(
                         band, self.wafers,
                         self.cal_frame, self.pipe_info_frame,
                         logfun=self.detail)
                self.log("* Done.")
                self.log("")
            
            
            # --- Calculate mean value of the pW/K conversion factors
            #     for each band of each wafer
            
            if self.calc_median_pW_to_K:
                
                if subtract_maps_in_this_frame:
                    for waf in self.wafers:
                        for fa in self.calfactors:
                            self.medians_of_temp_cal_factors[waf][fa] = \
                                self.remove_partial_mapmapdouble(
                                    self.medians_of_temp_cal_factors[waf][fa],
                                    obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering the median values of the pW/K")
                        self.log("* conversion factors that were")
                        self.log("* calculated previously ...")
                        meds_pwks_from_this_fr = {}
                        for key in frame.keys():
                            if self.key_prefix_pwk in key:
                                for waf in self.wafers:
                                    meds_pwks_from_this_fr[waf] = {}
                                    for fa in self.calfactors:
                                        ck = self.key_prefix_pwk + waf + fa
                                        meds_pwks_from_this_fr[waf][fa] = \
                                            {id_for_coadds: frame[ck]}
                    else:
                        self.log("")
                        self.log("* Gathering the median values of the pW/K")
                        self.log("* conversion factors ...")
                        meds_pwks_from_this_fr = {}
                        meds_to_be_reorganized = \
                            collect_medians_of_pW_per_K_factors(
                                self.cal_frame,
                                self.wafers_and_bolos, band,
                                logfun=self.detail)
                        for wafer in meds_to_be_reorganized.keys():
                            meds_pwks_from_this_fr[wafer] = {}
                            for factor in meds_to_be_reorganized[wafer].keys():
                                med = meds_to_be_reorganized[wafer][factor]
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, med)
                                meds_pwks_from_this_fr[wafer][factor] = \
                                    {id_for_coadds: mmd}
                    for waf in meds_pwks_from_this_fr.keys():
                        for fa in meds_pwks_from_this_fr[waf].keys():
                            self.medians_of_temp_cal_factors[waf][fa] = \
                                self.combine_mapmapdoubles(
                                    self.medians_of_temp_cal_factors[waf][fa],
                                    meds_pwks_from_this_fr[waf][fa])
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate changes in detectors' response to calibrator
            
            if self.calc_cal_vs_el:
                
                if subtract_maps_in_this_frame:
                    for waf in self.wafers:
                        for lc in self.cal_locations:
                            self.cal_resp_diffs_els[waf][lc] = \
                                self.remove_partial_mapmapdouble(
                                    self.cal_resp_diffs_els[waf][lc],
                                    obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering the fractional changes")
                        self.log("* in detectors' response to the calibrator")
                        self.log("* that were calculated previously ...")
                        frac_diffs_from_this_fr = {}
                        for key in frame.keys():
                            if self.key_prefix_calel in key:
                                for waf in self.wafers:
                                    frac_diffs_from_this_fr[waf] = {}
                                    for lc in self.cal_locations:
                                        ck = self.key_prefix_calel + waf + lc
                                        frac_diffs_from_this_fr[waf][lc] = \
                                            {id_for_coadds: frame[ck]}
                    else:
                        self.log("")
                        self.log("* Gathering the fractional changes")
                        self.log("* in detectors' response to the")
                        self.log("* calibrator observations taken at")
                        self.log("* the bottom and top of the field ...")
                        frac_diffs_from_this_fr = {}
                        fracs_to_be_reorganized = \
                            calculate_change_in_cal_response_vs_elevation(
                                oid,
                                center_dec, dec_height,
                                self.cal_obs_ids,
                                self.cal_ts_dir, self.cal_autoproc_dir,
                                self.wafers_and_bolos,
                                logfun=self.detail)
                        self.log("* ... the fractional change of")
                        self.log("* median cal. response at the top "
                                 "with respect to the bottom")
                        self.log("* was %7.3e.",
                                 fracs_to_be_reorganized\
                                 ["AllBolos"]["FractionalChangesTopToBottom"])
                        for waf in fracs_to_be_reorganized.keys():
                            frac_diffs_from_this_fr[waf] = {}
                            for lc in fracs_to_be_reorganized[waf].keys():
                                val = fracs_to_be_reorganized[waf][lc]
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, val)
                                frac_diffs_from_this_fr[waf][lc] = \
                                    {id_for_coadds: mmd}
                    for waf in frac_diffs_from_this_fr.keys():
                        for lc in frac_diffs_from_this_fr[waf].keys():
                            self.cal_resp_diffs_els[waf][lc] = \
                                self.combine_mapmapdoubles(
                                    self.cal_resp_diffs_els[waf][lc],
                                    frac_diffs_from_this_fr[waf][lc])
                    self.log("* Done.")
                    self.log("")
            
            
            # -- Perform analysis tasks actually using the map data
            
            # --- Decide whether it's time to perform some analysis on maps
            
            if self.analyze_maps:
                
                self.log("")
                self.log("* Preparing the maps needed for the "
                         "requested analysis ...")
                
                if frame_has_old_coadds:
                    time_to_analyze_maps = True
                    # * Set this to True just to collect 
                    #   the previous analysis results
                    if self.maps_split_by_scan_direction:
                        n_mp_ea_fld = {fld: len(obs_ids) for fld, obs_ids in \
                                       frame["CoaddedObservationIDs"].items()}
                        for f, n in n_mp_ea_fld.items():
                            self.map_fr_arrival_counters[frame["Id"]][f] += n
                
                else:
                    if self.maps_split_by_scan_direction:
                        time_to_analyze_maps = False
                        self.map_fr_arrival_counters[frame["Id"]][sbfd] += 1
                        self.map_cache_for_both_dirs[frame["Id"]][sbfd] = frame
                        for di_mid in self.direction_independent_map_ids:
                            if "SomeDirection" in di_mid:
                                idl = di_mid.replace("SomeDirection", "Left")
                                idr = di_mid.replace("SomeDirection", "Right")
                            else:
                                idl = "Left"  + di_mid
                                idr = "Right" + di_mid
                            n_mp_l = self.map_fr_arrival_counters[idl][sbfd]
                            n_mp_r = self.map_fr_arrival_counters[idr][sbfd]
                            self.detail("$$$ n_maps_left : %d, n_maps_right : %d"
                                        %(n_mp_l, n_mp_r))
                            if (n_mp_l > 0) and (n_mp_l == n_mp_r):
                                time_to_analyze_maps = True
                                break
                        
                        if time_to_analyze_maps:
                            self.log("* ... since both the %s and %s maps of\n"
                                     "* observation %s have been cached,",
                                     idl, idr, oid)
                            self.log("* the requested analyses can now start.")
                        else:
                            self.log("* ... since either %s or %s map of\n"
                                     "* observation %s has not been cached yet,",
                                     idl, idr, oid)
                            self.log("* the requested analyses need to wait.")
                        
                        if time_to_analyze_maps:
                            summ_map_frame_coadded = \
                                add_two_map_frames(
                                    self.coadded_map_frames[idl],
                                    self.coadded_map_frames[idr],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=False,
                                    multiply_2nd_maps_by=1.0,
                                    multiply_2nd_weights_by=1.0,
                                    remove_weights_afterward=True,
                                    divide_result_by=1.0,
                                    record_weights=True,
                                    logfun=self.detail,
                                    iy=center_y, ix=center_x)
                            """
                            summ_map_frame_coadded = None
                            """
                            """
                            summ_map_frame_individ = \
                                add_two_map_frames(
                                    self.map_cache_for_both_dirs[idl][sbfd],
                                    self.map_cache_for_both_dirs[idr][sbfd],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=False,
                                    multiply_2nd_maps_by=1.0,
                                    multiply_2nd_weights_by=1.0,
                                    remove_weights_afterward=True,
                                    divide_result_by=1.0,
                                    record_weights=True,
                                    logfun=self.detail,
                                    iy=center_y, ix=center_x)
                            """
                            summ_map_frame_individ = None
                            diff_map_frame_coadded = \
                                add_two_map_frames(
                                    self.coadded_map_frames[idl],
                                    self.coadded_map_frames[idr],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=True,
                                    multiply_2nd_maps_by=-1.0,
                                    multiply_2nd_weights_by=1.0,
                                    remove_weights_afterward=False,
                                    divide_result_by=2.0,
                                    record_weights=False,
                                    logfun=self.detail,
                                    iy=center_y, ix=center_x)
                            """
                            diff_map_frame_coadded = None
                            """
                            """
                            diff_map_frame_individ = \
                                add_two_map_frames(
                                    self.map_cache_for_both_dirs[idl][sbfd],
                                    self.map_cache_for_both_dirs[idr][sbfd],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=True,
                                    multiply_2nd_maps_by=-1.0,
                                    multiply_2nd_weights_by=1.0,
                                    remove_weights_afterward=False,
                                    divide_result_by=2.0,
                                    record_weights=False,
                                    logfun=self.detail,
                                    iy=center_y, ix=center_x)
                            """
                            diff_map_frame_individ = None
                            # * Depending on what analysis tasks are done,
                            #   not all of the four types of frames are needed.
                            #   Having no """'s would be ideal, though.
                            
                            summ_map_frame_individ_mini,  \
                            summ_map_frame_coadded_mini,  \
                            diff_map_frame_individ_mini,  \
                            diff_map_frame_coadded_mini = \
                                [create_new_map_frame_with_smaller_region(
                                     one_type_of_frame,
                                     center_ra, center_dec,
                                     minimap_length_x, minimap_length_y,
                                     t_only=self.t_only)
                                 for one_type_of_frame in
                                     (summ_map_frame_individ,
                                      summ_map_frame_coadded,
                                      diff_map_frame_individ,
                                      diff_map_frame_coadded)]
                    else:
                        time_to_analyze_maps   = True
                        summ_map_frame_individ = frame
                        tu = core.G3Units.mK
                        if self.t_only:
                            maps.RemoveWeights(frame, zero_nans=True)
                            self.detail("$$$ t of original map "
                                        "@ field center : %7.3e [mK]"
                                        %(numpy.asarray(
                                          summ_map_frame_individ["T"])\
                                          [center_y][center_x]/tu))
                            maps.CompactMaps(summ_map_frame_individ,
                                             zero_nans=True)
                        summ_map_frame_coadded = None
                        diff_map_frame_individ = frame
                        diff_map_frame_coadded = None
                        # * Depending on what analysis tasks are done,
                        #   not all of the four types of frames are needed.
                        #   Having no """'s would be ideal, though.
                        
                        summ_map_frame_individ_mini,  \
                        summ_map_frame_coadded_mini,  \
                        diff_map_frame_individ_mini,  \
                        diff_map_frame_coadded_mini = \
                            [create_new_map_frame_with_smaller_region(
                                 one_type_of_frame,
                                 center_ra, center_dec,
                                 minimap_length_x, minimap_length_y,
                                 t_only=self.t_only)
                             for one_type_of_frame in
                                 (summ_map_frame_individ,
                                  summ_map_frame_coadded,
                                  diff_map_frame_individ,
                                  diff_map_frame_coadded)]
                self.log("* Preparations done.")
                self.log("")
            
            else:
                time_to_analyze_maps = False



            # --- Calculate a standard deviation of the temperature map,
            #     a mean of the TT weight map, and the number of pixels in the
            #     weight map whose weights are above certain threshold
            #     for each applicable map frame
            
            if time_to_analyze_maps and self.calc_map_stddev_etc:
                
                if subtract_maps_in_this_frame:
                    for mt in self.map_types_for_flc:
                        for fk in self.fluctuation_keys:
                            self.map_fluctuation_metrics[mt][fk] = \
                                self.remove_partial_mapmapdouble(
                                    self.map_fluctuation_metrics[mt][fk],
                                    obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering map fluctuation metrics")
                        self.log("* that were calculated previously ...")
                        for mt in self.map_types_for_flc:
                            for fk in self.fluctuation_keys:
                                mtfk = self.key_prefix_flc + mt + fk
                                self.map_fluctuation_metrics[mt][fk] = \
                                    self.combine_mapmapdoubles(
                                        self.map_fluctuation_metrics[mt][fk],
                                        {id_for_coadds: frame[mtfk]})
                    else:
                        self.log("")
                        self.log("* Calculating map fluctuation metrics")
                        self.log("* for the following types of maps:")
                        self.log("* %s ...", self.map_types_for_flc)
                        if self.pixels_to_use_for_flc_calc[sbfd] is None:
                            self.log("* (Need to figure out what pixels "
                                     "to use for the calculations,")
                            if summ_map_frame_coadded is None:
                                fr_for_getting_pixels = summ_map_frame_individ
                            else:
                                fr_for_getting_pixels = summ_map_frame_coadded
                            self.pixels_to_use_for_flc_calc[sbfd] = \
                                identify_pixels_of_non_atypical_region(
                                    sbfd, fr_for_getting_pixels,
                                    self.point_source_list_file)
                            self.log("*  the results of which "
                                     "will be repeatedly used later.)")
                        for mt in self.map_types_for_flc:
                            if mt == "IndividualSignalMaps":
                                frame_to_use = summ_map_frame_individ
                            if mt == "CoaddedSignalMaps":
                                frame_to_use = summ_map_frame_coadded
                            if mt == "IndividualNoiseMaps":
                                frame_to_use = diff_map_frame_indvid
                            if mt == "CoaddedNoiseMaps":
                                frame_to_use = diff_map_frame_coadded
                            self.log("* Results for %s :", mt)
                            fluctuation_metrics = \
                                calculate_map_fluctuation_metrics(
                                    frame_to_use, band, sbfd,
                                    pixs=self.pixels_to_use_for_flc_calc[sbfd],
                                    t_only=self.t_only,
                                    logfun=self.detail)
                            tu = core.G3Units.mK
                            for k, v in fluctuation_metrics.items():
                                self.log("*  %35s : %s %s",
                                         k, v/self.flu_met_unis[k],
                                         self.flu_met_ustr[k])
                            for k, v in fluctuation_metrics.items():
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, v)
                                self.map_fluctuation_metrics[mt][k] = \
                                    self.combine_mapmapdoubles(
                                        self.map_fluctuation_metrics[mt][k],
                                        {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate pointing discrepancies for three bright sources
            #     in the sub-field
            
            if time_to_analyze_maps and self.calc_pointing_discrepancies:
                
                if subtract_maps_in_this_frame:
                    for r in self.ptsrc_ranks:
                        for d in self.ptsrc_offset_types:
                            self.ptsrc_delta_ras_and_decs[r][d] = \
                                self.remove_partial_mapmapdouble(
                                    self.ptsrc_delta_ras_and_decs[r][d],
                                    obs_ids_from_this_frame)
                        self.ptsrc_fluxes[r] = \
                            self.remove_partial_mapmapdouble(
                                self.ptsrc_fluxes[r],
                                obs_ids_from_this_frame)
                        self.ptsrc_snrs[r] = \
                            self.remove_partial_mapmapdouble(
                                self.ptsrc_snrs[r],
                                obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering data related to point sources")
                        self.log("* that were calculated previously ...")
                        for r in self.ptsrc_ranks:
                            for d in self.ptsrc_offset_types:
                                dk = self.key_template_ptsrc_d.\
                                     replace("Delta", d).replace("Nth", r)
                                self.ptsrc_delta_ras_and_decs[r][d] = \
                                    self.combine_mapmapdoubles(
                                        self.ptsrc_delta_ras_and_decs[r][d],
                                        {id_for_coadds: frame[dk]})
                            fk = self.key_template_ptsrc_f.replace("Nth", r)
                            self.ptsrc_fluxes[r] = \
                                self.combine_mapmapdoubles(
                                    self.ptsrc_fluxes[r],
                                    {id_for_coadds: frame[fk]})
                            sk = self.key_template_ptsrc_s.replace("Nth", r)
                            self.ptsrc_snrs[r] = \
                                self.combine_mapmapdoubles(
                                    self.ptsrc_snrs[r],
                                    {id_for_coadds: frame[sk]})
                    else:
                        self.log("")
                        self.log("* Calculating pointing offsets of")
                        self.log("* 3 bright point sources in the field ...")
                        if self.calc_map_stddev_etc and \
                           ("IndividualSignalMaps" in self.map_types_for_flc):
                            map_stddev = self.map_fluctuation_metrics   \
                                             ["IndividualSignalMaps"]   \
                                             ["TMapStandardDeviations"] \
                                             [id_for_coadds][sbfd][str(oid)]
                        else:
                            map_stddev = None
                        thumbnails_file = os.path.join(
                            self.field_autoproc_dir.replace("FIELD", sbfd), 
                            "bright_sources", str(oid)+".g3")
                        thumbnails_band = frame["Id"]
                        offsets, fluxes, snrs = \
                            calculate_pointing_discrepancies(
                                summ_map_frame_individ, sbfd,
                                map_stddev=map_stddev,
                                thumbnails_file=thumbnails_file,
                                thumbnails_band=thumbnails_band)
                        for rank, diffs in offsets.items():
                            u = core.G3Units
                            for coord, diff in diffs.items():
                                self.log("* %s Source - %9s : %+7.3e arcsec.",
                                         rank, coord, diff/u.arcsec)
                            flux  = fluxes[rank]
                            flux /= u.mK * u.arcmin * u.arcmin
                            self.log("* %s Source - %9s : %+7.3e mK.arcmin.^2",
                                     rank, "flux", flux)
                            snr = snrs[rank]
                            self.log("* %s Source - %9s : %+7.3e",
                                     rank, "SNR", snr)
                        for r in self.ptsrc_ranks:
                            for d in self.ptsrc_offset_types:
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, offsets[r][d])
                                self.ptsrc_delta_ras_and_decs[r][d] = \
                                    self.combine_mapmapdoubles(
                                        self.ptsrc_delta_ras_and_decs[r][d],
                                        {id_for_coadds: mmd})
                            mmd = self.create_mmd_for_one_value(
                                      sbfd, oid, fluxes[r])
                            self.ptsrc_fluxes[r] = \
                                self.combine_mapmapdoubles(
                                    self.ptsrc_fluxes[r],
                                    {id_for_coadds: mmd})
                            mmd = self.create_mmd_for_one_value(
                                      sbfd, oid, snrs[r])
                            self.ptsrc_snrs[r] = \
                                self.combine_mapmapdoubles(
                                    self.ptsrc_snrs[r],
                                    {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate the ratio of SPT x Planck spectrum to
            #     Planck x Planck spectrum
            
            if time_to_analyze_maps and self.calc_xspec_with_plck_map:
                
                if subtract_maps_in_this_frame:
                    for x_typ in self.xspec_types:
                        self.avgs_xspec_ratios[x_typ] = \
                            self.remove_partial_mapmapdouble(
                                self.avgs_xspec_ratios[x_typ],
                                obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering averages of ratios of")
                        self.log("* cross spectra (SPT x Planck / Pla x Pla)")
                        self.log("* that were calculated previously ...")
                        for x_typ in self.xspec_types:
                            xk = self.av_xspec_ratio_key + x_typ + "spectra"
                            self.avgs_xspec_ratios[x_typ] = \
                                self.combine_mapmapdoubles(
                                    self.avgs_xspec_ratios[x_typ],
                                    {id_for_coadds: frame[xk]})
                    else:
                        self.log("")
                        self.log("* Calculating averages of ratios of")
                        self.log("* SPT x Planck / Planck x Planck spectra")
                        self.log("* in a low ell range ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_of_small_region(
                                    summ_map_frame_individ_mini,
                                    self.point_source_list_file,
                                    center_ra, center_dec, dec_height)
                        for mission in self.planck_missions:
                            if self.mini_planck_map_frames\
                               [mission][band][sbfd] is None:
                                self.log("* (Need to make a %s "
                                         "mini Planck map,", mission)
                                planck_band = \
                                    self.spt_to_planck_bands[band]
                                fits_file = \
                                    self.planck_map_fits_file.\
                                    replace("BAND", planck_band+"GHz").\
                                    replace("MISSION", mission)
                                self.mini_planck_map_frames\
                                [mission][band][sbfd] = \
                                    create_mini_planck_map_frame(
                                        fits_file, summ_map_frame_individ,
                                        center_ra, center_dec,
                                        minimap_length_x, minimap_length_y,
                                        t_only=self.t_only)
                                self.log("*  which will be used repeatedly.")
                        avg_xspec_ratios = \
                            calculate_average_ratio_of_spt_planck_xspectra(
                                summ_map_frame_individ_mini,
                                {mission: self.mini_planck_map_frames \
                                          [mission][band][sbfd]       \
                                 for mission in self.planck_missions},
                                self.masks_for_powspec_calculations[sbfd],
                                t_only=self.t_only)
                        self.log("* ... the average was calculated to be %s.",
                                 avg_xspec_ratios)
                        for x_typ in self.xspec_types:
                            mmd = self.create_mmd_for_one_value(
                                      sbfd, oid, avg_xspec_ratios[x_typ])
                            self.avgs_xspec_ratios[x_typ] = \
                                self.combine_mapmapdoubles(
                                    self.avgs_xspec_ratios[x_typ],
                                    {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate noise levels from individual maps
            
            if time_to_analyze_maps and self.calc_noise_from_individ_maps:
                
                if subtract_maps_in_this_frame:
                    for mt in self.map_types_for_noise_calc:
                        self.noise_from_individual_maps[mt] = \
                            self.remove_partial_mapmapdouble(
                                self.noise_from_individual_maps[mt],
                                obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering noise levels of individual maps ")
                        self.log("* that were calculated previously ...")
                        for mt in self.map_types_for_noise_calc:
                            nk = self.indvid_noise_key.replace("{}", mt)
                            self.noise_from_individual_maps[mt] = \
                                self.combine_mapmapdoubles(
                                    self.noise_from_individual_maps[mt],
                                    {id_for_coadds: frame[nk]})
                    else:
                        self.log("")
                        self.log("* Calculating noise level "
                                 "from this observation ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_of_small_region(
                                    diff_map_frame_individ_mini,
                                    self.point_source_list_file,
                                    center_ra, center_dec, dec_height)
                        noises = calculate_noise_levels(
                                     diff_map_frame_individ_mini,
                                     self.masks_for_powspec_calculations[sbfd],
                                     t_only=self.t_only)
                        nu = core.G3Units.uK * core.G3Units.arcmin
                        self.log("* ... the noise level was calculated to be ")
                        self.log("* %s uK.arcmin.",
                                 {k: v/nu for k, v in noises.items()})
                        for mt in self.map_types_for_noise_calc:
                            mmd = self.create_mmd_for_one_value(
                                      sbfd, oid, noises[mt])
                            self.noise_from_individual_maps[mt] = \
                                self.combine_mapmapdoubles(
                                    self.noise_from_individual_maps[mt],
                                    {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate noise levels from coadded maps
            
            if time_to_analyze_maps and self.calc_noise_from_coadded_maps:
                
                if subtract_maps_in_this_frame:
                    for mt in self.map_types_for_noise_calc:
                        self.noise_from_coadded_maps[mt] = \
                            self.remove_partial_mapmapdouble(
                                self.noise_from_coadded_maps[mt],
                                obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering noise levels of coadded maps ")
                        self.log("* that were calculated previously ...")
                        for mt in self.map_types_for_noise_calc:
                            nk = self.coadds_noise_key.replace("{}", mt)
                            self.noise_from_coadded_maps[mt] = \
                                self.combine_mapmapdoubles(
                                    self.noise_from_coadded_maps[mt],
                                    {id_for_coadds: frame[nk]})
                    else:
                        self.log("")
                        self.log("* Calculating noise level "
                                 "from the running coadded maps ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_of_small_region(
                                    diff_map_frame_coadded_mini,
                                    self.point_source_list_file,
                                    center_ra, center_dec, dec_height)
                        noises = calculate_noise_levels(
                                     diff_map_frame_coadded_mini,
                                     self.masks_for_powspec_calculations[sbfd],
                                     t_only=self.t_only)
                        nu = core.G3Units.uK * core.G3Units.arcmin
                        self.log("* ... the noise level was calculated to be ")
                        self.log("* %s uK.arcmin.",
                                 {k: v/nu for k, v in noises.items()})
                        for mt in self.map_types_for_noise_calc:
                            mmd = self.create_mmd_for_one_value(
                                      sbfd, oid, noises[mt])
                            self.noise_from_coadded_maps[mt] = \
                                self.combine_mapmapdoubles(
                                    self.noise_from_coadded_maps[mt],
                                    {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Reinitialize some of the variables used when
            #     making preparations for the map analyses
            
            if time_to_analyze_maps:
                
                summ_map_frame_individ = None
                diff_map_frame_individ = None
                summ_map_frame_coadded = None
                diff_map_frame_coadded = None
                summ_map_frame_individ_mini = None
                diff_map_frame_individ_mini = None
                summ_map_frame_coadded_mini = None
                diff_map_frame_coadded_mini = None
                
                if (not frame_has_old_coadds) and \
                   self.maps_split_by_scan_direction:
                    self.map_fr_arrival_counters[idl][sbfd] = 0
                    self.map_fr_arrival_counters[idr][sbfd] = 0
                    self.map_cache_for_both_dirs[idl][sbfd] = None
                    self.map_cache_for_both_dirs[idr][sbfd] = None
            
            
            # - Drop the frame now that it was added and processed
            
            return []
        
        
        
        if frame.type == core.G3FrameType.EndProcessing:
            
            self.log("\n")
            self.log("* It's time to combine all the data gathered so far!")
            self.log("\n")
            
            meta_info_frames_to_return = []
            
            for map_id in self.map_ids:
                mp_fr = self.coadded_map_frames[map_id]
                if len(mp_fr.keys()) == 0:
                    continue
                
                self.log("")
                self.log("# --------------------------- #")
                self.log("# Map ID: %s", map_id )
                self.log("# --------------------------- #")
                self.log("")
                
                mp_fr["Id"] = map_id
                mp_fr["AreCoaddedMapsContained"] = True
                mp_fr["CoaddedMapIDs"] = self.coadded_map_ids[map_id]
                mp_fr["CoaddedObservationIDs"] = self.coadded_obs_ids[map_id]
                mp_fr["IgnoredObservationIDs"] = self.ignored_obs_ids[map_id]
                mp_fr["ObservationDurations"]  = self.obser_durations[map_id]
                
                self.log("# Observation IDs used:")
                self.log("")
                for sf, oids \
                in  mp_fr["CoaddedObservationIDs"].items():
                    mp_fr["NumberOfCoaddedMapsOf"+sf.upper()] = len(oids)
                    self.log("- %s", sf)
                    self.log("")
                    self.log(" "*3 +" %s", list(oids))
                    self.log("")
                self.log("\n")
                
                self.log("# Observation IDs ignored:")
                self.log("")
                for sf, oids \
                in  mp_fr["IgnoredObservationIDs"].items():
                    self.log("- %s", sf)
                    self.log("")
                    self.log(" "*3 +" %s", list(oids))
                    self.log("")
                self.log("\n")
                
                self.log("# Map IDs used:")
                self.log("")
                for sf, mids \
                in  mp_fr["CoaddedMapIDs"].items():
                    self.log("- %s", sf)
                    self.log("")
                    self.log(" "*3 + " %s", list(mids))
                    self.log("")
                self.log("\n")
                
                self.log("# Observation durations [minutes]:")
                du = core.G3Units.min   # * display units
                self.print_mapmapdouble(mp_fr["ObservationDurations"], du, 0)
                self.log("\n")
                
                
                if self.get_avgs_flagging_stats:
                    self.log("# Avg. numbers related to flagging statistics:")
                    self.log("")
                    
                    for waf in self.wafers:
                        for fr in self.flagging_reasons:
                            k_rec = self.key_prefix_flg + waf + fr
                            mp_fr[k_rec] = self.avgs_flagging_stats \
                                           [waf][fr][map_id]
                            if waf == "AllBolos":
                                self.log("- %s", fr)
                                self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_median_pW_to_K:
                    self.log("# Medians of pW/K conversion factors [pW/K]:")
                    self.log("")
                    
                    for waf in self.wafers:
                        self.log(" - %s: ", waf)
                        for fa in self.calfactors:
                            k_rec = self.key_prefix_pwk + waf + fa
                            mp_fr[k_rec] = self.medians_of_temp_cal_factors \
                                           [waf][fa][map_id]
                            if fa == "PicowattsPerKelvin":
                                du = self.cal_fat_unis[fa]
                                self.print_mapmapdouble(mp_fr[k_rec], du, 6)
                    self.log("\n")
                
                
                if self.calc_cal_vs_el:
                    self.log("# Fractional change in calibrator response:")
                    self.log("")
                    
                    for waf in self.wafers:
                        self.log(" - %s: ", waf)                        
                        for lc in self.cal_locations:
                            k_rec = self.key_prefix_calel + waf + lc
                            mp_fr[k_rec] = self.cal_resp_diffs_els \
                                           [waf][lc][map_id]
                            if lc  == "FractionalChangesTopToBottom":
                                self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_map_stddev_etc:
                    self.log("# Basic metrics for fluctuations of map values:")
                    self.log("")
                    
                    k_rec_root = self.key_prefix_flc
                    for mt in self.map_types_for_flc:
                        self.log("- Map type: %s", mt)
                        self.log("")
                        
                        for fk in self.fluctuation_keys:
                            k_rec = k_rec_root + mt + fk
                            mp_fr[k_rec] = \
                                self.map_fluctuation_metrics[mt][fk][map_id]
                            self.log(" "*3 + " - Metric: %s [%s]",
                                     fk, self.flu_met_ustr[fk])
                            self.print_mapmapdouble(
                                mp_fr[k_rec], self.flu_met_unis[fk], 6)
                    self.log("\n")
                
                
                if self.calc_pointing_discrepancies:
                    self.log("# Pointing offsets, fluxes, and SNRs of sources:")
                    self.log("")
                    
                    for r in self.ptsrc_ranks:
                        self.log("- %s brightest point sources", r)
                        self.log("")
                        
                        du = core.G3Units.arcsec
                        for o in self.ptsrc_offset_types:
                            k_rec = self.key_template_ptsrc_d.\
                                    replace("Delta", o).replace("Nth", r)
                            mp_fr[k_rec] = \
                                self.ptsrc_delta_ras_and_decs[r][o][map_id]
                            self.log(" "*3 + " - %s [arcseconds]", o)
                            self.print_mapmapdouble(mp_fr[k_rec], du, 6)
                        
                        k_rec = self.key_template_ptsrc_f.replace("Nth", r)
                        mp_fr[k_rec] = self.ptsrc_fluxes[r][map_id]
                        self.log(" "*3 + " - Flux [mK.arcmin^2]")
                        usys = core.G3Units
                        du   = usys.mK * usys.arcmin * usys.arcmin
                        self.print_mapmapdouble(mp_fr[k_rec], du, 6)
                        
                        k_rec = self.key_template_ptsrc_s.replace("Nth", r)
                        mp_fr[k_rec] = self.ptsrc_snrs[r][map_id]                      
                        self.log(" "*3 + " - SNR")
                        self.print_mapmapdouble(mp_fr[k_rec], 1.0, 6)
                    self.log("\n")
                
                
                if self.calc_xspec_with_plck_map:
                    self.log("# Average ratio of SxP / PxP cross spectra:")
                    self.log("")
                    
                    for x_typ in self.xspec_types:
                        k_rec = self.av_xspec_ratio_key + x_typ + "spectra"
                        mp_fr[k_rec] = self.avgs_xspec_ratios[x_typ][map_id]
                        self.log("- Spectra type: %s", x_typ)
                        self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_noise_from_individ_maps:
                    self.log("# Noise levels from individual maps [uK.arcmin]:")
                    self.log("")
                    
                    for m_typ in self.map_types_for_noise_calc:
                        k_rec = self.indvid_noise_key.replace("{}", m_typ)
                        mp_fr[k_rec] = \
                            self.noise_from_individual_maps[m_typ][map_id]
                        self.log("- Map type: %s", m_typ)
                        du = core.G3Units.uK * core.G3Units.arcmin
                        self.print_mapmapdouble(mp_fr[k_rec], du, 3)
                    self.log("\n")
                
                
                if self.calc_noise_from_coadded_maps:
                    self.log("# Noise levels from coadded maps [uK.arcmin]:")
                    self.log("")
                    
                    for m_typ in self.map_types_for_noise_calc:
                        k_rec = self.coadds_noise_key.replace("{}", m_typ)
                        mp_fr[k_rec] = \
                            self.noise_from_coadded_maps[m_typ][map_id]
                        self.log("- Map type: %s", m_typ)
                        du = core.G3Units.uK * core.G3Units.arcmin
                        self.print_mapmapdouble(mp_fr[k_rec], du, 3)
                    self.log("\n")
                
                
                meta_info_frame = core.G3Frame()
                for k in mp_fr.keys():
                    if k not in ["T", "Q", "U", "Wunpol", "Wpol",
                                 "AreCoaddedMapsContained"]:
                        meta_info_frame[k] = mp_fr[k]
                meta_info_frame["AreCoaddedMapsContained"] = False
                meta_info_frames_to_return.append(meta_info_frame)
                
                
                self.log("# Here is what the frame looks like:")
                self.log("")
                self.log(mp_fr)
                self.log("")
                
                
                self.log("# --------------------------- #")
            self.log("\n")
            
            return meta_info_frames_to_return             + \
                   list(self.coadded_map_frames.values()) + \
                   [frame]


# ==============================================================================




# ==============================================================================
# Define some miscellaneous modules
# ------------------------------------------------------------------------------


def drop_duplicate_analysis_results(frame):
    
    if "AreCoaddedMapsContained" in frame.keys():
        if frame["AreCoaddedMapsContained"] == False:
            return []




def do_weights_look_bad(mean_wt, max_wt, band, sub_field):
    
    wu = 1 / (core.G3Units.mK * core.G3Units.mK)
    thresholds = {"90" : {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":  90.0*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi": 105.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi": 128.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi": 174.0*wu},
                          "ra5hdec-24.5": {"lo": 0.0*wu, "hi": 170.0*wu},
                          "ra5hdec-31.5": {"lo": 0.0*wu, "hi": 188.0*wu},
                          "ra5hdec-38.5": {"lo": 0.0*wu, "hi": 166.0*wu},
                          "ra5hdec-45.5": {"lo": 0.0*wu, "hi": 170.0*wu},
                          "ra5hdec-52.5": {"lo": 0.0*wu, "hi": 220.0*wu},
                          "ra5hdec-59.5": {"lo": 0.0*wu, "hi": 240.0*wu},
                          "ra5hdec-29.75": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra5hdec-33.25": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra5hdec-36.75": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra5hdec-40.25": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra1h40dec-29.75": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra1h40dec-33.25": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra1h40dec-36.75": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra1h40dec-40.25": {"lo": 0.0*wu, "hi": 300.0*wu},
                          "ra12h30dec-29.75": {"lo": 0.0*wu, "hi": 200.0*wu},
                          "ra12h30dec-33.25": {"lo": 0.0*wu, "hi": 200.0*wu},
                          "ra12h30dec-36.75": {"lo": 0.0*wu, "hi": 200.0*wu},
                          "ra12h30dec-40.25": {"lo": 0.0*wu, "hi": 200.0*wu}},
                  "150": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi": 123.0*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi": 160.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi": 185.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi": 260.0*wu},
                          "ra5hdec-24.5": {"lo": 0.0*wu, "hi": 130.0*wu},
                          "ra5hdec-31.5": {"lo": 0.0*wu, "hi": 160.0*wu},
                          "ra5hdec-38.5": {"lo": 0.0*wu, "hi": 160.0*wu},
                          "ra5hdec-45.5": {"lo": 0.0*wu, "hi": 185.0*wu},
                          "ra5hdec-52.5": {"lo": 0.0*wu, "hi": 232.0*wu},
                          "ra5hdec-59.5": {"lo": 0.0*wu, "hi": 310.0*wu},
                          "ra5hdec-29.75": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra5hdec-33.25": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra5hdec-36.75": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra5hdec-40.25": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra1h40dec-29.75": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra1h40dec-33.25": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra1h40dec-36.75": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra1h40dec-40.25": {"lo": 0.0*wu, "hi": 350.0*wu},
                          "ra12h30dec-29.75": {"lo": 0.0*wu, "hi": 235.0*wu},
                          "ra12h30dec-33.25": {"lo": 0.0*wu, "hi": 235.0*wu},
                          "ra12h30dec-36.75": {"lo": 0.0*wu, "hi": 235.0*wu},
                          "ra12h30dec-40.25": {"lo": 0.0*wu, "hi": 235.0*wu}},
                  "220": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":  9.5*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi": 12.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi": 15.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi": 21.0*wu},
                          "ra5hdec-24.5": {"lo": 0.0*wu, "hi": 12.0*wu},
                          "ra5hdec-31.5": {"lo": 0.0*wu, "hi": 13.0*wu},
                          "ra5hdec-38.5": {"lo": 0.0*wu, "hi": 14.0*wu},
                          "ra5hdec-45.5": {"lo": 0.0*wu, "hi": 15.0*wu},
                          "ra5hdec-52.5": {"lo": 0.0*wu, "hi": 18.0*wu},
                          "ra5hdec-59.5": {"lo": 0.0*wu, "hi": 26.0*wu},
                          "ra5hdec-29.75": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra5hdec-33.25": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra5hdec-36.75": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra5hdec-40.25": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra1h40dec-29.75": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra1h40dec-33.25": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra1h40dec-36.75": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra1h40dec-40.25": {"lo": 0.0*wu, "hi": 30.0*wu},
                          "ra12h30dec-29.75": {"lo": 0.0*wu, "hi": 20.0*wu},
                          "ra12h30dec-33.25": {"lo": 0.0*wu, "hi": 20.0*wu},
                          "ra12h30dec-36.75": {"lo": 0.0*wu, "hi": 20.0*wu},
                          "ra12h30dec-40.25": {"lo": 0.0*wu, "hi": 20.0*wu}}}
    # * Thresholds were decided in the following ways:
    #       Winter subfields: around 1.5 x median weights from 2020 data
    #       Summer subfields: around 2.0 x median weights from 2019 data
    #       Summer-b subfields: around 2.0 x median weights from 2020 data (1st week)
    #       Summer-p subfields: same as the numbers from Summer-b subfields
    #         (Summer-p refers to the subdivided Summer el1 & el2 subfield)
    #       Summer-c subfields: 2/3 of the numbers from Summer-b subfields
    #         (Each Summer-c subfield is larger than
    #          the corresponding Summer-b subfield by 50%,
    #          but an observation of the former and an observation of the latter
    #          last for the same amount of time.)
    
    if (mean_wt < thresholds[band][sub_field]["lo"]) or \
       (mean_wt > thresholds[band][sub_field]["hi"]) or \
       (max_wt  > thresholds[band][sub_field]["hi"] * 5.0):
        return True
    else:
        return False




def record_bad_obs_id(text_file, map_id, sub_field, obs_id, bad_reason):
    
    def pad_by_space(string):
        return " "*2 + string + " "*2
    
    words_to_write = []
    words_to_write.append(pad_by_space(map_id.rjust(6)))
    words_to_write.append(pad_by_space(sub_field.rjust(8)))
    words_to_write.append(pad_by_space(str(obs_id).rjust(9)))
    words_to_write.append(" "*2 + bad_reason)
    line_to_write = "|".join(words_to_write) + "\n"
    
    with open(text_file, "a") as file_object:
        file_object.write(line_to_write)




class FlagBadMaps(object):
    
    def __init__(self, map_ids="", t_only=True,
                 bad_map_list_file=None, point_source_list_file=None,
                 auxiliary_files_directory=None,
                 logging_function=print):
        
        self.t_only  = t_only        
        self.sb_fld  = None
        self.obs_id  = None
        self.map_ids = map_ids
        self.p_list  = point_source_list_file
        self.x_list  = bad_map_list_file
        self.aux_dir = auxiliary_files_directory
        self.logfun  = logging_function
        self.pixels_to_use_for_flc_calc = {}
        assert self.x_list is not None, "File to record bad maps not supplied!"    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Observation:
            self.sb_fld = frame["SourceName"]
            self.obs_id = frame["ObservationID"]
            if self.sb_fld not in self.pixels_to_use_for_flc_calc.keys():
                pkl_path = "pixel_numbers_for_calculating_"+\
                           "fluctuation_metrics_of_{}_sub_field.pickle".\
                           format(self.sb_fld)
                pkl_path = os.path.join(self.aux_dir, pkl_path)
                if os.path.isfile(pkl_path):
                    with open(pkl_path, "rb") as f_obj:
                        indices = pickle.load(f_obj)
                    self.pixels_to_use_for_flc_calc[self.sb_fld] = indices
                else:
                    self.pixels_to_use_for_flc_calc[self.sb_fld] = None
        
        if frame.type == core.G3FrameType.Map:
            if frame["Id"] not in self.map_ids:
                return
            
            self.logfun("")
            self.logfun("* A relevant map frame has arrived.")
            self.logfun("* Checking whether the map is reasonable ...")
            
            for b in ["90GHz", "150GHz", "220GHz"]:
                if b in frame["Id"]:
                    band = b.replace("GHz", "")
                    break
            
            # - Remove weights and assume the map is not bad
            
            maps.RemoveWeights(frame, zero_nans=True)
            map_is_bad = False
            
            
            # -- One criterion for cutting a map is
            #    mean TT weight of a map
            
            if self.pixels_to_use_for_flc_calc[self.sb_fld] is None:
                self.logfun("### Figuring out what pixels to use")
                self.logfun("### to calculate mean weight etc...")
                self.pixels_to_use_for_flc_calc[self.sb_fld] = \
                    identify_pixels_of_non_atypical_region(
                        self.sb_fld, frame, self.p_list)
            
            flc_dct = calculate_map_fluctuation_metrics(
                          frame, band, self.sb_fld,
                          pixs=self.pixels_to_use_for_flc_calc[self.sb_fld],
                          t_only=self.t_only,
                          logfun=self.logfun)
            
            mean_wt = flc_dct["MeansOfTTWeights"]
            max_wt  = numpy.max(frame["Wunpol"].TT)
            bad_weights = do_weights_look_bad(
                              mean_wt, max_wt, band, self.sb_fld)
            if bad_weights:
                wu = 1.0 / (core.G3Units.mK * core.G3Units.mK)
                self.logfun("### The weights were determined to be bad!")
                self.logfun("### For the record:")
                self.logfun("###   mean weight: %f mK^-2", mean_wt/wu)
                self.logfun("###   max weight : %f mK^-2", max_wt/wu)
                record_bad_obs_id(self.x_list,
                                  frame["Id"], self.sb_fld, self.obs_id,
                                  "Anomalously large weights.")
            
            
            # -- Another criterion is presence of NaNs
            
            has_nans = False
            maps.ApplyWeights(frame)
            
            if self.t_only:
                map_types = ["T"]
            
            for map_type in map_types:
                if False in numpy.isfinite(numpy.asarray(frame[map_type])):
                    has_nans = True
                    self.logfun("### NaNs were found in the %s map!", map_type)
                    record_bad_obs_id(self.x_list,
                                      frame["Id"], self.sb_fld, self.obs_id,
                                      "NaNs in "+map_type+" map.")
            
            
            # - Report the badness and make the maps sparse again
            
            map_is_bad = bad_weights or has_nans
            
            if map_is_bad:
                frame["IgnoreThisMap"] = True
                self.logfun("* ... the map looks bad!")
            else:
                self.logfun("* ... the map doesn't look obviously bad!")
            self.logfun("")
            
            maps.CompactMaps(frame, zero_nans=True)




class AppendDirectionsToMapIDs(object):
    
    def __init__(self, sub_fields, ids_to_modify, logging_function):
        
        self.ids_to_modify = ids_to_modify
        self.mod_history = \
            {sub_field: {map_id: {"left": 0, "right": 0} \
                         for map_id in ids_to_modify}    \
             for sub_field in sub_fields}
        self.sb_fld = None
        self.logfun = logging_function
    
    def __call__(self, frame):
        
        if "CoaddedObservationIDs" in frame.keys():
            for sb_fld, obs_ids in frame["CoaddedObservationIDs"].items():
                n_prev_mod = len(obs_ids)
                if "Left" in frame["Id"]:
                    id_no_dir = frame["Id"].replace("Left", "")
                    self.mod_history[sb_fld][id_no_dir]["left"]  = n_prev_mod
                elif "Right" in frame["Id"]:
                    id_no_dir = frame["Id"].replace("Right", "")
                    self.mod_history[sb_fld][id_no_dir]["right"] = n_prev_mod
        
        if frame.type == core.G3FrameType.Observation:
            self.sb_fld = frame["SourceName"]
        
        if (frame.type == core.G3FrameType.Map) and \
           (frame["Id"] in self.ids_to_modify):
            self.logfun("")
            self.logfun("* Assigning a fake ID to the map ...")
            
            old_id = frame.pop("Id")
            if "IgnoreThisMap" in frame.keys():
                self.logfun("* (This map seems to be a bad one!")
                self.logfun("*  A fake ID will be assigned anyway,")
                self.logfun("*  but it will not be used for coadding.)")
                frame["Id"] = "Left" + old_id
                self.logfun("* Done.")
                self.logfun("")
            else:
                n_mp_l = self.mod_history[self.sb_fld][old_id]["left"]
                n_mp_r = self.mod_history[self.sb_fld][old_id]["right"]
                if n_mp_l > n_mp_r:
                    frame["Id"] = "Right" + old_id
                    self.mod_history[self.sb_fld][old_id]["right"] += 1
                else:
                    frame["Id"] = "Left" + old_id
                    self.mod_history[self.sb_fld][old_id]["left"]  += 1
                self.logfun("* Done.")
                self.logfun("")




def print_analysis_results(
        filepath, info_tuples, print_frame=True, output_path=None):
    
    iterator = core.G3File(filepath)
    
    while True:
        new_frame = iterator.next()
        if "CoaddedObservationIDs" in new_frame.keys():
            database = new_frame
            break
    
    contents  = "\n"
    if print_frame:
        contents += "# -------------- #\n"
        contents += "#  Keys stored:  #\n"
        contents += "# -------------- #\n"
        contents += "\n"
        keys = list(database.keys())
        keys = "\n".join(sorted(keys))
        contents += keys
        contents += "\n"
    
    contents += "\n\n"
    contents += "# ------------------------ #\n"
    contents += "#  Requested information:  #\n"
    contents += "# ------------------------ #\n"
    contents += "\n\n"
    
    all_oids = database["CoaddedObservationIDs"]
    
    quantity_names  = "   ".join([entry[2].ljust(15) for entry in info_tuples])
    table_first_row = "Observation ID   {}\n".format(quantity_names)
    
    for sub_field, oids in all_oids.items():
        contents += "* {}\n".format(sub_field)
        contents += "---------------\n"
        contents += "\n"
        
        contents += table_first_row
        contents += "\n"
        
        for oid in oids:
            data = []
            for entry in info_tuples:
                key   = entry[0]
                units = entry[1]
                data.append(database[key][sub_field][str(oid)] / units)
            data = "   ".join(["{:+15.2e}".format(number) for number in data])
            
            contents += "{:>15}   {}\n".format(oid, data)
        
        contents += "\n\n"
    
    print(contents)
    
    if output_path is not None:
        with open(output_path, "w") as f_obj:
            f_obj.writelines(contents)


# ==============================================================================




# ==============================================================================
# Construct and start the pipeline
# ------------------------------------------------------------------------------


def run(input_files=[], min_file_size=0.01, output_file='coadded_maps.g3',
        map_ids=["90GHz"], sub_fields=["ra0hdec-44.75"],
        min_obs_id=0, max_obs_id=999999999, bad_obs_ids=[],
        temperature_maps_only=False,
        trick_pipeline_to_receive_left_right_maps=False,
        ids_to_append_directions_to=[],
        maps_split_by_scan_direction=False,
        combine_left_right=False, combine_different_wafers=False,
        subtract_existing_maps=False,
        collect_averages_from_flagging_statistics=False,
        calculate_pW_to_K_conversion_factors=False,
        calculate_calibrator_response_vs_el=False,
        calibration_data_dir=None,
        bolo_timestreams_dir=None,
        calculate_map_rmss_and_weight_stats=False,
        rmss_and_wgts_from_coadds_or_individuals=None,
        rmss_and_wgts_from_signals_or_noises=None,
        calculate_noise_from_individual_maps=False,
        calculate_noise_from_coadded_maps=False,
        point_source_list_file=None,
        calculate_cross_spectra_with_planck_map=False,
        planck_map_fits_file=None,
        auxiliary_files_directory=None,
        calculate_pointing_discrepancies=False,
        logger_name="", less_verbose=False,
        bad_map_list_file=None, log_file=None):
    
    
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
    
    
    def is_a_good_obs_id(g3_file, min_obs_id, max_obs_id, bad_obs_ids=[]):
    
        meaningful_part_of_path = g3_file.split('/')[-1].split('.')[0]
        try:
            obs_id = int(meaningful_part_of_path.split('_')[0])
        except ValueError:
            # * Probably a g3 file that contains coadded maps is encountered
            #   if the name cannot be converted to an integer.
            return True
    
        if (obs_id >= min_obs_id) and \
           (obs_id <= max_obs_id) and \
           (obs_id not in bad_obs_ids):
            return True
        else:
            return False
    
    
    all_good_g3_files = []
    for g3_file in input_files:
        if os.path.isfile(g3_file) and \
           is_a_good_obs_id(g3_file, min_obs_id, max_obs_id, bad_obs_ids) and \
           (os.path.getsize(g3_file) >= min_file_size*2**30):
            all_good_g3_files.append(g3_file)
    
    if len(all_good_g3_files) == 0:
        log("")
        log("No applicable input files, so there is nothing to do!")
        log("")
        return 0
    if len(all_good_g3_files) == 1:
        output_file_another_name = output_file.replace(".g4", ".g3")
        # * A situation very specific to how cache_field_maps.py
        #   handle the input and output files.
        if output_file_another_name == all_good_g3_files[0]:
            log("")
            log("The output file will just be the same as the input, ")
            log("so there is nothing to do!")
            log("")
            return 0
    
    
    # - Show settings related to I/O
    
    log("")
    log("# ============================ #")
    log("#  Start making coadded maps!  #")
    log("# ============================ #")
    log("")
    
    log("- Input files to supply:")
    log("")
    for input_file in all_good_g3_files:
        log(input_file)
    log("\n")
    
    log("- File to which the coadded maps will be saved:")
    log("")
    log(output_file)
    log("\n")
    
    log("- The range of observation ID of interest:")
    log("")
    log("from %s to %s", min_obs_id, max_obs_id)
    log("\n")
    
    log("- Observation IDs that were excluded:")
    log("  (IDs that were either added or ignored previously)")
    log("")
    log(bad_obs_ids)
    log("")
    
    
    # - Construct a pipeline that analyzes and adds maps
    
    pipeline = core.G3Pipeline()
    
    for i, f in enumerate(all_good_g3_files):
        f = maybe_get_gfal_copied_file(
                f, outdir=os.path.dirname(output_file))
        all_good_g3_files[i] = f
    
    pipeline.Add(core.G3Reader,
                 filename=all_good_g3_files)
    
    pipeline.Add(drop_duplicate_analysis_results)
    
    if trick_pipeline_to_receive_left_right_maps:
        pipeline.Add(FlagBadMaps,
                     map_ids=ids_to_append_directions_to,
                     t_only=temperature_maps_only,
                     point_source_list_file=point_source_list_file,
                     auxiliary_files_directory=auxiliary_files_directory,
                     bad_map_list_file=bad_map_list_file,
                     logging_function=log)
        pipeline.Add(AppendDirectionsToMapIDs,
                     sub_fields=sub_fields,
                     ids_to_modify=ids_to_append_directions_to,
                     logging_function=log)
    
    def print_frame(frame):
        log(""); log(frame); log("")
    pipeline.Add(print_frame)
    
    pipeline.Add(AnalyzeAndCoaddMaps,
                 map_ids=map_ids, map_sources=sub_fields,
                 t_only=temperature_maps_only,
                 maps_split_by_scan_direction=maps_split_by_scan_direction,
                 combine_left_right=combine_left_right,
                 combine_different_wafers=combine_different_wafers,
                 allow_subtraction=subtract_existing_maps,
                 collect_averages_from_flagging_stats=\
                     collect_averages_from_flagging_statistics,
                 calculate_pW_to_K_conversion_factors=\
                     calculate_pW_to_K_conversion_factors,
                 calculate_calibrator_response_vs_el=\
                     calculate_calibrator_response_vs_el,
                 calibration_data_dir=calibration_data_dir,
                 bolo_timestreams_dir=bolo_timestreams_dir,
                 min_field_obs_id=min_obs_id,
                 max_field_obs_id=max_obs_id,
                 calculate_map_rmss_and_weight_stats=\
                     calculate_map_rmss_and_weight_stats,
                 rmss_and_wgts_from_coadds_or_individuals=\
                     rmss_and_wgts_from_coadds_or_individuals,
                 rmss_and_wgts_from_signals_or_noises=\
                     rmss_and_wgts_from_signals_or_noises,
                 calculate_pointing_discrepancies=\
                     calculate_pointing_discrepancies,
                 calculate_noise_from_individual_maps=\
                     calculate_noise_from_individual_maps,
                 calculate_noise_from_coadded_maps=\
                     calculate_noise_from_coadded_maps,
                 point_source_list_file=point_source_list_file,
                 calculate_cross_spectra_with_planck_map=\
                     calculate_cross_spectra_with_planck_map,
                 planck_map_fits_file=planck_map_fits_file,
                 auxiliary_files_directory=auxiliary_files_directory,
                 logging_function=log,
                 less_verbose=less_verbose)
    
    pipeline.Add(lambda frame: "AreCoaddedMapsContained" in frame)
    
    pipeline.Add(core.G3Writer,
                 filename=output_file)
    
    
    # - Time to start!
    
    log("\n")
    log("# -------------------------- #")
    log("#  Starting the pipeline...  #")
    log("# -------------------------- #")
    log("\n")
    log("Every frame passed to the pipeline will be printed.")
    log("\n")
    
    if log_file is None:
        profile = True
    else:
        profile = False
    
    pipeline.Run(profile=profile)
    
    log("")
    
    for f in all_good_g3_files:
        f = maybe_del_gfal_copied_file(f)
    
    return 1


# ==============================================================================





# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                 description="This script makes coadded maps of the "
                             "1500 sq. deg. field and can also perform some "
                             "analyses on the maps such as calculating noise.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-I", "--input_files",
                        type=str, action="store", nargs="+",
                        help="The paths to the g3 files "
                             "that contain maps to be analyzed and added. ")
    
    parser.add_argument("-m", "--min_file_size",
                        type=float, action="store", default=0.01,
                        help="The minimum size (in GB) a g3 file needs "
                             "to have for it to be considered as a good file. "
                             "This is to reduce the occurrence of a situation "
                             "where the script crashes due to reading a "
                             "problematic file.")
    
    parser.add_argument("-O", "--output_file",
                        type=str, action="store", default="coadded_maps.g3",
                        help="The path of the output g3 file "
                             "that will contain the coadded maps and "
                             "quantities obtained during the analyses.")
    
    parser.add_argument("-i", "--map_ids",
                        type=str, action="store", nargs="+",
                        default=["90GHz"],
                        help="Only those map frames whose 'Id' store values "
                             "that are in this list will be analyzed and "
                             "added. If the list is empty, then the bands "
                             "argument above will be used.")
    
    parser.add_argument("-S", "--sub_fields",
                        type=str, action="store", nargs="+",
                        default=["ra0hdec-44.75"],
                        help="The sub-field(s) of which we are interested in "
                             "analyzing and adding maps.")
    
    parser.add_argument("-s", "--min_obs_id",
                        type=int, action="store", default=0,
                        help="The smallest observation ID that will be "
                             "considered to be used in the coadded maps.")
    
    parser.add_argument("-l", "--max_obs_id",
                        type=int, action="store", default=999999999,
                        help="The largest observation ID that will be "
                             "considered to be used in the coadded map.")
    
    parser.add_argument("-b", "--bad_obs_ids",
                        type=int, action="store", nargs="+", default=[],
                        help="The observation IDs that will be excluded from "
                             "making the coadded maps. The script has some "
                             "simple cuts to decide whether to include "
                             "certain maps, but one can manually "
                             "specify observations to ignore, too. Despite "
                             "the name 'bad', these maps may not necessrily"
                             "be the ones that have low data quality; they can "
                             "be the ones that have already been added to the "
                             "coadded maps and thus don't need to be processed "
                             "again.")
    
    parser.add_argument("-t", "--temperature_maps_only",
                        action="store_true", default=False,
                        help="Whether the maps are temperature-only "
                             "maps or not.")
    
    parser.add_argument("-P", "--trick_pipeline_to_receive_left_right_maps",
                        action="store_true", default=False,
                        help="Whether to trick the pipeline module that "
                             "co-adds maps and performs analyses into "
                             "receiving maps whose IDs imply that "
                             "they are made from only one scan direction. "
                             "This option is useful in the case where "
                             "we want to calculate noise levels from "
                             "coadded maps but the individual maps"
                             "contain data from both scan directions, in which "
                             "case we can pretend that one half of the maps "
                             "were 'left-going' ones and the other half "
                             "'right-going' ones.")
    
    parser.add_argument("-d", "--ids_to_append_directions_to",
                        type=str, action="store", nargs="+", default=[],
                        help="The original map IDs to which fake directions "
                             "will be appended.") 
    
    parser.add_argument("-r", "--maps_split_by_scan_direction",
                        action="store_true", default=False,
                        help="Whether the maps stored in the data frames are "
                             "split by the scan direction, including the case "
                             "explained above, where we pretend the maps are "
                             "split by the direction.")
    
    parser.add_argument("-c", "--combine_left_right",
                        action="store_true", default=False,
                        help="Whether to add maps made from left-going scans "
                             "and those made from right-going scans together. "
                             "If True, a map whose 'Id' is 'Left'+map_id and "
                             "one whose 'Id' is 'Right'+map_id will be added, "
                             "where map_id is an entry specified in the list "
                             "map_ids above.")
    
    parser.add_argument("-C", "--combine_different_wafers",
                        action="store_true", default=False,
                        help="Whether to add maps made from different wafers "
                             "together. If True, maps whose 'Id's look like "
                             "map_id+a_wafer_name will be added, where"
                             "map_id is an entry specified in the list "
                             "map_ids above.")
    
    parser.add_argument("-u", "--subtract_existing_maps",
                        action="store_true", default=False,
                        help="Whether to subtract certain observation's map "
                             "from the coadded maps when the former is not "
                             "wanted anymore.")
    
    parser.add_argument("-a", "--collect_averages_from_flagging_statistics",
                        action="store_true", default=False,
                        help="Whether to calculate average number of "
                             "bolometers that were not flagged "
                             "(average over all scans of an observation), "
                             "and average number of bolometers that "
                             "were flagged by each reason.")
    
    parser.add_argument("-k", "--calculate_pW_to_K_conversion_factors",
                        action="store_true", default=False,
                        help="Whether to calculate the mean value of the "
                             "pW/K absolute calibration factors for each band "
                             "of each wafer.")
    
    parser.add_argument("-e", "--calculate_calibrator_response_vs_el",
                        action="store_true", default=False,
                        help="Whether to calculate a fractional change of "
                             "detectors' response to the calibrator at the "
                             "top of a field with respect to the bottom.")
    
    parser.add_argument("-L", "--calibration_data_dir",
                        type=str, action="store", default=None,
                        help="Path to the directory that contain "
                             "auto-processed calibration results, which should "
                             "contain a directory labelled 'calibrator', "
                             "where CalibratorResponse data will be obtained")
    
    parser.add_argument("-T", "--bolo_timestreams_dir",
                        type=str, action="store", default=None,
                        help="Path to the directory that contain "
                             "detectors' raw timestreams, which should contain "
                             "a directory labelled 'calibrator', where "
                             "elevations at which some calibrator observations "
                             "were taken will be checked.")
    
    parser.add_argument("-w", "--calculate_map_rmss_and_weight_stats",
                        action="store_true", default=False,
                        help="Whether to calculate an RMS (standard deviation) "
                             "from a temperature map, a mean from the TT "
                             "weight map, and the number of pixels in the "
                             "weight map whose weights are above certain "
                             "threshold.")
    
    parser.add_argument("-y", "--rmss_and_wgts_from_coadds_or_individuals",
                        type=str, action="store",
                        choices=["c", "i", "ci"], default="",
                        help="Whether to calculate the statistics from "
                             "individual maps, coadded maps, or both.")
    
    parser.add_argument("-Y", "--rmss_and_wgts_from_signals_or_noises",
                        type=str, action="store",
                        choices=["s", "n", "sn"], default="",
                        help="Whether to calculate the statistics from "
                             "signal maps, noise maps, or both. 'Signal maps' "
                             "just refer to either individual observations'"
                             "maps or coadded maps, and 'noise maps' refer to "
                             "'left-going' maps minus 'right-going' maps.")
    
    parser.add_argument("-p", "--calculate_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to calculate the difference between the "
                             "true positions of some point sources and their "
                             "measured positions.")
    
    parser.add_argument("-x", "--calculate_cross_spectra_with_planck_map",
                        action="store_true", default=False,
                        help="Whether to calculate the cross spectra between "
                             "individual maps and a Planck map.")
    
    parser.add_argument("-X", "--planck_map_fits_file",
                        type=str, action="store", default=None,
                        help="Path to FITS files that contain Planck maps. "
                             "The path here actually does not point to any "
                             "actual file. Instead, the path should contain "
                             "the words 'BAND' and 'MISSION' in it so that "
                             "the paths of the actual files can be obtained "
                             "by replacing 'BAND' with '100GHz', '143GHz', "
                             "or '217GHz' and 'MISSION' with 'fullmission', "
                             "halfmission-1', and halfmission-2'.")
    
    parser.add_argument("-n", "--calculate_noise_from_individual_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels from "
                             "each map that goes into the coadded maps.")
    
    parser.add_argument("-N", "--calculate_noise_from_coadded_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels from the "
                             "running coadded maps. In other words, the noise "
                             "is calculated from the coadded 'left-going' maps "
                             "minus the coadded 'right-going' maps whenever "
                             "the same number of maps have been added to both "
                             "of the coadds.")
    
    parser.add_argument("-M", "--point_source_list_file",
                        type=str, action="store", default="",
                        help="Path to a point source list, which will be used "
                             "when making a mask for calcuting power spectra "
                             "and map rms values.")
    
    parser.add_argument("-A", "--auxiliary_files_directory",
                        type=str, action="store", default=None,
                        help="The path to the directory that contains "
                             "auxiliary data such as masks for calculating "
                             "power spectra. If the appropriate files exist, "
                             "then the analysis module can just load them "
                             "instead of creating them.")
    
    parser.add_argument("-g", "--logger_name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "
                             "record log messages.")
    
    parser.add_argument("-V", "--less_verbose",
                        action="store_true", default=False,
                        help="Whether to make the logger less verbose "
                             "by turning off some messages. Search '$$$' "
                             "to see the nature of these messages.")
    
    parser.add_argument("-B", "--bad_map_list_file",
                        type=str, action="store", default="",
                        help="A text file in which maps that are found to be "
                             "bad during the coadding processes are recorded. "
                             "In other words, the script does not read a list "
                             "of observation IDs from the file and exclude "
                             "them from the analyses and coadding; rather, the "
                             "script only records maps that were found to be "
                             "bad when checking them as part of the analyses.")
    
    parser.add_argument("-G", "--log_file",
                        type=str, action="store", default=None,
                        help="The file to which the logger will send messages.")
    
    arguments = parser.parse_args()
    
    run(**vars(arguments))


# ==============================================================================

