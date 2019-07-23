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
#        each wafer.                                                           #
#                                                                              #
#    - Comparing SPT maps with Planck maps                                     #
#        The script can compare an individual SPT map with a Planck map by     #
#        calculating the cross spectrum between the two, calculating the       #
#        autospectrum of the Planck map, dividing the first spectrum by the    #
#        second one, and calculating the average of the ratios in the ell      #
#        range [750, 1500], i.e., average of                                   #
#        Cl_(SPT x Planck) / Cl_(Planck_x_Planck).                             #
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
#        spectra. In order to save time and memory, the 5 degree by 5 degree   #
#        patches mentioned above are used.                                     #
#                                                                              #
#    - Calculating pointing discrepancies                                      #
#        The script can calculate the differences between the true positions   #
#        and the measured ones of three bright point sources in each of the    #
#        four sub-fields. The measured positions are results of fitting a      #
#        2D gaussian to a small (4 arcmin. by 4 arcmin.) region of a map       #
#        that is roughly centered on each point source.                        #
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
#  The first part of this script consists of definitions of functions that     #
#  are relevent to performing these analyses, and the second part defines a    #
#  giant (~1200 lines...) pipeline module (AnalyzeAndCoaddMaps) that           #
#  adds maps together and performs the analyses.                               #
#                                                                              #
# ============================================================================ #


import argparse
import datetime
import logging
import os
import sys
import glob
import gc
import numpy

from spt3g import core
from spt3g import mapmaker
from spt3g import util
from spt3g import coordinateutils
from spt3g import std_processing
from spt3g import mapspectra
from scipy import signal
from spt3g.mapspectra import map_analysis



# ==============================================================================
# Define functions needed to perform the analyses
# ------------------------------------------------------------------------------


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
            t_map_one = mapmaker.mapmakerutils.remove_weight_t(
                            frame_one["T"], frame_one["Wunpol"])
            t_map_two = mapmaker.mapmakerutils.remove_weight_t(
                            frame_two["T"], frame_two["Wunpol"])
        else:
            t_map_one = frame_one["T"]
            t_map_two = frame_two["T"]
        
        logfun("$$$ t1 after potential weights removal : %+7.3e"
               %(numpy.asarray(t_map_one)[iy][ix]/tu))
        logfun("$$$ t2 after potential weights removal : %+7.3e"
               %(numpy.asarray(t_map_two)[iy][ix]/tu))
        
        # - Then, add field maps (and weight maps if applicable) together
        
        t_map = t_map_one + t_map_two * multiply_2nd_maps_by
        del t_map_one, t_map_two
        gc.collect()
        w_map  = frame_one["Wunpol"].TT + \
                 frame_two["Wunpol"].TT * multiply_2nd_weights_by
        w_maps = core.G3SkyMapWeights()
        w_maps.TT = w_map
        del w_map
        gc.collect()
        
        logfun("$$$ t resultant : %+7.3e"
               %(numpy.asarray(t_map)[iy][ix]/tu))
        logfun("$$$ w resultant : %+7.3e [1/mK^2]"
               %(numpy.asarray(w_maps.TT)[iy][ix]*tu*tu))
        
        # - After that, remove weights if necessary
        
        if remove_weights_afterward:
            t_map  = mapmaker.mapmakerutils.remove_weight_t(
                         t_map, w_maps)
        
        logfun("$$$ t after potential weights removal : %+7.3e [mK]"
               %(numpy.asarray(t_map)[iy][ix]))
        
        # - Furthur divide the map by some number if requested
        
        if divide_result_by != 1.0:
            t_map = t_map / divide_result_by
        
        logfun("$$$ t after division : %+7.3e [mK]"
               %(numpy.asarray(t_map)[iy][ix]/tu))
        
        # - Finally, record the results in a new frame
        
        new_frame = core.G3Frame(core.G3FrameType.Map)
        new_frame["T"] = t_map
        if record_weights:
            new_frame["Wunpol"] = w_maps
        
        logfun("$$$ t final : %+7.3e"
               %(numpy.asarray(new_frame["T"])[iy][ix]/tu))
        try:
            logfun("$$$ w final : %+7.3e [1/mK^2]"
                   %(numpy.asarray(new_frame["Wunpol"].TT)[iy][ix]*tu*tu))
        except:
            logfun("$$$ w final : N/A")
    
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
        map_frame, center_ra, center_dec, t_only=True):
    
    new_map_frame = core.G3Frame(core.G3FrameType.Map)

    if t_only:
        
        map_parameters = {"x_len": int(5.0*core.G3Units.deg/map_frame["T"].res),
                          "y_len": int(5.0*core.G3Units.deg/map_frame["T"].res),
                          "res"  : map_frame["T"].res,
                          "proj" : map_frame["T"].proj,
                          "alpha_center": center_ra,
                          "delta_center": center_dec,
                          "coord_ref"   : map_frame["T"].coord_ref,
                          "pol_type"    : map_frame["T"].pol_type}
        
        smaller_field_map = coordinateutils.FlatSkyMap(**map_parameters)
        coordinateutils.reproj_map(map_frame["T"], smaller_field_map, 1)
        new_map_frame["T"] = smaller_field_map
        
        if "Wunpol" in map_frame.keys():
            smaller_weight_map = coordinateutils.FlatSkyMap(**map_parameters)
            coordinateutils.reproj_map(
                map_frame["Wunpol"].TT, smaller_weight_map, 1)
            new_map_frame["Wunpol"] = core.G3SkyMapWeights()
            new_map_frame["Wunpol"].TT = smaller_weight_map
    
    return new_map_frame




def get_band_average(bolo_names_from_each_scan, bolo_props_map):
    
    # * This one is called by the function below
    
    n_bolos_from_each_scan = {"90": [], "150": [], "220": []}
    bolos_with_properties  = set(bolo_props_map.keys())
    
    for bolo_names_from_one_scan in bolo_names_from_each_scan.values():
        n_bolos_from_one_scan = {"90": 0, "150": 0, "220": 0}
        identifiable_bolos = set(bolo_names_from_one_scan) & \
                             bolos_with_properties
        for bn in identifiable_bolos:
            band = bolo_props_map[bn].band / core.G3Units.GHz
            if not numpy.isfinite(band):
                continue
            band = str(int(band))
            if band in n_bolos_from_one_scan.keys():
                n_bolos_from_one_scan[band] += 1
        
        for band, n_bolos in n_bolos_from_one_scan.items():
            n_bolos_from_each_scan[band].append(n_bolos)
    
    return {band: numpy.mean(ns_bolos) \
            for band, ns_bolos in n_bolos_from_each_scan.items()}




def collect_averages_from_flagging_info(
        pipe_info_frame, bolo_props_map, recognized_reasons):
    
    averages = {"90": {}, "150": {}, "220": {}}
    
    # - Collect average numbers of bolos not flagged
    
    bolos_not_flagged_each_scan = pipe_info_frame["SurvivingBolos"]
    # * This is a core.G3MapVectorString instance
    
    total_n_scans = numpy.max([len(scans) for reason, scans \
                               in pipe_info_frame["DroppedBolos"].items()])
    n_recorded_here = len(bolos_not_flagged_each_scan.keys())
    
    diff_n_scans = total_n_scans - n_recorded_here
    for i in range(diff_n_scans):
        bolos_not_flagged_each_scan["Dummy"+str(i)] = core.G3VectorString()
    
    avg_n_okay_bolos_each_band = \
        get_band_average(bolos_not_flagged_each_scan, bolo_props_map)
    for band, average in avg_n_okay_bolos_each_band.items():
        averages[band]["TotalNotFlagged"] = average
    
    # - Collect average numbers of bolos flagged for each reason
    
    # -- Gather the numbers from all reasons
    
    all_bolos_flagged_each_scan = pipe_info_frame["DroppedBolos"]
    # * This is a core.G3MapVectorVectorString instance
    
    for flagging_reason, bolos_flagged_each_scan \
    in  all_bolos_flagged_each_scan.items():
        bl_flgd_ea_sc_reorganized = core.G3MapVectorString()
        for scan_number, list_of_bolos in enumerate(bolos_flagged_each_scan):
            bl_flgd_ea_sc_reorganized[str(scan_number)] = list_of_bolos
        avg_bad_bolos_each_band = \
            get_band_average(bl_flgd_ea_sc_reorganized, bolo_props_map)
        for band, average in avg_bad_bolos_each_band.items():
            averages[band][flagging_reason] = average
    
    # -- Collectively call unrecognized reasons "Others"
    
    for band, reasons_and_numbers in averages.items():
        averages[band]["Others"] = 0.0
        unrecognized_reasons = []
        for reason, number in reasons_and_numbers.items():
            if reason not in recognized_reasons:
                averages[band]["Others"] += number
                unrecognized_reasons.append(reason)
        for reason in unrecognized_reasons:
            averages[band].pop(reason)
        for reason in recognized_reasons:
            if reason not in reasons_and_numbers.keys():
                averages[band][reason] = 0.0
    
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
    avg_n_removed_bolos_each_band = \
        get_band_average(bolos_removed_each_scan, bolo_props_map)
    for band, average in avg_n_removed_bolos_each_band.items():
        averages[band]["TotalRemoved"] = average
    
    # * "TotalNotFlagged", "TotalRemoved", and "Others" keys
    #   need to be changed accordingly if the corresponding names
    #   in the initialization of AnalyzeAndCoaddMaps are changed.
    
    return averages




def collect_medians_of_pW_per_K_factors(
        wafers_of_interest, band_of_interest, calframe, flagging_stats=None,
        logfun=print):
    
    bpm = calframe["BolometerProperties"]
    flu = core.G3Units.K * core.G3Units.rad * core.G3Units.rad
    flx = {"RCW38": { "90": 4.0549662e-07*flu,
                     "150": 2.5601153e-07*flu,
                     "220": 2.8025804e-07*flu},
           "MAT5A": { "90": 2.5738063e-07*flu,
                     "150": 1.7319235e-07*flu,
                     "220": 2.1451640e-07*flu}}
           # * Copied from spt3g_software/calibration/python/apply_t_cal.py
    
    values_by_wafer = \
         {wafer: {"PicowattsPerKelvin": [],
                  "CalibratorResponse": [],
                  "HIIFluxCalibration": [],
                  "HIIIntegratedFluux": [],
                  "HIISkyTransmission": []} \
          for wafer in wafers_of_interest}
    
    for key in calframe.keys():
        if "FluxCalibration" in key:
            flux_calibration_key = key
            if "RCW38" in flux_calibration_key:
                source = "RCW38"
            if "MAT5A" in flux_calibration_key:
                source = "MAT5A"
        if "IntegralFlux" in key:
            integral_flux_key = key
        if "SkyTransmission" in key:
            sky_transmission_key = key
        
    available_bolos = set(calframe["BolometerProperties"].keys()) & \
                      set(calframe["CalibratorResponse"].keys())  & \
                      set(calframe[flux_calibration_key].keys())  & \
                      set(calframe[integral_flux_key].keys())
    
    logfun("$$$ Recording calibration chain-related quantities ...")
    logfun("$$$ n available bolos from calframe : %d" %(len(available_bolos)))
    
    if flagging_stats is not None:
        bolos_used_each_scan = flagging_stats["SurvivingBolos"]
        bolos_used_at_least_once = set([])
        for scan_number, bolos_used_that_scan \
        in  bolos_used_each_scan.iteritems():
            bolos_used_at_least_once = \
                bolos_used_at_least_once | set(bolos_used_that_scan)
        available_bolos = available_bolos & bolos_used_at_least_once
    
    logfun("$$$ n_available bolos after using flagging stats : %d"
          %(len(available_bolos)))
    
    for bolo in available_bolos:
        band  = bpm[bolo].band / core.G3Units.GHz
        if not numpy.isfinite(band):
            continue
        band  = str(int(band))
        wafer = bpm[bolo].wafer_id.upper() 
        if (band != band_of_interest) or (wafer not in wafers_of_interest):
            continue
        
        cr = calframe["CalibratorResponse"][bolo]
        fc = calframe[flux_calibration_key][bolo]
        fi = calframe[integral_flux_key][bolo]
        if band in calframe[sky_transmission_key].keys():
            st = calframe[sky_transmission_key][band]
            pk = cr * fc * fi * st / flx[source][band]
        else:
            st = numpy.nan
            pk = numpy.nan
        
        values_by_wafer[wafer]["PicowattsPerKelvin"].append(pk)
        values_by_wafer[wafer]["CalibratorResponse"].append(cr)
        values_by_wafer[wafer]["HIIFluxCalibration"].append(fc)
        values_by_wafer[wafer]["HIIIntegratedFluux"].append(fi)
        values_by_wafer[wafer]["HIISkyTransmission"].append(st)
    
    for waf in values_by_wafer.keys():
        for quantity, vals in values_by_wafer[waf].items():
            try:
                
                if quantity == "PicowattsPerKelvin":
                    logfun("$$$ [%s] number of collected pW/K factors : %d"
                           %(waf, len(vals)))
                
                median = numpy.nanmedian(vals)
                
                if quantity == "PicowattsPerKelvin":
                    u = core.G3Units.pW / core.G3Units.K
                    logfun("$$$        median of pW/K : %7.3e" %(median/u))
            
            except:
                median = numpy.nan
            values_by_wafer[waf][quantity] = median
    
    logfun("$$$ ... that's all I wanted to check.")
    
    return values_by_wafer




def identify_pixels_of_non_atypical_region(
        map_frame, center_dec, point_source_list_file):
    
    ras, decs = coordinateutils.maputils.get_ra_dec_map(map_frame["T"])
    ras  = numpy.asarray(ras/core.G3Units.deg)
    decs = numpy.asarray(decs/core.G3Units.deg)
    center_dec /= core.G3Units.deg
    
    ptsrc_msk = numpy.asarray(mapspectra.apodmask.makeApodizedPointSourceMask(
                                  map_frame, point_source_list_file,
                                  apod_type="none", zero_border_arcmin=0.0))
        
    typical_pixels = numpy.where((ras > -48.0) & (ras < 48.0) &
                                 (decs < (center_dec + 3.0))  & 
                                 (decs > (center_dec - 3.0))  &
                                 (ptsrc_msk == 1.0))
    
    del ras, decs, ptsrc_msk
    gc.collect()
    
    """final_mask = numpy.zeros(map_frame["T"].shape)
    final_mask[typical_pixels] = 1.0
    from matplotlib import pyplot
    pyplot.imshow(final_mask, cmap="gray")
    pyplot.show()"""
    
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
        assert (not map_frame[fmk].is_weighted), fmk + " map is still weighted!"
    
    if t_only:
        
        wu = 1 / (core.G3Units.mK * core.G3Units.mK)
        nominal_weights_dict = \
            { "90": {"ra0hdec-44.75": 30.0 * wu,
                     "ra0hdec-52.25": 34.0 * wu,
                     "ra0hdec-59.75": 36.0 * wu,
                     "ra0hdec-67.25": 50.0 * wu},
             "150": {"ra0hdec-44.75": 45.0 * wu,
                     "ra0hdec-52.25": 51.0 * wu,
                     "ra0hdec-59.75": 54.0 * wu,
                     "ra0hdec-67.25": 75.0 * wu},
             "220": {"ra0hdec-44.75":  3.5 * wu,
                     "ra0hdec-52.25":  3.8 * wu,
                     "ra0hdec-59.75":  4.5 * wu,
                     "ra0hdec-67.25":  6.5 * wu}}
        
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
            better_t_vals = t_vals[pixs]
            if w_vals is None:
                better_w_vals = None
            else:
                better_w_vals = w_vals[pixs]
        
        map_stdddev = numpy.std(better_t_vals)
        if better_w_vals is None:
            mean_weight = numpy.nan
            n_pix_g_wgt = numpy.nan
        else:
            mean_weight = numpy.mean(better_w_vals)
            nominal_wgt = nominal_weights_dict[band][sub_field]
            n_pix_g_wgt = len(numpy.where(better_w_vals > nominal_wgt)[0])
        
        logfun("$$$ number of pixels with good weights : %s" %(n_pix_g_wgt))
        
        fluctuation_metrics["TMapStandardDeviations"] = map_stdddev
        fluctuation_metrics["MeansOfTTWeights"]       = mean_weight
        fluctuation_metrics["NumbersOfPixelsWithGoodTTWeights"] = n_pix_g_wgt
    
    logfun("$$$ ... that's all I wanted to know.")
    
    return fluctuation_metrics




def calculate_pointing_discrepancies(
        map_frame, sub_field, map_stddev=None):
    
    assert (not map_frame["T"].is_weighted), \
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
             {"1st": numpy.array([352.32358, -47.50531, 1408.089]),
              "2nd": numpy.array([314.06833, -47.24664, 1374.576]),
              "3rd": numpy.array([ 41.50037, -46.85467,  715.329])},
         "ra0hdec-52.25": \
             {"1st": numpy.array([ 32.69288, -51.01703, 3819.563]),
              "2nd": numpy.array([ 23.27404, -52.00094, 1026.665]),
              "3rd": numpy.array([359.47312, -53.18686,  864.243])},
         "ra0hdec-59.75": \
             {"1st": numpy.array([ 47.48363, -60.97761,  869.843]),
              "2nd": numpy.array([ 45.96104, -62.19042,  832.717]),
              "3rd": numpy.array([ 14.69433, -56.98650,  785.528])},
         "ra0hdec-67.25": \
             {"1st": numpy.array([329.27542, -69.68981, 1114.524]),
              "2nd": numpy.array([337.25092, -69.17492,  445.331]),
              "3rd": numpy.array([325.44375, -64.18742,  388.218])}}
    
    discrep_dict = {"1st": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan},
                    "2nd": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan},
                    "3rd": {"DeltaRa": numpy.nan, "DeltaDec": numpy.nan}}
    flux_dict = {"1st": numpy.nan, "2nd": numpy.nan, "3rd": numpy.nan}
    snr_dict  = {"1st": numpy.nan, "2nd": numpy.nan, "3rd": numpy.nan}
    
    
    relevant_point_sources = some_brightest_sources[sub_field]
    
    for point_source_rank, info in relevant_point_sources.items():
        
        true_right_ascension = info[0]
        if true_right_ascension > 180:
            true_right_ascension -= 360
        true_right_ascension *= usys.deg
        
        true_declination  = info[1]
        true_declination *= usys.deg
        
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
        
        measured_right_ascension, measured_declination = \
            map_frame["T"].xy_to_angle(measured_x, measured_y)
        
        delta_ra  = measured_right_ascension - true_right_ascension
        delta_dec = measured_declination     - true_declination
        effective_delta_ra = delta_ra * numpy.cos(true_declination/usys.rad)
        
        discrep_dict[point_source_rank]["DeltaRa"]  = effective_delta_ra
        discrep_dict[point_source_rank]["DeltaDec"] = delta_dec
        
        nearest_x = int(numpy.round(measured_x))
        nearest_y = int(numpy.round(measured_y))
        source_centered_mini_map = \
            numpy.asarray(map_frame["T"])\
                [nearest_y-half_y_length:nearest_y+half_y_length+1,
                 nearest_x-half_x_length:nearest_x+half_x_length+1]
        source_flux = numpy.sum(source_centered_mini_map) * \
                      map_frame["T"].x_res * map_frame["T"].y_res
        flux_dict[point_source_rank] = source_flux
        
        if map_stddev is None:
            source_snr = numpy.nan
        else:
            source_snr = numpy.max(source_centered_mini_map) / map_stddev
        snr_dict[point_source_rank] = source_snr
    
    return discrep_dict, flux_dict, snr_dict




def create_mask_for_powspec_calc_of_small_region(map_frame, point_source_file):
    
    ptsrc_mask  = mapspectra.apodmask.makeApodizedPointSourceMask(
                    map_frame, point_source_file)
    
    edge_mask = coordinateutils.FlatSkyMap(
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




def create_mini_planck_map_frame(
        fits_file, spt_map_frame, center_ra, center_dec, t_only=True):
    
    planck_healpix_maps = mapmaker.load_spt3g_map(fits_file)
    
    if t_only:
        planck_flatsky_map = coordinateutils.maputils.healpix_to_flatsky(
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
                            new_map_frame, center_ra, center_dec)
    
    del new_map_frame
    gc.collect()
    
    return map_frame_smaller




def calculate_average_ratio_of_spt_planck_xspectra(
        spt_map_frame, planck_map_frame, mask, t_only=True):
    
    assert (not spt_map_frame["T"].is_weighted), \
           "The SPT T map is still weighted!"
    assert (not planck_map_frame["T"].is_weighted), \
           "The Planck T map is still weighted!"
    
    average_ratios = {}
    
    ell_lo = 750
    ell_hi = 1250
    
    if t_only:
        planck_x_planck = map_analysis.calculate_powerspectra(
                              planck_map_frame, input2=planck_map_frame,
                              delta_l=50, l_min=300, l_max=6000,
                              apod_mask=mask, realimag="real", flatten=False)
        
        spt_x_planck = map_analysis.calculate_powerspectra(
                           spt_map_frame, input2=planck_map_frame,
                           delta_l=50, l_min=300, l_max=6000,
                           apod_mask=mask, realimag="real", flatten=False)
        
        cl_ratios   = spt_x_planck["TT"] / planck_x_planck["TT"]
        bin_centers = planck_x_planck["TT"].bin_centers
        idx = numpy.where((bin_centers > ell_lo) & (bin_centers < ell_hi))[0]
        average_ratio_tt = numpy.mean(cl_ratios[idx])
        
        average_ratios["TT"] = average_ratio_tt
    
    return average_ratios




def calculate_noise_levels(map_frame, mask, t_only=True):
    
    assert (not map_frame["T"].is_weighted), \
            "The T maps is still weighted!"
    
    noises = {}
    
    ell_lo = 3000
    ell_hi = 5000
    
    if t_only:
        cls = map_analysis.calculate_powerspectra(
                  map_frame, input2=None,
                  delta_l=50, l_min=300, l_max=6000,
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
                 calculate_map_rmss_and_weight_stats=False,
                 rmss_and_wgts_from_coadds_or_individuals="",
                 rmss_and_wgts_from_signals_or_noises="",
                 calculate_noise_from_individual_maps=False,
                 calculate_noise_from_coadded_maps=False,
                 point_source_list_file=None,
                 calculate_cross_spectra_with_planck_map=False,
                 planck_map_fits_file=None,
                 calculate_pointing_discrepancies=False,
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
        
        
        # - Initialize variables related to collecting average numbers
        #   of various bolometer flagging statistics
        
        self.get_avgs_flagging_stats = collect_averages_from_flagging_stats
        
        if self.get_avgs_flagging_stats:
            self.bolo_props       = None
            self.pipe_info_frame  = None
            self.flagging_reasons = \
                ["BadCalSn", "BadWeight",  "Glitchy", "Oscillating",
                 "Latched",  "Overbiased", "BadHk",
                 "PostCalibrationNaNs",    "UnphysicalLowVariance",
                 "MissingFluxCalibration",
                 "Others",  "TotalNotFlagged", "TotalRemoved"]
            self.key_prefix_flg = "FlaggingStatisticsAverageNumbersOf"
            
            self.avgs_flagging_stats =  \
                {flagging_related_name: \
                 {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
                 for flagging_related_name in self.flagging_reasons}
        
        
        # - Initialize variables related to calculating the median value of
        #   the pW/K conversion factors for each wafer
        
        self.calc_median_pW_to_K = calculate_pW_to_K_conversion_factors
        
        if self.calc_median_pW_to_K:
            self.calframe   = None
            self.wafers     = ["W172", "W174", "W176", "W177", "W180",
                               "W181", "W188", "W203", "W204", "W206"]
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
        
        
        # - Initialize variables related to calculations of cross spectra
        #   with Planck maps
        
        self.calc_xspec_with_plck_map = calculate_cross_spectra_with_planck_map
        
        if self.calc_xspec_with_plck_map:
            self.planck_map_fits_file = planck_map_fits_file
            self.spt_to_planck_bands  = \
                {"90": "100", "150": "143", "220": "217"}
            self.mini_planck_map_frames = \
                {band: {sub_field: None for sub_field in self.map_sources} \
                 for band in self.spt_to_planck_bands.keys()}
            
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
        
        
        # - Initialize variables needed for most of the four analyses,
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
        
        if not numpy.isfinite(numpy.nanmean(numpy.asarray(map_frame["T"]))):
            self.log("")
            self.log("* There seem to be only NaNs in the map!")
            self.log("")
            return False
        
        return True
    
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Calibration:
            try:
                self.bolo_props = frame["BolometerProperties"]
                self.calframe   = frame
            except KeyError:
                pass
            return []
        
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
        
        if frame.type == core.G3FrameType.PipelineInfo:
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
                    center_ra  = core.G3Units.deg *   0.0
                    center_dec = core.G3Units.deg * -56.0
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
                    center_ra  = core.G3Units.deg * 0.0
                    center_dec = core.G3Units.deg * float(sbfd[-6:])
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
                                "the center of the full field: %7.3e"
                                %(numpy.asarray(
                                  self.coadded_map_frames[id_for_coadds]["T"])\
                                  [center_y][center_x]/tu))
                    self.detail("$$$ w of original map @ "
                                "the center of the full field: %7.3e"
                                %(numpy.asarray(
                                  self.coadded_map_frames\
                                  [id_for_coadds]["Wunpol"].TT)\
                                  [center_y][center_x]*tu*tu))                    
                    self.log("* ... field and weight maps were added "
                             "to the empty cache.")
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
            gc.collect()
            
            
            # -- Collect averages related to detector flagging
            
            if self.get_avgs_flagging_stats:
                
                if subtract_maps_in_this_frame:
                    for flg_typ in self.avgs_flagging_stats.keys():
                        self.avgs_flagging_stats[flg_typ] = \
                            self.remove_partial_mapmapdouble(
                                self.avgs_flagging_stats[flg_typ],
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
                                flg_typ = key.replace(self.key_prefix_flg, "")
                                avgs_flg_stats_from_this_fr[flg_typ] = \
                                    {id_for_coadds: frame[key]}
                    else:
                        self.log("")
                        self.log("* Gathering average numbers related to")
                        self.log("* detector flagging ...")
                        avgs_flg_stats_from_this_fr = {}
                        avgs_from_each_band = \
                            collect_averages_from_flagging_info(
                                self.pipe_info_frame,
                                self.bolo_props,
                                self.flagging_reasons)
                        for band, flg_typs_avgs in avgs_from_each_band.items():
                            if (band+"GHz") in id_for_coadds:
                                for flg_typ, avg in flg_typs_avgs.items():
                                    mmd = self.create_mmd_for_one_value(
                                              sbfd, oid, avg)
                                    avgs_flg_stats_from_this_fr[flg_typ] = \
                                        {id_for_coadds: mmd}
                                break
                    for flg_typ in avgs_flg_stats_from_this_fr.keys():
                        self.avgs_flagging_stats[flg_typ] = \
                            self.combine_mapmapdoubles(
                                self.avgs_flagging_stats[flg_typ],
                                avgs_flg_stats_from_this_fr[flg_typ])
                    self.log("* Done.")
                    self.log("")
            
            
            # -- Calculate mean value of the pW/K conversion factors
            #    for each band of each wafer
            
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
                                self.wafers, band, self.calframe,
                                flagging_stats=self.pipe_info_frame,
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
            
            
            # -- Decide whether it's time to perform some analysis on maps
            
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
                            summ_map_frame_individ = None
                            diff_map_frame_individ = None
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
                                    record_weights=False,
                                    logfun=self.detail)
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
                                    logfun=self.detail)
                            summ_map_frame_individ_mini = None
                            diff_map_frame_individ_mini = None
                            summ_map_frame_coadded_mini = None
                            diff_map_frame_coadded_mini = \
                                create_new_map_frame_with_smaller_region(
                                    diff_map_frame_coadded,
                                    center_ra, center_dec,
                                    t_only=self.t_only)
                            # * There are conceivable situations where
                            #   individual maps are needed, but these frames
                            #   suffice for the time being.
                    else:
                        time_to_analyze_maps   = True
                        summ_map_frame_individ = \
                            core.G3Frame(core.G3FrameType.Map)
                        tu = core.G3Units.mK
                        if self.t_only:
                            summ_map_frame_individ["T"] = \
                                mapmaker.mapmakerutils.remove_weight_t(
                                    frame["T"], frame["Wunpol"])
                            summ_map_frame_individ["Wunpol"] = \
                                frame["Wunpol"]
                            self.detail("$$$ t of original map "
                                        "@ field center: %7.3e [mK]"
                                        %(numpy.asarray(
                                          summ_map_frame_individ["T"])\
                                          [center_y][center_x]/tu))
                        diff_map_frame_individ = None
                        summ_map_frame_coadded = None
                        diff_map_frame_coadded = None
                        summ_map_frame_individ_mini = \
                            create_new_map_frame_with_smaller_region(
                                summ_map_frame_individ,
                                center_ra, center_dec,
                                t_only=self.t_only)
                        self.detail("$$$ t of smaller  map "
                                    "@ field center: %7.3e [mK]"
                                    %(numpy.asarray(
                                      summ_map_frame_individ["T"])\
                                      [center_y][center_x]/tu))
                        diff_map_frame_individ_mini = \
                            summ_map_frame_individ_mini
                        summ_map_frame_coadded_mini = None
                        diff_map_frame_coadded_mini = None
                        # * There are conceivable situations where
                        #   coadded maps are needed, but these frames
                        #   suffice for the time being.
                        gc.collect()
                if self.t_only:
                    del frame["T"], frame["Wunpol"]
                    gc.collect()
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
                            self.pixels_to_use_for_flc_calc[sbfd] = \
                                identify_pixels_of_non_atypical_region(
                                    summ_map_frame_coadded, center_dec,
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
                                    t_only=self.t_only)
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
                        offsets, fluxes, snrs = \
                            calculate_pointing_discrepancies(
                                summ_map_frame_individ, sbfd,
                                map_stddev=map_stddev)
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
                                    self.point_source_list_file)
                        if self.mini_planck_map_frames[band][sbfd] is None:
                            self.log("* (Need to make a Planck mini map,")
                            planck_band = self.spt_to_planck_bands[band]
                            fits_file = \
                                self.planck_map_fits_file.replace(
                                    "BAND", planck_band+"GHz")
                            self.mini_planck_map_frames[band][sbfd] = \
                                create_mini_planck_map_frame(
                                    fits_file, summ_map_frame_individ,
                                    center_ra, center_dec,
                                    t_only=self.t_only)
                            self.log("*  which will be used repeatedly.")
                        avg_xspec_ratios = \
                            calculate_average_ratio_of_spt_planck_xspectra(
                                summ_map_frame_individ_mini,
                                self.mini_planck_map_frames[band][sbfd],
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
                                    self.point_source_list_file)
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
                                    self.point_source_list_file)
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
            
            
            # - Finally try to free up some memory
            
            del frame
            gc.collect()
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
                    
                    for flg_typ in self.avgs_flagging_stats.keys():
                        k_rec = self.key_prefix_flg + flg_typ
                        mp_fr[k_rec] = self.avgs_flagging_stats[flg_typ][map_id]
                        self.log("- %s", flg_typ)
                        self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_median_pW_to_K:
                    self.log("# Medians of pW/K conversion factors:")
                    self.log("")
                    
                    for waf in self.wafers:
                        self.log(" - Wafer: %s", waf)
                        self.log("")
                        
                        for fa in self.calfactors:
                            k_rec = self.key_prefix_pwk + waf + fa
                            mp_fr[k_rec] = self.medians_of_temp_cal_factors \
                                           [waf][fa][map_id]
                            self.log(" "*3 + " - %s [%s]",
                                     fa, self.cal_fat_ustr[fa])
                            du = self.cal_fat_unis[fa]
                            self.print_mapmapdouble(mp_fr[k_rec], du, 6)
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




def do_weights_look_bad(mean_wt, band, sub_field):
    
    wu = 1 / (core.G3Units.mK * core.G3Units.mK)
    
    thresholds = {"90" : {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":  65.0*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi":  80.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi":  80.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi":  95.0*wu}},
                  "150": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi": 120.0*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi": 140.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi": 145.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi": 165.0*wu}},
                  "220": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":   8.0*wu},
                          "ra0hdec-52.25": {"lo": 0.0*wu, "hi":  10.0*wu},
                          "ra0hdec-59.75": {"lo": 0.0*wu, "hi":  12.0*wu},
                          "ra0hdec-67.25": {"lo": 0.0*wu, "hi":  14.0*wu}}}
    
    if (mean_wt < thresholds[band][sub_field]["lo"]) or \
       (mean_wt > thresholds[band][sub_field]["hi"]):
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
                 logging_function=print):
        
        self.t_only = t_only
        
        self.sb_fld    = None
        self.center_dc = None
        
        self.obs_id  = None
        self.map_ids = map_ids
        self.p_list  = point_source_list_file
        self.x_list  = bad_map_list_file
        self.logfun  = logging_function
        
        self.pixels_to_use_for_flc_calc = \
            {sub_field: None for sub_field in \
             ["ra0hdec-44.75", "ra0hdec-52.25",
              "ra0hdec-59.75", "ra0hdec-67.25"]}
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Observation:
            self.sb_fld = frame["SourceName"]
            self.center_dec = \
                core.G3Units.deg * float(self.sb_fld.replace("ra0hdec", ""))
        
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
            
            new_frame = core.G3Frame(core.G3FrameType.Map)
            
            if self.t_only:
                t_map_no_weights = mapmaker.mapmakerutils.remove_weight_t(
                                       frame["T"], frame["Wunpol"])
                new_frame["T"] = t_map_no_weights
                new_frame["Wunpol"] = frame["Wunpol"]
            
            
            # - Assume the map is not bad
            
            map_is_bad = False
                        
            # -- One criterion for cutting a map is
            #    mean TT weight of a map
            
            if self.pixels_to_use_for_flc_calc[self.sb_fld] is None:
                self.pixels_to_use_for_flc_calc[self.sb_fld] = \
                    identify_pixels_of_non_atypical_region(
                        frame, self.center_dec, self.p_list)
            
            flc_dct = calculate_map_fluctuation_metrics(
                          new_frame, band, self.sb_fld,
                          pixs=self.pixels_to_use_for_flc_calc[self.sb_fld],
                          t_only=self.t_only,
                          logfun=self.logfun)
            
            mean_wt = flc_dct["MeansOfTTWeights"]
            bad_weights = do_weights_look_bad(
                              mean_wt, band, self.sb_fld)
            if bad_weights:
                record_bad_obs_id(self.x_list,
                                  frame["Id"], self.sb_fld, self.obs_id,
                                  "Anomalously large weights.")
            
            # -- Another criterion is presence of NaNs
            
            has_nans = False
            
            if self.t_only:
                map_types = ["T"]
            
            for map_type in map_types:
                if False in numpy.isfinite(numpy.asarray(new_frame[map_type])):
                    has_nans = True
                    record_bad_obs_id(self.x_list,
                                      frame["Id"], self.sb_fld, self.obs_id,
                                      "NaNs in "+map_type+" map.")
            
            
            # - Check the badness again
            
            map_is_bad = bad_weights or has_nans
            
            if map_is_bad:
                frame["IgnoreThisMap"] = True
                self.logfun("* ... the map looks bad!")
            else:
                self.logfun("* ... the map doesn't look obviously bad!")
            self.logfun("")




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
        calculate_map_rmss_and_weight_stats=False,
        rmss_and_wgts_from_coadds_or_individuals=None,
        rmss_and_wgts_from_signals_or_noises=None,
        calculate_noise_from_individual_maps=False,
        calculate_noise_from_coadded_maps=False,
        point_source_list_file=None,
        calculate_cross_spectra_with_planck_map=False,
        planck_map_fits_file=None,
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
            # * if the name cannot be converted to an integer.
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
    
    if (len(all_good_g3_files) <= 1):
        log("")
        log("No applicable input files (or there is only one input file), ")
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
    log("")
    log(bad_obs_ids)
    log("")
    
    
    # - Construct a pipeline that analyzes and adds maps
    
    pipeline = core.G3Pipeline()
    
    pipeline.Add(core.G3Reader,
                 filename=all_good_g3_files)
    
    pipeline.Add(drop_duplicate_analysis_results)
    
    if trick_pipeline_to_receive_left_right_maps:
        pipeline.Add(FlagBadMaps,
                     map_ids=ids_to_append_directions_to,
                     t_only=temperature_maps_only,
                     point_source_list_file=point_source_list_file,
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
                 calculate_map_rmss_and_weight_stats=\
                     calculate_map_rmss_and_weight_stats,
                 rmss_and_wgts_from_coadds_or_individuals=\
                     rmss_and_wgts_from_coadds_or_individuals,
                 rmss_and_wgts_from_signals_or_noises=\
                     rmss_and_wgts_from_signals_or_noises,
                 calculate_noise_from_individual_maps=\
                     calculate_noise_from_individual_maps,
                 calculate_noise_from_coadded_maps=\
                     calculate_noise_from_coadded_maps,
                 point_source_list_file=point_source_list_file,
                 calculate_cross_spectra_with_planck_map=\
                     calculate_cross_spectra_with_planck_map,
                 planck_map_fits_file=planck_map_fits_file,
                 calculate_pointing_discrepancies=\
                     calculate_pointing_discrepancies,
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
                             "the word 'BAND' in it so that the paths of the "
                             "actual files can be obtained by replacing "
                             "'BAND' with '100GHz', '143GHz', or '217GHz'.")
    
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
                             "of observations IDs from the file and exclude "
                             "them from the analyses and coadding; rather, the "
                             "script only records maps that were found to be "
                             "bad when checking them as part of the analyses.")
    
    parser.add_argument("-G", "--log_file",
                        type=str, action="store", default=None,
                        help="The file to which the logger will send messages.")
    
    arguments = parser.parse_args()
    
    run(**vars(arguments))


# ==============================================================================

