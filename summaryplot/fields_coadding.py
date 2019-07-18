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
        multiply_2nd_frame_by=1.0,
        remove_weights_afterward=False,
        divide_results_by=1.0,
        record_weights=True):
    
    print("%%% PARAMETERS:")
    print("%%%", locals().items())
    
    if t_only:
        
        # (* Handle the case where a new frame (frame two) is being added to
        #    an empty cache (frame one))
        
        if "T" not in frame_one.keys():
            
            print("%%% adding a frame to the empty cache")
            
            new_frame = core.G3Frame(core.G3FrameType.Map)
            new_frame["T"] = frame_two["T"]
            new_frame["Wunpol"] = frame_two["Wunpol"]
            return new_frame
        
        iy = 9000
        ix = 9000
        print("%%% use", ix, "and", iy, "for the indices")
        print("%%% t1 initial:", numpy.asarray(frame_one["T"])[iy][ix])
        print("%%% t2 initial:", numpy.asarray(frame_two["T"])[iy][ix])
        print("%%% w1 initial:", numpy.asarray(frame_one["Wunpol"].TT)[iy][ix])
        print("%%% w2 initial:", numpy.asarray(frame_two["Wunpol"].TT)[iy][ix])
        
        # - First, remove weights if necessary
        
        if remove_weights_beforehand:
            t_map_one = mapmaker.mapmakerutils.remove_weight_t(
                            frame_one["T"], frame_one["Wunpol"])
            t_map_two = mapmaker.mapmakerutils.remove_weight_t(
                            frame_two["T"], frame_two["Wunpol"])
        else:
            t_map_one = frame_one["T"]
            t_map_two = frame_two["T"]
        
        print("%%% t1 after w. rm.:", numpy.asarray(t_map_one)[iy][ix])
        print("%%% t2 after w. rm.:", numpy.asarray(t_map_two)[iy][ix])
        
        # - Then, add field maps (and weight maps if applicable) together
        
        t_map = t_map_one + t_map_two * multiply_2nd_frame_by
        del t_map_one, t_map_two
        gc.collect()
        w_map = frame_one["Wunpol"].TT + \
                frame_two["Wunpol"].TT * multiply_2nd_frame_by
        
        print("%%% t", numpy.asarray(t_map)[iy][ix])
        print("%%% w", numpy.asarray(w_map)[iy][ix])
        
        # - After that, remove weights if necessary
        
        if remove_weights_afterward:
            t_map = mapmaker.mapmakerutils.remove_weight_t(
                        t_map, w_map)
        
        print("%%% t after w. rm.:", numpy.asarray(t_map)[iy][ix])
        
        # - Furthur divide the map by some number if requested
        
        if divide_results_by != 1.0:
            t_map = t_map / divide_result_by
        
        print("%%% t after div.:", numpy.asarray(t_map)[iy][ix])
        
        # - Finally, record the results in a new frame
        
        new_frame = core.G3Frame(core.G3FrameType.Map)
        new_frame["T"] = t_map
        if record_weights:
            new_frame["Wunpol"] = core.G3SkyMapWeights()
            new_frame["Wunpol"].TT = w_map
        else:
            del w_map
        
        print("%%% t final:", numpy.asarray(new_frame["T"])[iy][ix])
        try:
            print("%%% w final:", numpy.asarray(new_frame["Wunpol"].TT)[iy][ix])
        except:
            print("%%% w final: N/A")
        print("%%%", new_frame, "%%%")
    
    else:
        raise Exception("The case where map frames contain temperature "
                        "as well as polarization maps still needs to be "
                        "handled!")
    
    return new_frame
    
    
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    




def get_band_average(bolo_names_from_each_scan, bolo_props_map):

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
                n_bolos_from_one_scan[band] = \
                    n_bolos_from_one_scan.get(band) + 1

        for band, n_bolos in n_bolos_from_one_scan.items():
            n_bolos_from_each_scan[band].append(n_bolos)

    return {band: numpy.mean(ns_bolos) \
            for band, ns_bolos in n_bolos_from_each_scan.items()}




def collect_averages_from_flagging_info(
        pipe_info_frame, bolo_props_map, recognized_reasons):    
    
    averages = {"90": {}, "150": {}, "220": {}}
    
    # - Collect average numbers of bolos not flagged
        
    bolos_not_flagged_each_scan = pipe_info_frame["SurvivingBolos"]
    max_n_scans = numpy.max([len(scans) for reason, scans \
                             in pipe_info_frame["DroppedBolos"].items()])
    n_success_scans = len(bolos_not_flagged_each_scan.keys())
    
    diff_n_scans = max_n_scans - n_success_scans
    for i in range(diff_n_scans):
        bolos_not_flagged_each_scan["Dummy"+str(i)] = core.G3VectorString()
    
    avg_n_okay_bolos_each_band = \
        get_band_average(bolos_not_flagged_each_scan, bolo_props_map)
    for band, average in avg_n_okay_bolos_each_band.items():
        averages[band]["TotalNotFlagged"] = average
    
    # - Collect average numbers of bolos flagged for each reason
    
    all_bolos_flagged_each_scan = pipe_info_frame["DroppedBolos"]
    for flagging_reason, bolos_flagged_each_scan \
    in  all_bolos_flagged_each_scan.items():    
        bl_flgd_ea_sc_reorganized = core.G3MapVectorString()
        for scan_number, list_of_bolos in enumerate(bolos_flagged_each_scan):
            bl_flgd_ea_sc_reorganized[str(scan_number)] = list_of_bolos    
        avg_bad_bolos_each_band = \
            get_band_average(bl_flgd_ea_sc_reorganized, bolo_props_map)
        for band, average in avg_bad_bolos_each_band.items():
            averages[band][flagging_reason] = average
    
    for band, reasons_and_numbers in averages.items():
        averages[band]["Others"] = 0.0
        unrecognized_reasons = []
        for reason, number in reasons_and_numbers.items():
            if reason not in recognized_reasons:
                averages[band]["Others"] = \
                    averages[band].get("Others") + number
                unrecognized_reasons.append(reason)
        for reason in unrecognized_reasons:
            averages[band].pop(reason)
        for reason in recognized_reasons:
            if reason not in reasons_and_numbers.keys():
                averages[band][reason] = 0.0
    
    # - Collect average numbers of bolos removed
    
    all_bolos_flagged_each_scan = pipe_info_frame["DroppedBolos"]
    n_scans = numpy.min([len(bl_flgd_ea_sc) for flg_rsn, bl_flgd_ea_sc \
                         in all_bolos_flagged_each_scan.items()])
    bolos_removed_each_scan = core.G3MapVectorString()
    for i in range(n_scans):
        bolos_removed_one_scan = set([])
        for flagging_reason, bolos_flagged_each_scan \
        in  all_bolos_flagged_each_scan.items():
            bolos_removed_one_scan = \
                bolos_removed_one_scan | set(bolos_flagged_each_scan[i])
        bolos_removed_each_scan[str(i)] = \
            core.G3VectorString(bolos_removed_one_scan)
    avg_n_removed_bolos_each_band = \
        get_band_average(bolos_removed_each_scan, bolo_props_map)
    for band, average in avg_n_removed_bolos_each_band.items():
        averages[band]["TotalRemoved"] = average
    
    return averages




def collect_averages_of_pw_to_K_factors(wafers, calframe):
    
    bpm = calframe["BolometerProperties"]
    flu = core.G3Units.K * core.G3Units.rad * core.G3Units.rad
    flx = {"RCW38": { "90": 4.0549662e-07*flu,
                     "150": 2.5601153e-07*flu,
                     "220": 2.8025804e-07*flu},
           "MAT5A": { "90": 2.5738063e-07*flu,
                     "150": 1.7319235e-07*flu,
                     "220": 2.1451640e-07*flu}}
    
    values_by_band_and_wafer = \
        {band: {wafer: [] for wafer in wafers} for band in ["90", "150", "220"]}
    
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
    
    for bolo in available_bolos:
        band  = bpm[bolo].band / core.G3Units.GHz
        if not numpy.isfinite(band):
            continue
        band  = str(int(band))
        wafer = bpm[bolo].wafer_id.upper() 
        if (band  not in values_by_band_and_wafer.keys()) or \
           (wafer not in wafers):
            continue
        if band in calframe[sky_transmission_key].keys():
            pw_k  = calframe["CalibratorResponse"][bolo]
            pw_k *= calframe[flux_calibration_key][bolo]
            pw_k *= calframe[integral_flux_key][bolo]
            pw_k *= calframe[sky_transmission_key][band]
            pw_k /= flx[source][band]
        else:
            pw_k = numpy.nan
        values_by_band_and_wafer[band][wafer].append(pw_k)
    
    for band in values_by_band_and_wafer.keys():
        for wafer in wafers:
            vals = values_by_band_and_wafer[band][wafer]
            try:
                values_by_band_and_wafer[band][wafer] = numpy.nanmean(vals)
            except:
                values_by_band_and_wafer[band][wafer] = numpy.nan
    
    return values_by_band_and_wafer




def calculate_pointing_discrepancies(
        fr_one, fr_two, sub_field,
        temp_only=True, map_stddev=None):

    u = core.G3Units

    if temp_only:
        if fr_two != None:
            total_weight = fr_one["Wunpol"].TT + fr_two["Wunpol"].TT
            sum_map      = fr_one["T"] + fr_two["T"]
            map_to_use   = sum_map / total_weight
        else:
            map_to_use = mapmaker.mapmakerutils.remove_weight_t(
                             fr_one["T"], fr_one["Wunpol"])
    
    some_brightest_sources = \
        {"ra0hdec-44.75": {"1": numpy.array([352.32358, -47.50531, 1408.089]),
                           "2": numpy.array([314.06833, -47.24664, 1374.576]),
                           "3": numpy.array([ 41.50037, -46.85467,  715.329])},
         "ra0hdec-52.25": {"1": numpy.array([ 32.69288, -51.01703, 3819.563]),
                           "2": numpy.array([ 23.27404, -52.00094, 1026.665]),
                           "3": numpy.array([359.47312, -53.18686,  864.243])},
         "ra0hdec-59.75": {"1": numpy.array([ 47.48363, -60.97761,  869.843]),
                           "2": numpy.array([ 45.96104, -62.19042,  832.717]),
                           "3": numpy.array([ 14.69433, -56.98650,  785.528])},
         "ra0hdec-67.25": {"1": numpy.array([329.27542, -69.68981, 1114.524]),
                           "2": numpy.array([337.25092, -69.17492,  445.331]),
                           "3": numpy.array([325.44375, -64.18742,  388.218])}}    
    discrep_dict = {"1": {"delta_ra": numpy.nan, "delta_dec": numpy.nan},
                    "2": {"delta_ra": numpy.nan, "delta_dec": numpy.nan},
                    "3": {"delta_ra": numpy.nan, "delta_dec": numpy.nan}}
    flux_dict = {"1": numpy.nan, "2": numpy.nan, "3": numpy.nan}
    snr_dict  = {"1": numpy.nan, "2": numpy.nan, "3": numpy.nan}
    
    info_on_point_sources = some_brightest_sources[sub_field]
    
    for point_source_number, info in info_on_point_sources.items():
        
        true_right_ascension = info[0]
        if true_right_ascension > 180:
            true_right_ascension -= 360
        true_right_ascension *= u.deg
        
        true_declination  = info[1]
        true_declination *= u.deg
        
        true_x, true_y = \
            map_to_use.angle_to_xy(true_right_ascension,
                                   true_declination)
        nominal_pixel_number = \
            map_to_use.angle_to_pixel(true_right_ascension,
                                      true_declination)
        nominal_y_index, nominal_x_index = \
            numpy.unravel_index(
                nominal_pixel_number, map_to_use.shape)
        
        nominal_x_index = int(nominal_x_index)
        nominal_y_index = int(nominal_y_index)
        
        mini_map_width    = 4 * core.G3Units.arcmin
        mini_map_height   = 4 * core.G3Units.arcmin
        mini_map_x_length = int(mini_map_width  / map_to_use.x_res)
        mini_map_y_length = int(mini_map_height / map_to_use.y_res)
        half_x_length     = mini_map_x_length // 2
        half_y_length     = mini_map_y_length // 2
        mini_map_x_left   = nominal_x_index - half_x_length
        mini_map_x_right  = nominal_x_index + half_x_length
        mini_map_y_bottom = nominal_y_index - half_y_length
        mini_map_y_top    = nominal_y_index + half_y_length
        mini_map = numpy.asarray(map_to_use)\
                   [mini_map_y_bottom:mini_map_y_top+1,
                    mini_map_x_left:mini_map_x_right+1]
        
        try:
            fit_parameters = util.fitting.fit_gaussian2d(mini_map)
        except ValueError:
            discrep_dict[point_source_number]["delta_ra"]  = numpy.nan
            discrep_dict[point_source_number]["delta_dec"] = numpy.nan
            continue
            
        gaussian_center_x = fit_parameters[1]
        gaussian_center_y = fit_parameters[2]
        
        measured_x = gaussian_center_x - half_x_length + nominal_x_index
        measured_y = gaussian_center_y - half_y_length + nominal_y_index
        
        measured_right_ascension, measured_declination = \
            map_to_use.xy_to_angle(measured_x, measured_y)
        
        delta_ra  = measured_right_ascension - true_right_ascension
        delta_dec = measured_declination     - true_declination
        effective_delta_ra = delta_ra * numpy.cos(true_declination/u.rad)
        
        discrep_dict[point_source_number]["delta_ra"]  = effective_delta_ra
        discrep_dict[point_source_number]["delta_dec"] = delta_dec
        
        nearest_x = int(numpy.round(measured_x))
        nearest_y = int(numpy.round(measured_y))
        source_centered_mini_map =   \
            numpy.asarray(map_to_use)\
                [nearest_y-half_y_length:nearest_y+half_y_length+1,
                 nearest_x-half_x_length:nearest_x+half_x_length+1]
        source_flux = numpy.sum(source_centered_mini_map) * \
                      map_to_use.x_res * map_to_use.y_res
        flux_dict[point_source_number] = source_flux
        
        if map_stddev is None:
            source_snr = numpy.nan
        else:
            source_snr = numpy.max(source_centered_mini_map) / map_stddev
        snr_dict[point_source_number] = source_snr
        
    return discrep_dict, flux_dict, snr_dict




def calculate_map_fluctuation_metrics(
        map_frame, temperature_only, return_nans=False):
    
    keys = ["MapStdDevs",
            "MeansOfWeights",
            "NumbersOfPixelsWithGoodWeights",
            "VarsOfProductsOfMapValuesAndSqrtWeights"]
    
    fluctuation_metrics = {key: numpy.nan for key in keys}
    
    if return_nans:
        return fluctuation_metrics
    else:
        if temperature_only:
            if map_frame["T"].is_weighted:
                t_vals = numpy.asarray(mapmaker.mapmakerutils.remove_weight_t(
                                       map_frame["T"], map_frame["Wunpol"]))
            else:
                t_vals = numpy.asarray(map_frame["T"])
            w_vals = numpy.asarray(map_frame["Wunpol"].TT)
        
        w_cut = numpy.percentile(w_vals[numpy.where((w_vals>0.0) & numpy.isfinite(w_vals))], 25)
        igw   = numpy.where((w_vals > w_cut) & numpy.isfinite(w_vals))
        npnzw = len(numpy.where(w_vals > 0.0)[0])
        if len(igw[0]) == 0:
            return fluctuation_metrics
        
        good_w_vals = w_vals[igw]
        good_t_vals = t_vals[igw]
        
        pctl_hi = numpy.nanpercentile(good_t_vals, 99.5)
        pctl_lo = numpy.nanpercentile(good_t_vals,  0.5)
        inx = numpy.where(numpy.isfinite(good_t_vals) & (good_t_vals < pctl_hi) & (good_t_vals > pctl_lo))
        better_t_vals = good_t_vals[inx]
        better_w_vals = good_w_vals[inx]
        
        map_stddev  = numpy.std(better_t_vals)
        mean_weight = numpy.mean(better_w_vals)
        var_product = numpy.var(better_t_vals * numpy.sqrt(better_w_vals))
        
        fluctuation_metrics[keys[0]] = map_stddev
        fluctuation_metrics[keys[1]] = mean_weight
        fluctuation_metrics[keys[2]] = npnzw
        fluctuation_metrics[keys[3]] = var_product
        
        return fluctuation_metrics




def create_new_map_frame_with_smaller_map_region(
        map_frame, temperature_only,
        center_ra, center_dec,
        map_key_suffix="", divisive_factor=1.0):

    new_mp_fr = core.G3Frame(core.G3FrameType.Map)

    if temperature_only:
        original_sky_map = map_frame["T"+map_key_suffix] / divisive_factor
        
        smaller_sky_map = \
            coordinateutils.FlatSkyMap(
                x_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                y_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                res=original_sky_map.res,
                proj=original_sky_map.proj,
                alpha_center=center_ra,
                delta_center=center_dec,
                pol_type=original_sky_map.pol_type,
                coord_ref=original_sky_map.coord_ref)
        coordinateutils.reproj_map(original_sky_map,  smaller_sky_map,  1)
        new_mp_fr["T"] = smaller_sky_map
        
        if "Wunpol" in map_frame.keys():
            original_wght_map = map_frame["Wunpol"].TT
            
            smaller_wght_map = coordinateutils.FlatSkyMap(
                x_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                y_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                res=original_sky_map.res,
                proj=original_sky_map.proj,
                alpha_center=center_ra,
                delta_center=center_dec,
                pol_type=original_sky_map.pol_type,
                coord_ref=original_sky_map.coord_ref)
            coordinateutils.reproj_map(original_wght_map, smaller_wght_map, 1)    
            new_mp_fr["Wunpol"]    = core.G3SkyMapWeights()
            new_mp_fr["Wunpol"].TT = smaller_wght_map

    return new_mp_fr




def create_apodization_mask(map_frame, point_source_file):

    new_mp_fr = core.G3Frame(core.G3FrameType.Map)

    ptsrc_msk = mapspectra.apodmask.makeApodizedPointSourceMask(
                    map_frame, point_source_file)

    t_map     = map_frame["T"]
    brdr_msk  = coordinateutils.FlatSkyMap(
                    x_len=t_map.shape[1], y_len=t_map.shape[0],
                    res=t_map.res, proj=t_map.proj,
                    alpha_center=t_map.alpha_center,
                    delta_center=t_map.delta_center,
                    pol_type=t_map.pol_type,
                    coord_ref=t_map.coord_ref, units=None)

    twod_window  = numpy.ones(brdr_msk.shape)
    twod_window *= signal.tukey(brdr_msk.shape[1], alpha=0.2)
    twod_window  = (twod_window.transpose() *\
                    signal.tukey(brdr_msk.shape[0], alpha=0.2)).\
                   transpose()
    numpy.asarray(brdr_msk)[:,:] = twod_window

    full_mask = ptsrc_msk * brdr_msk

    return full_mask




def decide_operation_to_do_with_new_map(
        frame, past_operations, past_oids,
        temp_only, center_ra, center_dec):
    
    mp_fr = create_new_map_frame_with_smaller_map_region(
                frame, temp_only, center_ra, center_dec)
    
    if (False in numpy.isfinite(numpy.asarray(mp_fr["T"]))) or \
       (0.0   in numpy.asarray(mp_fr["Wunpol"].TT)):
        return "Ignore"
        
    past_operations = [past_operations[str(oid)] for oid in past_oids[:-1]]
    
    if (len(past_operations) == 0) or \
       (len(past_operations) == past_operations.count(0.0)):
        return "Copy"
    else:
        n_add = len([o for o in past_operations if o ==  1.0])
        n_sub = len([o for o in past_operations if o == -1.0])
        if n_add > n_sub:
            return "Subtract"
        else:
            return "Add"




def calculate_noise_level(
        fr_one, fr_two, ptsrc_lst, temp_only=False,
        smaller_region=False, center_ra=None, center_dec=None,
        map_key_suffix="", divisive_factor=1.0):

    if fr_two != None:
        mp_fr = map_analysis.subtract_two_maps(
                    fr_one, fr_two, divide_by_two=True, in_place=False)
    else:
        mp_fr = fr_one
    
    if smaller_region:
        mp_fr = create_new_map_frame_with_smaller_map_region(
                    mp_fr, temp_only,
                    center_ra, center_dec,
                    map_key_suffix=map_key_suffix,
                    divisive_factor=divisive_factor)
    
    mask = create_apodization_mask(mp_fr, ptsrc_lst)
    
    cls = map_analysis.calculateCls(
              mp_fr, cross_map=None, t_only=temp_only,
              apod_mask=mask, kspace_filt=None, tf_2d=None,
              ell_bins=None, ell_min=300, ell_max=6000, delta_ell=50,
              return_2d=False, realimag="real", in_place=False)
    
    idx = numpy.where((cls["ell"]>3000) & (cls["ell"]<5000))[0]
        
    if temp_only:
        noise = numpy.sqrt(numpy.mean(cls["TT"][idx]))
    
    return noise




def calculate_average_cl_of_cross_spectrum(
        frame_one, frame_two, ptsrc_lst, temp_only=False,
        smaller_region=False, center_ra=None, center_dec=None):
    
    if smaller_region:
        frame_one = create_new_map_frame_with_smaller_map_region(
                        frame_one, temp_only,
                        center_ra, center_dec)
        frame_two = create_new_map_frame_with_smaller_map_region(
                        frame_two, temp_only,
                        center_ra, center_dec)
    
    mask = create_apodization_mask(frame_one, ptsrc_lst)

    cls = map_analysis.calculateCls(
              frame_one, cross_map=frame_two, t_only=temp_only,
              apod_mask=mask, kspace_filt=None, tf_2d=None,
              ell_bins=None, ell_min=300, ell_max=6000, delta_ell=50,
              return_2d=False, realimag="real", in_place=False)
    
    idx = numpy.where((cls["ell"]>1000) & (cls["ell"]<2500))[0]
        
    if temp_only:
        rtcl = numpy.sqrt(numpy.mean(cls["TT"][idx]))
    
    return rtcl


# ==============================================================================




# ==============================================================================
# Define the pipeline module that calls the functions defined above to
# add maps together and perform the anlyses
# ------------------------------------------------------------------------------


class AnalyzeAndCoaddMaps(object):
    
    def __init__(self, map_ids, map_sources, temperature_only,
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
                 planck_map_fits_files=None,
                 calculate_pointing_discrepancies=False,
                 logging_function=logging.info):
        
        self.log = logging_function
        
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
        self.t_only      = temperature_only
        
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
            self.pipe_line_info   = None
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
        #   the pW/K conversion factors for each band of each wafer
        
        self.calc_mean_pW_to_K = calculate_pW_to_K_conversion_factors
        
        if self.calc_mean_pW_to_K:
            self.calframe = None
            self.wafers   = ["W172", "W174", "W176", "W177", "W180",
                             "W181", "W188", "W203", "W204", "W206"]
            self.key_prefix_pwk = "MeansOfTemperatureCalibrationFactors"
            
            self.means_of_temp_cal_factors = \
                {wafer: {map_id: core.G3MapMapDouble() \
                         for map_id in self.map_ids}   \
                 for wafer in self.wafers}
        
        
        # - Initialize variables related to calculating a Std. dev. of
        #   a map, a mean TT weight, and number of pixels in the weight map
        #   whose weights are above some thresholds.
        
        self.calc_map_stddev_etc = calculate_map_rmss_and_weight_stats
        
        if self.calc_map_stddev_etc:
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
                self.signals_of_noises_for_flc.append("Noise")
            
            self.map_types_for_flc = []
            for c_or_i in self.coadds_or_individs_for_flc:
                for s_or_n in self.signals_of_noises_for_flc:
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
            self.planck_map_fits_files = planck_map_fits_files
            self.spt_to_planck_bands  = \
                {"90GHz": "100GHz", "150GHz": "143GHz", "220GHz": "217GHz"}
            self.mini_planck_map_frames = \
                {band: {sub_field: None for sub_field in self.map_sources} \
                 for band in self.spt_to_planck_bands.keys()}
            
            self.av_xspec_ratio_key = "AveragesOfRatiosOfSPTxPlancktoPlckxPlck"
            
            self.avgs_xspec_ratios = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        
        
        # - Initialize variables related to noise calculations
        
        self.calc_noise_from_individ_maps = calculate_noise_from_individual_maps
        self.calc_noise_from_coadded_maps = calculate_noise_from_coadded_maps
        
        if self.calc_noise_from_individ_maps:
            self.indvid_noise_key = "NoiseFromIndividualMaps"
            self.noise_from_individual_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        if self.calc_noise_from_coadded_maps:
            self.coadds_noise_key = "NoiseFromCoaddedMaps"
            self.noise_from_coadded_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        
        
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
            self.log(" "*(indlev+3) + "%s", oids)
            self.log(" "*(indlev+3) + "%4.1e", [data[oid]/u for oid in oids])
            self.log("")
    
    
    def map_seems_fine(self, map_frame):
        
        if "IgnoreThisMap" in map_frame.keys():
            self.log("")
            self.log("* This map was determined to be bad from ")
            self.log("* a previous processing step!")
            return False
        
        if map_frame["T"].units != core.G3TimestreamUnits.Tcmb:
            self.log("")
            self.log("* The units of this map are not in Tcmb!")
            return False
        
        if not numpy.isfinite(numpy.nanmean(numpy.asarray(map_frame["T"]))):
            self.log("")
            self.log("* There seem to be only NaNs in the map!")
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
                self.pipe_info = frame
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
                    self.log("\n")
                    return []
                else:
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
                    self.log("\n")
                    return []
                elif self.obs_info is None:
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because it is unclear what type of observation")
                    self.log("* the map is from.")
                    self.log("\n")
                    return []
                elif self.obs_info["SourceName"] not in self.map_sources:
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because of the wrong sub-field.")
                    self.log("\n")
                    return []
                else:
                    id_for_coadds = self.id_mapping[frame["Id"]]
                    oid  = self.obs_info["ObservationID"]
                    sbfd = self.obs_info["SourceName"]
                    dmid = str(oid) + "_" + frame["Id"]   # * a more detailed map id
                    t_i  = self.obs_info["ObservationStart"]
                    t_f  = self.obs_info["ObservationStop"]
                    d_i  = std_processing.time_to_obsid(t_i)
                    d_f  = std_processing.time_to_obsid(t_f)
                    dura = {str(oid): (d_f - d_i) * core.G3Units.s}
                    center_ra  = core.G3Units.deg * 0.0
                    center_dec = core.G3Units.deg * float(sbfd[-6:])
                    for b in ["90GHz", "150GHz", "220GHz"]:
                        if b in id_for_coadds:
                            band = b
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
                self.log("* so, it will be skipped ...")
                self.log("\n")
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
                self.log("* so, it will be skipped ...")
                self.log("\n")
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
                        multiply_2nd_frame_by=-1.0,
                        remove_weights_afterward=False,
                        divide_results_by=1.0,
                        record_weights=True)
                self.log("* Done.")
                self.log("")
            else:
                self.log("")
                self.log("* Adding a set of ")
                self.log("* sky map and weight map (ID: %s)", frame["Id"])
                self.log("* to the coadded maps (ID: %s) ...", id_for_coadds)
                self.coadded_map_frames[id_for_coadds] = \
                    add_two_map_frames(
                        self.coadded_map_frames[id_for_coadds],
                        frame,
                        t_only=self.t_only,
                        remove_weights_beforehand=False,
                        multiply_2nd_frame_by=1.0,
                        remove_weights_afterward=False,
                        divide_results_by=1.0,
                        record_weights=True)
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
                                flg_typ = key.replace(prefix, "")
                                avgs_flg_stats_from_this_fr[flg_typ] = \
                                    {id_for_coadds: frame[key]}
                    else:
                        self.log("")
                        self.log("* Gathering average numbers related to")
                        self.log("* detector flagging  ...")
                        avgs_flg_stats_from_this_fr = {}
                        avgs_from_each_band = \
                            collect_averages_from_flagging_info(
                                self.pipe_info,
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
                    self.log("* Done.")
                    self.log("")
                    for flg_typ in avgs_flg_stats_from_this_fr.keys():
                        self.avgs_flagging_stats[flg_typ] = \
                            combine_mapmapdoubles(
                                self.avgs_flagging_stats[flg_typ],
                                avgs_flg_stats_from_this_fr[flg_typ])
            
            
            # -- Calculate mean value of the pW/K conversion factors
            #    for each band of each wafer
            
            if self.calc_mean_pW_to_K:
                
                if subtract_maps_in_this_frame:
                    for wafer in self.means_of_temp_cal_factors.keys():
                        self.means_of_temp_cal_factors[wafer] = \
                            self.remove_partial_mapmapdouble(
                                    self.means_of_temp_cal_factors[wafer],
                                    obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering the mean values of the pW/K")
                        self.log("* coversion factors that were ")
                        self.log("* calculated previously ...")
                        means_of_convs_from_this_fr = {}
                        for key in frame.keys():
                            if self.key_prefix_pwk in key:
                                wafer = key.replace(prefix, "")
                                means_of_convs_from_this_fr[wafer] = \
                                    {id_for_coadds: frame[key]}
                    else:
                        self.log("")
                        self.log("* Gathering the mean values of the pW/K")
                        self.log("* conversion factors ...")
                        means_of_convs_from_this_fr = {}
                        avgs_to_be_reorganized = \
                            collect_averages_of_pw_to_K_factors(
                                self.wafers,
                                self.calframe)
                        for band, wfs_avgs in avgs_to_be_reorganized.items():
                            if (band+"GHz") in id_for_coadds:
                                for wafer, avg in wfs_avgs.items():
                                    mmd = self.create_mmd_for_one_value(
                                              sbfd, oid, avg)
                                    means_of_convs_from_this_fr[wafer] = \
                                        {id_for_coadds: mmd}
                                break
                    self.log("* Done.")
                    self.log("")
                    for wafer in means_of_convs_from_this_fr.keys():
                        self.means_of_temp_cal_factors[wafer] = \
                            combine_mapmapdoubles(
                                self.means_of_temp_cal_factors[wafer],
                                means_of_convs_from_this_fr[wafer])
            
            
            # -- Decide whether it's time to perform some analysis on maps
            
            if self.analyze_maps:
                
                if frame_has_old_coadds:
                    time_to_analyze_maps = True
                    # * Set this to True just to collect 
                    #   the previous analysis results
                    if self.maps_split_by_direction:
                        n_mp_ea_fld = {fld: len(obs_ids) for fld, obs_ids in \
                                       frame["CoaddedObservationIDs"].items()}
                        for f, n in n_mp_ea_fld.items():
                            self.map_fr_arrival_counters[frame["Id"]][f] += n
                
                else:
                    if self.maps_split_by_direction:
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
                            n_mp_r = self.map_cache_for_both_dirs[idr][sbfd]
                            if (n_mp_l > 0) and (n_mp_l == n_mp_r):
                                time_to_analyze_maps = True
                                break
                        
                        if time_to_analyze_maps:
                            self.log("")
                            self.log("* Since both the %s and %s maps of "
                                     "observation %s have been cached,",
                                     id_l_ver, id_r_ver, oid)
                            self.log("* the requested analyses can now start.")
                            self.log("")
                        else:
                            self.log("")
                            self.log("* Since either %s or %s map of "
                                     "observation %s has not been cached yet,")
                            self.log("* the requested analyses need to wait.")
                            self.log("")
                        
                        if time_to_analyze_maps:
                            summ_map_frame_individ = None
                            diff_map_frame_individ = None
                            summ_map_frame_coadded = \
                                add_two_map_frames(
                                    self.coadded_map_frames[id_l_ver],
                                    self.coadded_map_frames[id_r_ver],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=False,
                                    multiply_2nd_frame_by=1.0,
                                    remove_weights_afterward=True,
                                    divide_results_by=1.0,
                                    record_weights=False)
                            diff_map_frame_coadded = \
                                add_two_map_frames(
                                    self.coadded_map_frames[id_l_ver],
                                    self.coadded_map_frames[id_r_ver],
                                    t_only=self.t_only,
                                    remove_weights_beforehand=True,
                                    multiply_2nd_frame_by=-1.0,
                                    remove_weights_afterward=False,
                                    divide_results_by=2.0,
                                    record_weights=False)
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
                            gc.collect()
                    else:
                        time_to_analyze_maps   = True
                        summ_map_frame_individ = \
                            add_two_map_frames(
                                frame,
                                None,
                                t_only=self.t_only,
                                remove_weights_beforehand=False,
                                multiply_2nd_frame_by=0.0,
                                remove_weights_afterward=True,
                                divide_result_by=1.0,
                                record_weights=False)
                        diff_map_frame_individ = None
                        summ_map_frame_coadded = None
                        diff_map_frame_coadded = None
                        summ_map_frame_individ_mini = \
                            create_new_map_frame_with_smaller_region(
                                summ_map_frame_individ,
                                center_ra, center_dec,
                                t_only=self.t_only)
                        diff_map_frame_individ_mini = \
                            summ_map_frame_indvid_mini
                        summ_map_frame_coadded_mini = None
                        diff_map_frame_coadded_mini = None
                        # * There are conceivable situations where
                        #   coadded maps are needed, but these frames
                        #   suffice for the time being.
                        gc.collect()
            
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
                        if self.pixels_to_use_for_flc_calc[src] is None:
                            self.pixels_to_use_for_flc_calc[src] = \
                                identify_pixels_of_non_atypical_region(
                                    frame, center_ra, center_dec,
                                    self.point_source_list_file)
                        for mt in self.map_types_for_flc:
                            if map_type == "IndividualSignalMaps":
                                frame_to_use = summ_map_frame_individ
                            if map_type == "CoaddedSignalMaps":
                                frame_to_use = summ_map_frame_coadded
                            if map_type == "IndividualNoiseMaps":
                                frame_to_use = diff_map_frame_indvid
                            if map_type == "CoaddedNoiseMaps":
                                frame_to_use = diff_map_frame_coadded
                            fluctuation_metrics = \
                                calculate_map_fluctuation_metrics(
                                    frame_to_use,
                                    pixs=self.pixels_to_use_for_flc_calc[src],
                                    t_only=self.t_only)
                            self.log("* %s", mt)
                            tu = core.G3Units.mK
                            for k, v in fluctuation_metrics.items():
                                self.log("*  %35s %s %s",
                                         k, v/self.flu_met_unis[k],
                                         self.flu_met_ustr[k])
                            for k, v in fluctuation_metrics.items():
                                mmd = self.create_mmd_for_one_value(
                                          sbfd, oid, v)
                                self.map_fluctuation_metrics[mt][k] = \
                                    self.combine_mapmapdoubles(
                                        self.map_fluctuation_metrics[mt][fk],
                                        {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate pointing discrepancies for three bright sources
            #     in the sub-field
            
            if time_to_analyze_maps and self.calc_pointing_discrepancies:
                
                if subtract_maps_in_this_frame:
                    for rank in self.ptsrc_ranks:
                        for delta_type in self.ptsrc_offset_types:
                            self.ptsrc_delta_ras_and_decs[rank][delta_type] = \
                                self.remove_partial_mapmapdouble(
                                    self.delta_ras_and_decs[rank][delta_type],
                                    obs_ids_from_this_frame)
                        self.ptsrc_fluxes[rank] = \
                            self.remove_partial_mapmapdouble(
                                self.ptsrc_fluxes[rank],
                                obs_ids_from_this_frame)
                        self.snrs[n] = \
                            self.remove_partial_mapmapdouble(
                                self.ptsrc_snrs[rank],
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
                            map_stddev = self.map_fluctuation_metrics  \
                                             ["IndividualSignalMaps"]  \
                                             ["MapStdDevs"]            \
                                             [id_for_coadds][src][str(oid)]
                        else:
                            map_stddev = None
                        offsets, fluxes, snrs = \
                            calculate_pointing_discrepancies(
                                summ_map_frame_individ, sbfd,
                                map_stddev=map_stddev)
                        for rank, diffs in offsets.items():
                            u = core.G3Units
                            for coord, diff in diffs.items():
                                self.log("* Source %s %9s %s arcseconds",
                                         rank, coord, diff/u.arcsec)
                            flux  = fluxes[rank]
                            flux /= u.mK * u.arcmin * u.arcmin
                            self.log("* Source %s %9s %s mK.arcmin^2",
                                     "flux", rank, flux)
                            snr = snrs[rank]
                            self.log("* Source %s %9s %s",
                                     "SNR", rank, snr)
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
                    self.avgs_xspec_ratios = \
                        self.remove_partial_mapmapdouble(
                            self.avgs_xspec_ratios,
                            obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering averages of ratios of")
                        self.log("* cross spectra (SPT x Planck / Pla x Pla)")
                        self.log("* that were calculated previously ...")
                        self.avgs_xspec_ratios = \
                            self.combine_mapmapdoubles(
                                self.avgs_xspec_ratios,
                                {id_for_coadds: frame[self.av_xspec_ratio_key]})
                    else:
                        self.log("")
                        self.log("* Calculating the average of the ratio of")
                        self.log("* SPT x Planck / Planck x Planck in a")
                        self.log("* low ell range ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_for_small_region(
                                    summ_map_frame_individ_mini,
                                    self.point_source_list_file)
                        if self.mini_planck_maps[band][sbfd] is None:
                            self.mini_planck_maps[band][sbfd] = \
                                create_a_mini_planck_map_frame(
                                    fits_file, frame, center_ra, center_dec,
                                    t_only=self.t_only)
                        avg_xspec_ratio = \
                            calculate_average_ratio_of_spt_planck_xspectra(
                                summ_map_frame_individ_mini,
                                self.mini_planck_map_frames[band][sbfd],
                                self.masks_for_powspec_calculations[sbfd],
                                t_only=self.t_only)
                        self.log("* ... the average was calculated to be %s.",
                                 avg_xspec_ratio)
                        mmd = self.create_mmd_for_one_value(
                                  sbfd, oid, avg_xspec_ratio)
                        self.avgs_xspec_ratios = \
                            self.combine_mapmapdoubles(
                                self.avgs_xspec_ratios,
                                {id_for_coadds: mmd})
                    self.log("* Done.")
                    self.log("")
            
            
            # --- Calculate noise levels from individual maps
            
            if time_to_analyze_maps and self.calc_noise_from_individ_maps:
                
                if subtract_maps_in_this_frame:
                    self.noise_from_individual_maps = \
                        self.remove_partial_mapmapdouble(
                            self.noise_from_individual_maps,
                            obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering noise levels of individual maps ")
                        self.log("* that were calculated previously ...")
                        self.noise_from_individual_maps = \
                            self.combine_mapmapdoubles(
                                self.noise_from_individual_maps,
                                {id_for_coadds: frame[self.indvid_noise_key]})
                    else:
                        self.log("")
                        self.log("* Calculating noise level "
                                 "from this observation ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_for_small_region(
                                    diff_map_frame_indvid_mini,
                                    self.point_source_list_file)
                        noise = calculate_noise_level(
                                    diff_map_frame_individ_mini,
                                    self.masks_for_powspec_calculations[sbfd],
                                    t_only=self.t_only)
                        nu = core.G3Units.uK * core.G3Units.arcmin
                        self.log("* ... the noise level was calculated to be ")
                        self.log("* %s uK.arcmin.", noise/nu)
                        mmd = self.create_mmd_for_one_value(
                                  sbfd, oid, noise)
                        self.noise_from_individual_maps = \
                            self.combine_mapmapdoubles(
                                self.noise_from_individual_maps,
                                {id_for_coadds: mmd})
            
            
            # --- Calculate noise levels from coadded maps
            
            if time_to_analyze_maps and self.calc_noise_from_coadded_maps:
                
                if subtract_maps_in_this_frame:
                    self.noise_from_coadded_maps = \
                        self.remove_partial_mapmapdouble(
                            self.noise_from_coadded_maps,
                            obs_ids_from_this_frame)
                
                else:
                    if frame_has_old_coadds:
                        self.log("")
                        self.log("* Gathering noise levels of coadded maps ")
                        self.log("* that were calculated previously ...")
                        self.noise_from_coadded_maps = \
                            self.combine_mapmapdoubles(
                                self.noise_from_coadded_maps,
                                {id_for_coadds: frame[self.coadds_noise_key]})
                    else:
                        self.log("")
                        self.log("* Calculating noise level "
                                 "from the running coadded maps ...")
                        if self.masks_for_powspec_calculations[sbfd] is None:
                            self.masks_for_powspec_calculations[sbfd] = \
                                create_mask_for_powspec_calc_for_small_region(
                                    diff_map_frame_coadded_mini,
                                    self.point_source_list_file)
                        noise = calculate_noise_level(
                                    diff_map_frame_coadded_mini,
                                    self.masks_for_powspec_calculations[sbfd],
                                    t_only=self.t_only)
                        nu = core.G3Units.uK * core.G3Units.arcmin
                        self.log("* ... the noise level was calculated to be ")
                        self.log("* %s uK.arcmin.", noise/nu)
                        mmd = self.create_mmd_for_one_value(
                                  sbfd, oid, noise)
                        self.noise_from_individual_maps = \
                            self.combine_mapmapdoubles(
                                self.noise_from_individual_maps,
                                {id_for_coadds: mmd})
            
            
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
                
                if self.maps_split_by_direction:
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
                
                self.log("# Observation durations [minutes.]:")
                self.log("")
                du = core.G3Units.min   # * display units
                for sf, data in mp_fr["ObservationDurations"].items():
                    self.log(" "*3 + "- %s", sf)
                    self.log("")
                    oids = data.keys()
                    self.log(" "*6 + "%s", oids)
                    self.log(" "*6 + "%s", [data[oid]/du for oid in oids])
                    self.log("")
                self.log("\n")
                
                
                if self.get_avgs_flagging_stats:
                    self.log("# Averages related to flagging statistics:")
                    self.log("")
                    
                    for flg_typ in self.avgs_flagging_stats.keys():
                        k_rec = self.key_prefix_flg + flg_typ
                        mp_fr[k_rec] = self.avgs_flagging_stats[flg_typ][map_id]
                        self.log("- %s", flg_typ)
                        self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_mean_pW_to_K:
                    self.log("# Averages of pW/K conversion factors:")
                    self.log("")
                    
                    du = core.G3Units.pW / core.G3Units.K
                    for w in self.means_of_temp_cal_factors.keys():
                        k_rec = self.key_prefix_pwk + wafer 
                        mp_fr[k_rec] = self.means_of_temp_cal_factors[w][map_id]
                        self.log("- %s", w)
                        self.print_mapmapdouble(mp_fr[k_rec], du, 3)
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
                                self.map_fluctuation_metrics[mp][fk][map_id]
                            self.log(" "*3 + " - Metric: %s [%s]",
                                     fk, self.flu_met_unis[fk])
                            self.print_mapmapdouble(
                                mp_fr[k_rec], self.flu_met_unis[fk], 6)
                    self.log("\n")
                
                
                if self.calc_pointing_discrepancies:
                    self.log("# Pointing offsets, fluxes, and SNRs of sources:")
                    self.log("")
                    
                    for r in self.ptsrc_ranks:
                        self.log("- Point source %s", rank)
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
                        mp_fr[k_rec] = self.ptsrc_snrs[rank][map_id]                      
                        self.log(" "*3 + " - SNR")
                        self.print_mapmapdouble(mp_fr[k_rec], 1.0, 6)
                    self.log("\n")
                
                
                if self.calc_xspec_with_plck_map:
                    self.log("# Average ratio of SxP / PxP cross spectra:")
                    self.log("")
                    
                    k_rec = self.av_xspec_ratio_key
                    mp_fr[k_rec] = self.avgs_xspec_ratios[map_id]
                    self.print_mapmapdouble(mp_fr[k_rec], 1.0, 3)
                    self.log("\n")
                
                
                if self.calc_noise_from_individ_maps:
                    self.log("# Noise level from individual maps [uK.arcmin]:")
                    self.log("")
                    
                    k_rec = self.indvid_noise_key
                    mp_fr[k_rec] = self.noise_from_individual_maps
                    du = core.G3Units.uK * core.G3Units.arcmin
                    self.print_mapmapdouble(mp_fr[k_rec], du, 3)
                
                
                if self.calc_noise_from_coadded_maps:
                    self.log("# Noise level from coadded maps [uK.arcmin]:")
                    self.log("")
                    
                    k_rec = self.coadds_noise_key
                    mp_fr[k_rec] = self.noise_from_individual_maps
                    du = core.G3Units.uK * core.G3Units.arcmin
                    self.print_mapmapdouble(mp_fr[k_rec], du, 3)
                
                
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
    
    thresholds = {"90GHz" : {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":  65.0*wu},
                             "ra0hdec-52.25": {"lo": 0.0*wu, "hi":  80.0*wu},
                             "ra0hdec-59.75": {"lo": 0.0*wu, "hi":  80.0*wu},
                             "ra0hdec-67.25": {"lo": 0.0*wu, "hi":  95.0*wu}},
                  "150GHz": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi": 120.0*wu},
                             "ra0hdec-52.25": {"lo": 0.0*wu, "hi": 140.0*wu},
                             "ra0hdec-59.75": {"lo": 0.0*wu, "hi": 145.0*wu},
                             "ra0hdec-67.25": {"lo": 0.0*wu, "hi": 165.0*wu}},
                  "220GHz": {"ra0hdec-44.75": {"lo": 0.0*wu, "hi":   8.0*wu},
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
    
    def __init__(self, map_ids, t_only,
                 bad_map_list_file, logging_function):
        
        self.t_only = t_only
        self.sb_fld = None
        self.obs_id = None
        self.mp_ids = map_ids
        self.x_list = bad_map_list_file
        self.logfun = logging_function
    
    def __call__(frame):
        
        if frame.type == core.G3FrameType.Observation:
            self.sb_field = frame["SourceName"]
        
        if (frame.type == core.G3FrameType.Map) and \
           (frame["Id"] in self.mp_ids):
            self.logfun("")
            self.logfun("* A relevant map frame has arrived.")
            self.logfun("* Checking whether the map is reasonable ...")
            
            map_is_bad = False
            
            flc_dct = calculate_map_fluctuation_metrics(frame, self.t_only)
            mean_wt = flc_dct["MeansOfWeights"]
            bad_weights = do_weights_look_bad(
                              mean_wt, frame["Id"], self.sb_fld)
            if bad_weights:
                record_bad_obs_id(self.x_list,
                                  frame["Id"], self.sb_fld, self.obs_id,
                                  "Anomalously large weights.")
            
            map_is_bad = bad_weights
            
            if map_is_bad:
                self.logfun("* ... the map looks bad!")
            else:
                self.logfun("* ... the map doesn't look obviously bad!")




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
                    self.mod_history[sb_fld][frame["Id"]]["left"]  = n_prev_mod
                elif "Right" in frame["Id"]:
                    self.mod_history[sb_fld][frame["Id"]]["right"] = n_prev_mod
        
        if frame.type == core.G3FrameType.Observation:
            self.sb_fld = frame["SourceName"]
        
        if (frame.type == core.G3FrameType.Map) and \
           (frame["Id"] in self.ids_to_modify):
            self.logfun("")
            self.logfun("Assigning a fake ID to the map ...")
            
            old_id = frame.pop("Id")
            if "IgnoreThisMap" in frame.keys():
                self.logfun("* (This map seems to be a bad one!")
                self.logfun("*  A fake ID will be assigned anyway,")
                self.logfun("*  but it will not be used for coadding.)")
                frame["Id"] = "Left" + old_id
                self.logfun("* Done.")
                self.logfun("")
            else:
                n_mp_l = self.mod_history[self.sb_fld][frame["Id"]]["left"]
                n_mp_r = self.mod_history[self.sb_fld][frame["Id"]]["right"]
                if n_mp_l > n_mp_r:
                    frame["Id"] = "Right" + old_id
                    self.mod_history[self.sb_fld][frame["Id"]]["right"] += 1
                else:
                    frame["Id"] = "Left" + old_id
                    self.mod_history[self.sb_fld][frame["Id"]]["left"] += 1
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
        rmss_and_wgts_from_coadds_or_individuals="",
        rmss_and_wgts_from_signals_or_noises="",
        calculate_noise_from_individual_maps=False,
        calculate_noise_from_coadded_maps=False,
        point_source_list_file="",
        calculate_cross_spectra_with_planck_map=False,
        planck_map_fits_files=[],
        calculate_pointing_discrepancies=False,
        logger_name="", bad_map_list_file='', log_file=None):
    
    
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
        pipeline.Add(FlagBadMaps)
        pipeline.Add(AppendDirectionsToMapIDs,
                     sub_fields=sub_fields,
                     ids_to_modify=ids_to_append_directions_to,
                     t_only=temperature_maps_only,
                     bad_map_list_file=bad_map_list_file,
                     logging_function=log)
    
    pipeline.Add(lambda frame: log(frame))
    
    pipeline.Add(AnalyzeAndCoaddMaps,
                 map_ids=map_ids, map_sources=sub_fields,
                 temperature_only=temperature_maps_only,
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
                 planck_map_fits_files=planck_map_fits_files,
                 calculate_pointing_discrepancies=\
                     calculate_pointing_discrepancies,
                 logging_function=log)
    
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
    
    parser.add_argument("-x", "--calculate_cross_spectra_with_planck_map",
                        action="store_true", default=False,
                        help="Whether to calculate the cross spectra between "
                             "individual maps and a Planck map.")
    
    parser.add_argument("-X", "--planck_map_fits_files",
                        type=str, action="store", nargs="+", default=[],
                        help="Path to a FITS file that contains a Planck map.")
    
    parser.add_argument("-p", "--calculate_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to calculate the difference between the "
                             "true positions of some point sources and their "
                             "measured positions.")
    
    parser.add_argument("-g", "--logger_name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "
                             "record log messages.")
    
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

