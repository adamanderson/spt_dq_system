## *******  Always need to import modules!  ******* ##

from    spt3g              import  core
from    spt3g              import  mapmaker
from    spt3g              import  util
from    spt3g              import  coordinateutils
from    spt3g              import  std_processing
from    spt3g              import  mapspectra
from    spt3g.mapspectra   import  map_analysis
from    scipy              import  signal
import  os
import  gc
import  sys
import  glob
import  numpy
import  argparse
import  datetime
import  logging




# ==============================================================================
# Define modules and functions needed for the pipeline
# ------------------------------------------------------------------------------


def map_seems_fine(map_frame, log_fun):
    
    if map_frame["T"].units != core.G3TimestreamUnits.Tcmb:
        log_fun("The units of this map are not in Tcmb!")
        return False
    
    map_values = map_frame["T"]
    if not numpy.isfinite(numpy.nanmean(numpy.asarray(map_values))):
        log_fun("")
        log_fun("There seem to be only NaNs in the map!")
        return False
        
    return True




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




def add_map_frames(
        frame_of_coadded_maps, new_frame,
        map_keys, operation, ignore_weights=True):
    
    for from_this_key, to_that_key in map_keys.items():
        if to_that_key not in frame_of_coadded_maps.keys():
            frame_of_coadded_maps[to_that_key] = new_frame[from_this_key]
        else:
            existing_map = frame_of_coadded_maps.pop(to_that_key)
            frame_of_coadded_maps[to_that_key] = \
                existing_map + new_frame[from_this_key] * operation
        
    if ignore_weights:
        pass
    elif "Wunpol" in new_frame.keys():
        if "Wunpol" not in frame_of_coadded_maps.keys():
            frame_of_coadded_maps["Wunpol"] = new_frame["Wunpol"]
        else:
            frame_of_coadded_maps["Wunpol"].TT += \
                new_frame["Wunpol"].TT * operation
    
    return frame_of_coadded_maps




def calculate_map_fluctuation_metrics(
        map_frame, temperature_only,
        center_ra, center_dec,
        return_nans=False):
    
    keys = ["MapStdDevs",
            "MeansOfWeights",
            "NumbersOfPixelsWithGoodWeights",
            "VarsOfProductsOfMapValuesAndSqrtWeights"]
    
    fluctuation_metrics = {key: numpy.nan for key in keys}
    
    if return_nans:
        return fluctuation_metrics
    else:
        if temperature_only:
            t_vals = numpy.asarray(mapmaker.mapmakerutils.remove_weight_t(
                                   map_frame["T"], map_frame["Wunpol"]))
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
    



def combine_ids(old_dict, new_dict):
    
    this_frame_already_added = False
    for map_id in new_dict.keys():
        for source, new_ids in new_dict[map_id].items():
            if source not in old_dict[map_id].keys():
                old_dict[map_id][source] = new_ids
            else:
                old_ids    = old_dict[map_id][source]
                common_ids = set(new_ids) & set(old_ids)
                if len(common_ids) > 0:
                    this_frame_already_added = True
                else:
                    for new_id in new_ids:
                        old_dict[map_id][source].append(new_id)
    
    return old_dict, this_frame_already_added




def remove_ids(old_dict, new_dict):
    
    for map_id in new_dict.keys():
        for source, ids_to_remove in new_dict[map_id].items():
            updated_ids = []
            for old_id in old_dict[map_id][source]:
                if old_id in ids_to_remove:
                    pass
                else:
                    updated_ids.append(old_id)
            updated_ids = type(old_dict[map_id][source])(updated_ids)
            old_dict[map_id].pop(source)
            old_dict[map_id][source] = updated_ids
    
    return old_dict




def combine_mapmapdoubles(old_dict, new_dict):
    
    for map_id in new_dict.keys():
        for source, some_map in new_dict[map_id].items():
            if source not in old_dict[map_id].keys():
                old_dict[map_id][source] = some_map
            else:
                for obs_id, some_value in some_map.items():
                    if obs_id not in old_dict[map_id][source].keys():
                        old_dict[map_id][source][obs_id] = some_value
    
    return old_dict




def remove_partial_mapmapdouble(old_dict, new_dict):
    
    for map_id in new_dict.keys():
        for source, obs_ids_to_remove in new_dict[map_id].items():
            for obs_id in obs_ids_to_remove:
                old_dict[map_id][source].pop(str(obs_id))
    
    return old_dict




class CoaddMapsAndDoSomeMapAnalysis(object):
    
    def __init__(self, map_sources, map_ids, temperature_only,
                 maps_split_by_scan_direction=False,
                 combine_left_right=False,
                 combine_different_wafers=False,
                 allow_subtraction=False,
                 collect_averages_from_flagging_stats=False,
                 calculate_pW_to_K_conversion_factors=False,
                 calculate_map_rms_and_weight_stat=False,
                 calculate_noise_from_individual_maps=False,
                 calculate_noise_from_coadded_maps=False,
                 calculate_xspec_with_coadded_maps=False,
                 point_source_file=None,
                 bad_map_list_file=None,
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
            self.log("- An ID stored in a map frame will be regarded")
            self.log("  as a different ID according to the relations below:")
            for id_from in sorted(self.id_mapping.keys()):
                id_to = self.id_mapping[id_from]
                self.log("      %s ==> %s", id_from, id_to)
        
        self.bad_map_list_file = bad_map_list_file

        
        # - Initialize variables related to storing coadded maps and
        # - their basic information
        
        self.map_sources = map_sources
        self.temp_only   = temperature_only
        
        self.coadded_map_frames = {map_id: core.G3Frame(core.G3FrameType.Map) \
                                   for map_id in self.map_ids}
        self.coadded_obs_ids    = {map_id: core.G3MapVectorInt()    \
                                   for map_id in self.map_ids}
        self.coadded_map_ids    = {map_id: core.G3MapVectorString() \
                                   for map_id in self.map_ids}
        self.observat_durations = {map_id: core.G3MapMapDouble()    \
                                   for map_id in self.map_ids}
        self.allow_subtraction  = allow_subtraction
        self.obs_info = None
        self.calframe = None
        
        
        # - Initialize variables related to collecting average numbers
        # - of various bolometer flagging statistics
        
        self.get_avgs_flagging_stats = collect_averages_from_flagging_stats
        
        if self.get_avgs_flagging_stats:
            self.bolo_props       = None
            self.pipe_line_info   = None
            self.flagging_reasons = \
                ["BadCalSn", "BadWeight",  "Glitchy", "Oscillating",
                 "Latched",  "Overbiased", "BadHk",
                 "PostCalibrationNaNs",
                 "UnphysicalLowVariance",
                 "MissingFluxCalibration",
                 "Others",  "TotalNotFlagged", "TotalRemoved"]
            
            self.avgs_flagging_stats =  \
                {flagging_related_name: \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
                 for flagging_related_name in self.flagging_reasons}
        
        
        # - Initialize variables related to calculating the median value of
        # - the pW/K conversion factors for each band of each wafer
        
        self.calc_mean_pW_to_K = calculate_pW_to_K_conversion_factors
        
        if self.calc_mean_pW_to_K:
            self.wafers = ["W172", "W174", "W176", "W177", "W180",
                           "W181", "W188", "W203", "W204", "W206"]
            self.means_of_temp_cal_factors = \
                {wafer: {mid: core.G3MapMapDouble() \
                         for mid in self.map_ids} for wafer in self.wafers}
        
        
        # - Initialize variables related to calculating a Std. dev. of
        # - a map, a mean weight, and a variance of a product between
        # - map values and square root of the weights
        
        self.calc_map_stddev_etc = calculate_map_rms_and_weight_stat
        
        if self.calc_map_stddev_etc:
            self.fluctuation_keys = ["MapStdDevs",
                                     "MeansOfWeights",
                                     "NumbersOfPixelsWithGoodWeights",
                                     "VarsOfProductsOfMapValuesAndSqrtWeights"]
            self.map_types = ["IndividualMaps",
                              "CoaddedMaps"]
            self.map_fluctuation_metrics =  \
                {map_type: {fluc_type: {mid: core.G3MapMapDouble()  \
                                        for mid in self.map_ids}    \
                            for fluc_type in self.fluctuation_keys} \
                 for map_type in self.map_types}
        
        
        # - Initialize variables related to pointing dicrepancy calculations
        # - but not noise nor cross spectra calculations
        
        self.calc_pointing_discrepancies = calculate_pointing_discrepancies
        
        if self.calc_pointing_discrepancies:
            self.delta_ras_and_decs = \
                {rank: {delta_type: {mid: core.G3MapMapDouble()       \
                                     for mid in self.map_ids}         \
                        for delta_type in ["delta_ra", "delta_dec"]}  \
                 for rank in ["1", "2", "3"]}
            self.fluxes = \
                {rank: {mid: core.G3MapMapDouble() for mid in self.map_ids} \
                 for rank in ["1", "2", "3"]}
            self.snrs = \
                {rank: {mid: core.G3MapMapDouble() for mid in self.map_ids} \
                 for rank in ["1", "2", "3"]}
        
        
        # - Initialize variables related to noise and croass spectra 
        # - calculations but not pointing discrepancy calculations
        
        self.calc_noise_from_individ_maps = calculate_noise_from_individual_maps
        self.calc_noise_from_coadded_maps = calculate_noise_from_coadded_maps
        self.calc_xspec_with_coadded_maps = calculate_xspec_with_coadded_maps
        self.point_source_file            = point_source_file
        self.maps_split_by_scan_direction = maps_split_by_scan_direction
        
        if self.calc_noise_from_individ_maps:
            self.noise_from_individual_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        elif self.calc_noise_from_coadded_maps:
            self.noise_from_coadded_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            if not self.maps_split_by_scan_direction:
                self.operations_done_to_maps = \
                    {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        if self.calc_xspec_with_coadded_maps:
            self.xspec_cls_with_coadds = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        
        
        # - Initialize variables related to
        # - noise, cross spectra, and pointing discrepancy calculations
                
        if (self.maps_split_by_scan_direction and \
            (self.calc_noise_from_individ_maps or \
             self.calc_noise_from_coadded_maps or \
             self.calc_pointing_discrepancies)):
            
            self.direction_independent_map_ids = []
            self.map_fr_arrival_counters       = {}
            for map_id in self.map_ids:
                if "Left" in map_id:
                    self.direction_independent_map_ids.append(
                        map_id.replace("Left", "SomeDirection"))
                    self.map_fr_arrival_counters[map_id] = 0
                elif "Right" in map_id:
                    self.direction_independent_map_ids.append(
                        map_id.replace("Right", "SomeDirection"))
                    self.map_fr_arrival_counters[map_id] = 0
                else:
                    self.direction_independent_map_ids.append(map_id)
                    self.map_fr_arrival_counters["Left"+map_id]  = 0
                    self.map_fr_arrival_counters["Right"+map_id] = 0
            self.direction_independent_map_ids = \
                set(self.direction_independent_map_ids)
            
            if self.calc_noise_from_individ_maps or \
               self.calc_pointing_discrepancies:
                self.map_frame_cache = \
                    {map_id: None for map_id in self.map_ids}
    
    
    
    def __call__(self, frame):
        
        if ("CoaddedMapsContained" in frame.keys()) and \
           (frame["CoaddedMapsContained"] == False):
            return []   # Some clean-up of non-useful data created previously
        
        if frame.type == core.G3FrameType.Calibration:
            try:
                self.bolo_props = frame["BolometerProperties"]
                self.calframe   = frame
            except KeyError:
                pass
            return
        
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
            return
        
        if frame.type == core.G3FrameType.PipelineInfo:
            if ["DroppedBolos", "SurvivingBolos"] <= frame.keys():
                self.pipe_info = frame
            return
        
        
        if frame.type == core.G3FrameType.Map:
            
            # - Load a map frame and figure out
            # - whether its data should go into the output file
            # - based on the ID
                        
            id_for_coadds = None
            
            # ** The above will be the map ID that identifies
            # ** what this frame is about
            # ** when coaddition and other processing is performed
            
            frame_has_coadds = ("CoaddedMapsContained" in frame.keys())
            
            if frame_has_coadds:
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
                    self.log("* because the observation type is")
                    self.log("* not of interest.")
                    self.log("\n")
                    return []
                else:
                    id_for_coadds = self.id_mapping[frame["Id"]]
                    oid  = self.obs_info["ObservationID"]
                    src  = self.obs_info["SourceName"]
                    dmid = str(oid) + frame["Id"]   # ** more detailed map id
                    t_i  = self.obs_info["ObservationStart"]
                    t_f  = self.obs_info["ObservationStop"]
                    d_i  = std_processing.time_to_obsid(t_i)
                    d_f  = std_processing.time_to_obsid(t_f)
                    dura = (d_f - d_i) * core.G3Units.s
                    center_ra  = core.G3Units.deg * 0.0
                    center_dec = core.G3Units.deg * float(src[-6:])
                    obs_ids_from_this_frame = \
                        {id_for_coadds: {src: core.G3VectorInt([oid])}}
                    map_ids_from_this_frame = \
                        {id_for_coadds: {src: core.G3VectorString([dmid])}}
                    obs_tms_from_this_frame = \
                        {id_for_coadds: {src: core.G3MapMapDouble()}}
                    obs_tms_from_this_frame \
                        [id_for_coadds][src] = core.G3MapDouble()
                    obs_tms_from_this_frame \
                        [id_for_coadds][src][str(oid)] = dura
                    
            
            
            # - Perform a few basic checks on the frame
            # - to furthur figure out whether it can go into the coadds
            
            if not map_seems_fine(frame, self.log):
                self.log("")
                self.log("* Well, the map doesn't look good,")
                self.log("* so, this one will be skipped ...")
                self.log("\n")
                return []
            
            self.coadded_obs_ids, this_frame_already_added = \
                combine_ids(self.coadded_obs_ids,
                            obs_ids_from_this_frame)
            self.coadded_map_ids, this_frame_already_added = \
                combine_ids(self.coadded_map_ids,
                            map_ids_from_this_frame)
            self.observat_durations = \
                combine_mapmapdoubles(self.observat_durations,
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

            
            # - Calculate mean value of the pW/K conversion factors
            # - for each band of each wafer
            
            if subtract_maps_in_this_frame:
                if self.calc_mean_pW_to_K:
                    for wafer in self.means_of_temp_cal_factors.keys():
                        self.means_of_temp_cal_factors[wafer] = \
                            remove_partial_mapmapdouble(
                                self.means_of_temp_cal_factors[wafer],
                                obs_ids_from_this_frame)
            
            elif self.calc_mean_pW_to_K:
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering the mean values of the pW/K")
                    self.log("* coversion factors calculated previously ...")
                    prefix = "MeansOfTemperatureCalibrationFactors"
                    means_of_convs_from_this_fr = {}
                    for key in frame.keys():
                        if prefix in key:
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
                    for band, wfs_and_avgs in avgs_to_be_reorganized.items():
                        if (band+"GHz") in id_for_coadds:
                            for wafer, avg in wfs_and_avgs.items():
                                mmd = core.G3MapMapDouble()
                                mmd[src] = core.G3MapDouble()
                                mmd[src][str(oid)] = avg
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
            
            
            # - Calculate cross spectrum between coadded maps and
            # - one observation's map
            
            if self.calc_xspec_with_coadded_maps:
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering avg.(sqrt(Cl)) of cross spectra")
                    self.log("* that were calculated previously ...")
                    cl_key = "ClsFromCrossSpectraBetweenCoaddsAndIndividuals"
                    try:
                        xspec_sqrt_cls_this_fr = {id_for_coadds: frame[cl_key]}
                    except KeyError:
                        xspec_sqrt_cls_this_fr = {id_for_coadds: core.G3MapMapDouble()}
                elif len(self.coadded_obs_ids[id_for_coadds][src]) <= 1:
                    self.log("")
                    self.log("* Currently, the maps in the cache")
                    self.log("* (or maps to be added to the cache) are")
                    self.log("* the same as the maps in this frame,")
                    self.log("* so the cross spectrum calculation")
                    self.log("* will be skipped for this frame.")
                    xspec_sqrt_cls = core.G3MapMapDouble()
                    xspec_sqrt_cls[src] = core.G3MapDouble()
                    xspec_sqrt_cls[src][str(oid)] = numpy.nan
                    xspec_sqrt_cls_this_fr = {id_for_coadds: xspec_sqrt_cls}
                else:
                    self.log("")
                    self.log("* Calculating avg.(sqrt(Cl)) of the")
                    self.log("* cross spectrum between the coadded maps")
                    self.log("* and the map from this observation ...")
                    xspec_sqrt_cl = calculate_average_cl_of_cross_spectrum(
                                        self.coadded_map_frames[id_for_coadds],
                                        frame,
                                        self.point_source_file,
                                        temp_only=self.temp_only,
                                        smaller_region=True,
                                        center_ra=center_ra,
                                        center_dec=center_dec)
                    xsc = xspec_sqrt_cl/(core.G3Units.uK*core.G3Units.arcmin)
                    self.log("* ... the average was calculated to be")
                    self.log("* %s uK.arcmin.", xsc)
                    xspec_sqrt_cls = core.G3MapMapDouble()
                    xspec_sqrt_cls[src] = core.G3MapDouble()
                    xspec_sqrt_cls[src][str(oid)] = xspec_sqrt_cl
                    xspec_sqrt_cls_this_fr = {id_for_coadds: xspec_sqrt_cls}
                self.log("* Done.")
                self.log("")
                
                self.xspec_cls_with_coadds = \
                    combine_mapmapdoubles(
                        self.xspec_cls_with_coadds,
                        xspec_sqrt_cls_this_fr)
            
            
            # - Add (including -1 of) maps together
                        
            if subtract_maps_in_this_frame:
                self.log("")
                self.log("* Removing a set of ")
                self.log("* sky map and weight map (%s)", frame["Id"])
                self.log("* from the coadded maps for %s ...", id_for_coadds)
                self.coadded_obs_ids = \
                    remove_ids(self.coadded_obs_ids, obs_ids_from_this_frame)
                self.coadded_map_ids = \
                    remove_ids(self.coadded_map_ids, map_ids_from_this_frame)
                self.observat_durations = \
                    remove_partial_mapmapdouble(self.observat_durations,
                                                obs_ids_from_this_frame)
                if self.temp_only:
                    self.coadded_map_frames[id_for_coadds] = \
                        add_map_frames(self.coadded_map_frames[id_for_coadds],
                                       frame,
                                       {"T": "T"}, -1.0, ignore_weights=False)
                self.log("* Done.")
                self.log("")
            else:
                self.log("")
                self.log("* Adding a set of ")
                self.log("* sky map and weight map (ID: %s)", frame["Id"])
                self.log("* to the coadded maps (ID: %s) ...", id_for_coadds)
                if self.temp_only:
                    self.coadded_map_frames[id_for_coadds] = \
                        add_map_frames(self.coadded_map_frames[id_for_coadds],
                                       frame,
                                       {"T": "T"}, 1.0, ignore_weights=False)
                self.log("* Done.")
                self.log("")
            
            
            # - Calculate a standard deviation of the map values,
            # - a mean of the weights, and a variance of the product
            # - of the map values and square root of the weights
            
            if subtract_maps_in_this_frame:
                if self.calc_map_stddev_etc:
                    for mtp in self.map_types:
                        for fk in self.fluctuation_keys:
                            try:
                                self.map_fluctuation_metrics[mtp][fk] = \
                                    remove_partial_mapmapdouble(
                                        self.map_fluctuation_metrics[mtp][fk],
                                        obs_ids_from_this_frame)
                            except KeyError:
                                pass
            
            elif self.calc_map_stddev_etc:
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering map standard deviations,")
                    self.log("* mean values of weights, and so on")
                    self.log("* that were calculated previously ...")
                    stddevs_etc = {}
                    prefix = "BasicMapFluctuationMetrics"
                    for mtp in self.map_types:
                        stddevs_etc[mtp] = {}
                        for fk in self.fluctuation_keys:
                            try:
                                stddevs_etc[mtp][fk] = \
                                    {id_for_coadds: frame[prefix+mtp+fk]}
                            except KeyError:
                                stddevs_etc[mtp][fk] = \
                                    {id_for_coadds: core.G3MapMapDouble()}
                            
                else:
                    self.log("")
                    self.log("* Gathering a map standard deviation,")
                    self.log("* a mean value of weights, and so on")
                    self.log("* from this observation's map and also ")
                    self.log("* the coadded maps ...")
                    stddevs_etc = {}
                    for mtp in self.map_types:
                        if mtp == "IndividualMaps":
                            fr_to_use = frame
                        elif mtp == "CoaddedMaps":
                            fr_to_use = self.coadded_map_frames[id_for_coadds]
                        stddevs_etc[mtp] = {}
                        flc_mtrs = calculate_map_fluctuation_metrics(
                                       fr_to_use,
                                       self.temp_only,
                                       center_ra, center_dec,
                                       return_nans=(self.allow_subtraction and\
                                                    mtp=="CoaddedMaps"))
                        for flc_k, flc_vl in flc_mtrs.items():
                            mmd = core.G3MapMapDouble()
                            mmd[src] = core.G3MapDouble()
                            mmd[src][str(oid)] = flc_vl
                            stddevs_etc[mtp][flc_k] = \
                                {id_for_coadds: mmd}
                self.log("* Done.")
                self.log("")
                for mtp in self.map_types:
                    for fk in self.fluctuation_keys:
                        self.map_fluctuation_metrics[mtp][fk] = \
                            combine_mapmapdoubles(
                                self.map_fluctuation_metrics[mtp][fk],
                                stddevs_etc[mtp][fk])
            
            
            # - Collect averages related to flagging statistics
            
            if subtract_maps_in_this_frame:
                if self.get_avgs_flagging_stats:
                    for flg_typ in self.avgs_flagging_stats.keys():
                        self.avgs_flagging_stats[flg_typ] = \
                            remove_partial_mapmapdouble(
                                self.avgs_flagging_stats[flg_typ],
                                obs_ids_from_this_frame)
            
            elif self.get_avgs_flagging_stats:
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering average numbers related to")
                    self.log("* flagging statistics that were")
                    self.log("* calculated previously ...")
                    prefix = "FlaggingStatisticsAverageNumberOf"
                    avgs_flg_stats_from_this_fr = {}
                    for key in frame.keys():
                        if prefix in key:
                            flg_typ = key.replace(prefix, "")
                            avgs_flg_stats_from_this_fr[flg_typ] = \
                                {id_for_coadds: frame[key]}
                else:
                    self.log("")
                    self.log("* Gathering average numbers related to")
                    self.log("* flagging statistics ...")
                    avgs_flg_stats_from_this_fr = {}
                    avgs_to_be_reorganized = \
                        collect_averages_from_flagging_info(
                            self.pipe_info,
                            self.bolo_props,
                            self.flagging_reasons)
                    for band, flg_typs_avgs in avgs_to_be_reorganized.items():
                        if (band+"GHz") in id_for_coadds:
                            for flg_typ, avg in flg_typs_avgs.items():
                                mmd = core.G3MapMapDouble()
                                mmd[src] = core.G3MapDouble()
                                mmd[src][str(oid)] = avg
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
            
            
            # - Prepare data to be used for pointing discrepancy
            # - and noise level calculations
            
            if subtract_maps_in_this_frame:
                if self.calc_pointing_discrepancies:
                    for rank in ["1", "2", "3"]:
                        for delta_type in ["delta_ra", "delta_dec"]:
                            self.delta_ras_and_decs[rank][delta_type] = \
                                remove_partial_mapmapdouble(
                                    self.delta_ras_and_decs[rank][delta_type],
                                    obs_ids_from_this_frame)
                        self.fluxes[rank] = \
                            remove_partial_mapmapdouble(
                                self.fluxes[rank], obs_ids_from_this_frame)
                        self.snrs[rank] = \
                            remove_partial_mapmapdouble(
                                self.snrs[rank], obs_ids_from_this_frame)
                
                if self.calc_noise_from_individ_maps:
                    self.noise_from_individual_maps = \
                        remove_partial_mapmapdouble(
                            self.noise_from_individual_maps,
                            obs_ids_from_this_frame)
                
                time_to_calculate_pointing_and_noise = False
            
            elif self.calc_noise_from_coadded_maps or \
                 self.calc_noise_from_individ_maps or \
                 self.calc_pointing_discrepancies:
                
                if frame_has_coadds:
                    time_to_calculate_pointing_and_noise = True
                    # ** The assignment above is just a dummy one
                
                elif self.maps_split_by_scan_direction:
                    if self.calc_noise_from_individ_maps or \
                       self.calc_pointing_discrepancies:
                        self.map_frame_cache[frame["Id"]] = frame
                    self.map_fr_arrival_counters[frame["Id"]] += 1
                    
                    time_to_calculate_pointing_and_noise = False
                    for di_mid in self.direction_independent_map_ids:
                        if "SomeDirection" in di_mid:
                            id_l_ver = di_mid.replace("SomeDirection", "Left")
                            id_r_ver = di_mid.replace("SomeDirection", "Right")
                        else:
                            id_l_ver = "Left"  + di_mid
                            id_r_ver = "Right" + di_mid
                        if self.map_fr_arrival_counters[id_l_ver] == \
                           self.map_fr_arrival_counters[id_r_ver]:
                            time_to_calculate_pointing_and_noise = True
                            break
                else:
                    time_to_calculate_pointing_and_noise = True
                
                if time_to_calculate_pointing_and_noise and \
                   (not frame_has_coadds):
                    if self.combine_left_right:
                        ids_for_recording = [id_for_coadds]
                    elif (not self.maps_split_by_scan_direction):
                        ids_for_recording = [id_for_coadds]
                    else:
                        ids_for_recording = [id_r_ver, id_l_ver]
            else:
                time_to_calculate_pointing_and_noise = False
            
            
            # - Calculate pointing discrepancies
            
            if time_to_calculate_pointing_and_noise and \
               self.calc_pointing_discrepancies:
                
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering pointing discrepancy and flux data")
                    self.log("* that were calculated previously ...")
                    for rank in ["1", "2", "3"]:
                        for coordinate in ["Ra", "Dec"]:
                            delta_type   = "delta_" + coordinate.lower()
                            relevant_key = "Delta"  + coordinate + "s" + \
                                           "FromSources" + rank
                            self.delta_ras_and_decs[rank][delta_type] = \
                                combine_mapmapdoubles(
                                    self.delta_ras_and_decs[rank][delta_type],
                                    {id_for_coadds: frame[relevant_key]})
                        flux_key = "IntegratedFluxesFromSources" + rank
                        self.fluxes[rank] = \
                            combine_mapmapdoubles(
                                self.fluxes[rank],
                                {id_for_coadds: frame[flux_key]})
                        snr_key = "PointSourceSNRs" + rank
                        self.snrs[rank] = \
                            combine_mapmapdoubles(
                                self.snrs[rank],
                                {id_for_coadds: frame[snr_key]})
                    self.log("* Done.")
                    self.log("")
                else:
                    if self.calc_map_stddev_etc:
                        map_stddev = self.map_fluctuation_metrics\
                                     ["IndividualMaps"]["MapStdDevs"]\
                                     [id_for_coadds][src][str(oid)]
                    else:
                        map_stddev = None
                    if self.maps_split_by_scan_direction:
                        self.log("")
                        self.log("* Pointing discrepancy and flux calculations")
                        self.log("* are about to start because both")
                        self.log("* %s and %s", id_l_ver, id_r_ver)
                        self.log("* of obs. %s", oid)
                        self.log("* have passed the pipeline ...")
                        
                        pt_offsets_to_be_reorganized,   \
                        fluxes_to_be_reorganized,       \
                        source_snrs_to_be_reorganized = \
                            calculate_pointing_discrepancies(
                                self.map_frame_cache[id_r_ver],
                                self.map_frame_cache[id_l_ver],
                                src,
                                temp_only=True,
                                map_stddev=map_stddev)
                    else:
                        self.log("")
                        self.log("* Pointing discrepancy and flux calculations")
                        self.log("* are about to start for")
                        self.log("* obs. %s %s ...", oid, id_for_coadds)
                        
                        pt_offsets_to_be_reorganized,   \
                        fluxes_to_be_reorganized,       \
                        source_snrs_to_be_reorganized = \
                            calculate_pointing_discrepancies(
                                frame, None, src,
                                temp_only=True, map_stddev=map_stddev)
                        
                    for rank, diffs in pt_offsets_to_be_reorganized.items():
                        for delta_coord, diff in diffs.items():
                            self.log("* Source %s " + \
                                     delta_coord.ljust(9) + \
                                     " %s arcseconds",
                                     rank, diff/core.G3Units.arcsec)
                        flux  = fluxes_to_be_reorganized[rank]
                        usys  = core.G3Units
                        flux /= (usys.mK * usys.arcmin * usys.arcmin)
                        snr   = source_snrs_to_be_reorganized[rank]
                        self.log("* Source %s "  + \
                                 "flux".ljust(9) + \
                                 " %s mK.arcmin.arcmin",
                                 rank, flux)
                        self.log("* Source %s " + \
                                 "SNR".ljust(9) + \
                                 " %s.", rank, snr)
                    self.log("* Done.")
                    self.log("")
                    
                    for idr in ids_for_recording:
                        pt_offsets_from_this_fr = \
                            {rank: {d_type: {idr: core.G3MapMapDouble()} \
                             for d_type in ["delta_ra", "delta_dec"]}    \
                             for rank   in ["1", "2", "3"]}
                        fluxes_from_this_fr = \
                            {rank: {idr: core.G3MapMapDouble()} \
                             for rank in ["1", "2", "3"]}
                        snrs_from_this_fr = \
                            {rank: {idr: core.G3MapMapDouble()} \
                             for rank in ["1", "2", "3"]}
                        
                        for rank in ["1", "2", "3"]:
                            for d_type in ["delta_ra", "delta_dec"]:
                                pt_offsets_from_this_fr \
                                [rank][d_type][idr][src] = core.G3MapDouble()
                                pt_offsets_from_this_fr \
                                [rank][d_type][idr][src][str(oid)] = \
                                    pt_offsets_to_be_reorganized[rank][d_type]
                                
                                self.delta_ras_and_decs[rank][d_type] = \
                                    combine_mapmapdoubles(
                                        self.delta_ras_and_decs[rank][d_type],
                                        pt_offsets_from_this_fr[rank][d_type])
                        
                        for rank in ["1", "2", "3"]:
                            fluxes_from_this_fr \
                            [rank][idr][src] = core.G3MapDouble()
                            fluxes_from_this_fr \
                            [rank][idr][src][str(oid)] = \
                                fluxes_to_be_reorganized[rank]
                            
                            snrs_from_this_fr \
                            [rank][idr][src] = core.G3MapDouble()
                            snrs_from_this_fr \
                            [rank][idr][src][str(oid)] = \
                                source_snrs_to_be_reorganized[rank]
                            
                            self.fluxes[rank] = \
                                combine_mapmapdoubles(
                                    self.fluxes[rank],
                                    fluxes_from_this_fr[rank])
                            
                            self.snrs[rank] = \
                                combine_mapmapdoubles(
                                    self.snrs[rank],
                                    snrs_from_this_fr[rank])
            
            
            # - Calculate noise levels
            
            if (time_to_calculate_pointing_and_noise and \
                (self.calc_noise_from_coadded_maps or \
                 self.calc_noise_from_individ_maps)):
                
                if frame_has_coadds:
                    self.log("")
                    self.log("* Gathering noise data")
                    self.log("* that were calculated previously ...")
                    if self.calc_noise_from_individ_maps:
                        noise_key = "NoiseFromIndividualMaps"
                        self.noise_from_individual_maps = \
                            combine_mapmapdoubles(
                                self.noise_from_individual_maps,
                                {id_for_coadds: frame[noise_key]})
                    else:
                        noise_key = "NoiseFromCoaddedMaps"
                        self.noise_from_coadded_maps = \
                            combine_mapmapdoubles(
                                self.noise_from_coadded_maps,
                                {id_for_coadds: frame[noise_key]})
                        if self.maps_split_by_scan_direction:
                            self.map_fr_arrival_counters[id_for_coadds] = \
                                len(frame["CoaddedObservationIDs"].values()[0])
                        else:
                            if self.temp_only:
                                self.coadded_map_frames[id_for_coadds] = \
                                    add_map_frames(
                                        self.coadded_map_frames[id_for_coadds],
                                        frame,
                                        {"T_Noise": "T_Noise"}, 1.0)
                            ops_key = "NoiseCalculationsOperationsDoneToMaps"
                            self.operations_done_to_maps = \
                                combine_mapmapdoubles(
                                    self.operations_done_to_maps,
                                    {id_for_coadds: frame[ops_key]})
                    self.log("* Done.")
                    self.log("")
                else:
                    noise_units = core.G3Units.uK*core.G3Units.arcmin
                    if self.maps_split_by_scan_direction:
                        self.log("")
                        self.log("* Noise calculation is about to start")
                        self.log("* b/c both %s and %s", id_l_ver, id_r_ver)
                        self.log("* of obs. %s", oid)
                        self.log("* have passed the pipeline ...")
                        
                        if self.calc_noise_from_individ_maps:
                            noise = calculate_noise_level(
                                        self.map_frame_cache[id_r_ver],
                                        self.map_frame_cashe[id_l_ver],
                                        self.point_source_file,
                                        temp_only=self.temp_only,
                                        smaller_region=True,
                                        center_ra=center_ra,
                                        center_dec=center_dec)
                        else:
                            noise = calculate_noise_level(
                                        self.coadded_map_frames[id_r_ver],
                                        self.coadded_map_frames[id_l_ver],
                                        self.point_source_file,
                                        temp_only=self.temp_only,
                                        smaller_region=True,
                                        center_ra=center_ra,
                                        center_dec=center_dec)
                    else:
                        self.log("")
                        self.log("* Noise calculation is about to start for")
                        self.log("* obs. %s %s ...", oid, id_for_coadds)
                        
                        if self.calc_noise_from_individ_maps:
                            noise = calculate_noise_level(
                                        frame,
                                        None,
                                        self.point_source_file,
                                        temp_only=self.temp_only,
                                        smaller_region=True,
                                        center_ra=center_ra,
                                        center_dec=center_dec)

                        else:
                            if src not in \
                            self.operations_done_to_maps[id_for_coadds].keys():
                                self.operations_done_to_maps\
                                [id_for_coadds][src] = core.G3MapDouble()
                            
                            new_frame = core.G3Frame(core.G3FrameType.Map)
                            if self.temp_only:
                                new_frame["T"] = \
                                    mapmaker.mapmakerutils.remove_weight_t(
                                        frame["T"], frame["Wunpol"])
                                new_frame["Wunpol"] = frame["Wunpol"]
                            
                            operation = decide_operation_to_do_with_new_map(
                                            new_frame,
                                            self.operations_done_to_maps\
                                            [id_for_coadds][src],
                                            self.coadded_obs_ids\
                                            [id_for_coadds][src],
                                            self.temp_only,
                                            center_ra, center_dec)
                            
                            op_dict = {"Copy"    :  1.0,
                                       "Add"     :  1.0,
                                       "Subtract": -1.0,
                                       "Ignore"  :  0.0}
                            
                            self.operations_done_to_maps\
                            [id_for_coadds][src][str(oid)] = op_dict[operation]                            
                            self.log("* (The operation to be performed")
                            self.log("*  on the map that is stored")
                            self.log("*  in this frame is %s.)", operation)
                            
                            if operation == "Ignore":
                                pass
                            else:
                                if self.temp_only:
                                    self.coadded_map_frames[id_for_coadds] = \
                                        add_map_frames(
                                            self.coadded_map_frames\
                                            [id_for_coadds],
                                            new_frame,
                                            {"T": "T_Noise"},
                                            op_dict[operation])
                                self.log("* (The operation was done.)")
                            
                            n_added = list(self.operations_done_to_maps\
                                         [id_for_coadds][src].values()).\
                                      count(1.0)
                            n_subed = list(self.operations_done_to_maps\
                                         [id_for_coadds][src].values()).\
                                      count(-1.0)
                            
                            if operation == "Ignore":
                                self.log("* (Since this map was bad "
                                         "in some way,")
                                self.log("*  and it was not added to "
                                         "nor subtracted from the "
                                         "running noise map,")
                                self.log("*  the noise calc. will not occur.")
                                self.log("*  -1 will be recorded as a dummy.)")
                                """record_bad_obs_id(
                                    self.bad_map_list_file,
                                    id_for_coadds, src, oid,
                                    "Bad values around field center.")"""
                                noise = -1.0 * noise_units                                
                            else:
                                divide_map_by = n_added + n_subed
                                noise = calculate_noise_level(
                                            self.coadded_map_frames\
                                            [id_for_coadds],
                                            None,
                                            self.point_source_file,
                                            temp_only=self.temp_only,
                                            smaller_region=True,
                                            center_ra=center_ra,
                                            center_dec=center_dec,
                                            map_key_suffix="_Noise",
                                            divisive_factor=divide_map_by)
                                self.log("* (So far, %s and %s maps have been "
                                         "added and subtracted, respectively,",
                                         n_added, n_subed)
                                self.log("*  and the running noise map was "
                                         "divided by %s when calculating "
                                         "the noise.)", divide_map_by)
                            
                            if operation in ["Ignore", "Copy"]:
                                pass
                            elif divide_map_by < 11:
                                pass
                            else:
                                past_obs_ids = \
                                    self.coadded_obs_ids \
                                    [id_for_coadds][src][:-1]
                                past_ops = \
                                    numpy.asarray(
                                        [self.operations_done_to_maps    \
                                         [id_for_coadds][src][str(poid)] \
                                         for poid in past_obs_ids])
                                valid_indices = \
                                    numpy.where(past_ops!=0.0)[0]
                                last_good_obs_id = \
                                    str(past_obs_ids[int(valid_indices[-1])])
                                last_noise  = \
                                    self.noise_from_coadded_maps\
                                    [id_for_coadds][src][last_good_obs_id]
                                if noise < 1.0 * last_noise:
                                    pass
                                else:
                                    self.log("* (The noise calculated "
                                             "this time is larger than the one "
                                             "calculated last time,")
                                    self.log("*  which indicates that "
                                             "this map might be too noisy,")
                                    self.log("*  so the operation done "
                                             "will be reversed.")
                                    self.log("*  For the record,")
                                    self.log("*  the noise is %s uK.arcmin.",
                                             noise/noise_units)
                                    self.log("*  It was %s last time.)",
                                             last_noise/noise_units)
                                    self.log("* (This map will be recorded in "
                                             "the list of bad maps.)")
                                    """record_bad_obs_id(
                                        self.bad_map_list_file,
                                        id_for_coadds, src, oid,
                                        "Made running noise map noisier.")"""
                                    reverse_op = op_dict[operation] * (-1)
                                    self.coadded_map_frames[id_for_coadds] = \
                                        add_map_frames(
                                            self.coadded_map_frames\
                                            [id_for_coadds],
                                            new_frame,
                                            {"T": "T_Noise"},
                                            reverse_op)
                                    self.operations_done_to_maps\
                                        [id_for_coadds][src].pop(str(oid))
                                    self.operations_done_to_maps\
                                        [id_for_coadds][src][str(oid)] = 0.0
                    
                    self.log("* ... the noise level was calculated to be")
                    self.log("* %s uK.arcmin.", noise/noise_units)
                    self.log("* Done.")
                    self.log("")
                    
                    for idr in ids_for_recording:
                        noise_info_to_add = {idr: core.G3MapMapDouble()}
                        noise_info_to_add[idr][src] = core.G3MapDouble()
                        noise_info_to_add[idr][src][str(oid)] = noise
                        if self.calc_noise_from_individ_maps:
                            self.noise_from_individual_maps = \
                                combine_mapmapdoubles(
                                    self.noise_from_individual_maps,
                                    noise_info_to_add)
                        else:
                            self.noise_from_coadded_maps = \
                                combine_mapmapdoubles(
                                    self.noise_from_coadded_maps,
                                    noise_info_to_add)
            
            
            # - Reset the preparations made for
            # - calculating pointing discrepancies and noise levels
            
            if time_to_calculate_pointing_and_noise:
                if self.maps_split_by_scan_direction:
                    if self.calc_noise_from_individ_maps or \
                       self.calc_pointing_discrepancies:
                        for each_id in [id_r_ver, id_l_ver]:
                            self.map_frame_cache[each_id] = None
                        self.map_fr_arrival_counters[id_r_ver] = 0
                        self.map_fr_arrival_counters[id_l_ver] = 0
            
            
            # - Free up some memory
            
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
                mp_fr["CoaddedMapsContained"] = True
                mp_fr["CoaddedMapIDs"] = self.coadded_map_ids[map_id]
                mp_fr["CoaddedObservationIDs"] = self.coadded_obs_ids[map_id]
                mp_fr["ObservationDurations"] = self.observat_durations[map_id]
                
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
                
                self.log("# Map IDs used:")
                self.log("")
                for sf, mids \
                in  mp_fr["CoaddedMapIDs"].items():
                    self.log("- %s", sf)
                    self.log("")
                    self.log(" "*3 + " %s", list(mids))
                    self.log("")
                self.log("\n")
                
                self.log("# Observation durations:")
                self.log("")
                for sub_fld, data in mp_fr["ObservationDurations"].items():
                    self.log(" "*3 + "- %s", sub_fld)
                    self.log("")
                    self.log(" "*6 + "%s",
                             data.keys())
                    self.log(" "*6 + "%s",
                             [t/core.G3Units.min for t in data.values()])
                    self.log("")
                self.log("\n")
                
                
                if self.get_avgs_flagging_stats:
                    mp_fr["FlaggingStatisticsContained"] = True
                    self.log("# Averages related to flagging statistics:")
                    self.log("")

                    prefix = "FlaggingStatisticsAverageNumberOf"
                    for flg_typ in self.avgs_flagging_stats.keys():
                        mp_fr[prefix+flg_typ] = \
                            self.avgs_flagging_stats[flg_typ][map_id]
                    
                        self.log("- %s", flg_typ)
                        self.log("")
                        for sub_fld, data in mp_fr[prefix+flg_typ].items():
                            self.log(" "*3 + "- %s", sub_fld)
                            self.log("")
                            self.log(" "*6 + "%s",
                                     data.keys())
                            self.log(" "*6 + "%s",
                                     [int(avg) for avg in data.values()])
                            self.log("")
                    self.log("\n")
                
                
                if self.calc_mean_pW_to_K:
                    mp_fr["MeansOfTempCalibrationFactorsContained"] = True
                    self.log("# Averages of pW/K conversion factors:")
                    self.log("")
                    
                    prefix = "MeansOfTemperatureCalibrationFactors"
                    uni    = core.G3Units.pW / core.G3Units.K
                    for wafer in self.means_of_temp_cal_factors.keys():
                        mp_fr[prefix+wafer] = \
                            self.means_of_temp_cal_factors[wafer][map_id]
                        
                        self.log("- %s", wafer)
                        self.log("")
                        for sub_fld, data in mp_fr[prefix+wafer].items():
                            self.log(" "*3 + "- %s", sub_fld)
                            self.log("")
                            self.log(" "*6 + "%s",
                                     data.keys())
                            self.log(" "*6 + "%s",
                                     [cf/uni for cf in data.values()])
                            self.log("")
                    self.log("\n")
                
                
                if self.calc_map_stddev_etc:
                    mp_fr["BasicMapFluctuationMetricsContained"] = True
                    self.log("# Basic metrics for fluctuations of map values:")
                    self.log("")
                    
                    prefix = "BasicMapFluctuationMetrics"
                    for mp in self.map_types:
                        self.log("- Map type: %s", mp)
                        self.log("")
                        
                        for flc in self.fluctuation_keys:
                            key_name = prefix + mp + flc
                            mp_fr[key_name] = \
                                self.map_fluctuation_metrics[mp][flc][map_id]
                            
                            self.log(" "*3 + " - Metric: %s", flc)
                            self.log("")
                            if flc == "MapStdDevs":
                                u = core.G3Units.mK
                            elif flc == "MeansOfWeights":
                                u = 1 / (core.G3Units.mK*core.G3Units.mK)
                            else:
                                u = 1.0
                            for sub_fld, data in mp_fr[key_name].items():
                                self.log(" "*6 + " - %s", sub_fld)
                                self.log("")
                                self.log(" "*9 + " %s", data.keys())
                                self.log(" "*9 + " %s",
                                         [str(d/u)[0:6] for d in data.values()])
                                self.log("")
                    self.log("\n")
                
                
                if self.calc_pointing_discrepancies:
                    mp_fr["DeltaRasAndDecsContained"]  = True
                    mp_fr["IntegratedFluxesContained"] = True
                    mp_fr["PointSourceSNRsContained"]  = True
                    self.log("# Pointing discrepancies, fluxes, and SNRs:")
                    self.log("")

                    for rank in ["1", "2", "3"]:
                        self.log("- Point source %s", rank)
                        self.log("")
                        
                        for coordinate in ["Ra", "Dec"]:
                            d_typ    = "delta_" + coordinate.lower()
                            key_name = "Delta"  + coordinate + "s" + \
                                       "FromSources" + rank
                            mp_fr[key_name] = \
                                self.delta_ras_and_decs[rank][d_typ][map_id]
                            
                            self.log(" "*3 + " - %s", d_typ)
                            self.log("")
                            for sub_fld, data in mp_fr[key_name].items():
                                self.log(" "*6 + " - %s", sub_fld)
                                self.log("")
                                self.log(" "*9 + " %s", data.keys())
                                self.log(" "*9 + " %s",
                                         str([str(d/core.G3Units.arcsec)[0:6] \
                                              for d in data.values()]))
                                self.log("")
                        
                        flux_key = "IntegratedFluxesFromSources" + rank
                        mp_fr[flux_key] = self.fluxes[rank][map_id]
                        usys = core.G3Units
                        flu  = usys.mK * usys.arcmin * usys.arcmin
                        
                        self.log(" "*3 + " - Flux")
                        self.log("")
                        for sub_fld, data in mp_fr[flux_key].items():
                            self.log(" "*6 + " - %s", sub_fld)
                            self.log("")
                            self.log(" "*9 + " %s", data.keys())
                            self.log(" "*9 + " %s",
                                     [str(d/flu)[0:6] for d in data.values()])
                            self.log("")
                        
                        snr_key = "PointSourceSNRs" + rank
                        mp_fr[snr_key] = self.snrs[rank][map_id]
                        
                        self.log(" "*3 + " - SNR")
                        self.log("")
                        for sub_fld, data in mp_fr[snr_key].items():
                            self.log(" "*6 + " - %s", sub_fld)
                            self.log("")
                            self.log(" "*9 + " %s", data.keys())
                            self.log(" "*9 + " %s",
                                     [str(d)[0:3] for d in data.values()])
                            self.log("")
                    self.log("\n")
                
                
                if self.calc_xspec_with_coadded_maps:
                    mp_fr["ClsFromCrossSpectraContained"] = True
                    self.log("# Cls from cross spectra:")
                    self.log("")
                    
                    cl_key = "ClsFromCrossSpectraBetweenCoaddsAndIndividuals"
                    mp_fr[cl_key] = self.xspec_cls_with_coadds[map_id]
                    
                    for sub_field, data in mp_fr[cl_key].items():
                        self.log("- %s", sub_field)
                        self.log("")
                        for obs_id in mp_fr["CoaddedObservationIDs"][sub_field]:
                            try:
                                sqrt_cl = data[str(obs_id)]
                            except KeyError:
                                sqrt_cl = numpy.nan
                            sqrt_cl /= (core.G3Units.uK*core.G3Units.arcmin)
                            self.log(" "*3 + "%s %s uK.arcmin", obs_id, sqrt_cl)
                        self.log("")
                    self.log("\n")
                
                
                if self.calc_noise_from_individ_maps or \
                   self.calc_noise_from_coadded_maps:
                    mp_fr["NoiseLevelsContained"] = True
                    ops_k = None
                    self.log("# Noise levels:")
                    self.log("")
                    
                    if self.calc_noise_from_individ_maps:
                        noise_key = "NoiseFromIndividualMaps"
                        mp_fr[noise_key] = \
                            self.noise_from_individual_maps[map_id]
                    else:
                        noise_key = "NoiseFromCoaddedMaps"
                        mp_fr[noise_key] = \
                            self.noise_from_coadded_maps[map_id]
                        if not self.maps_split_by_scan_direction:
                            ops_k = "NoiseCalculationsOperationsDoneToMaps"
                            mp_fr[ops_k] = self.operations_done_to_maps[map_id]
                    
                    for sub_field, data in mp_fr[noise_key].items():
                        self.log("- %s", sub_field)
                        self.log("")
                        for obs_id in mp_fr["CoaddedObservationIDs"][sub_field]:
                            try:
                                noise = data[str(obs_id)]
                            except KeyError:
                                noise = numpy.nan
                            noise /= (core.G3Units.uK*core.G3Units.arcmin)
                            if ops_k is not None:
                                operation = mp_fr[ops_k][sub_field][str(obs_id)]
                            else:
                                operation = ""
                            self.log(" "*3 +"%s, %s uK.arcmin %s",
                                     obs_id, noise, operation)
                        self.log("")
                    self.log("\n")
                
                meta_info_frame = core.G3Frame()
                for k in mp_fr.keys():
                    if k not in ["T", "Q", "U",
                                 "T_Noise", "Q_Noise", "U_Noise",
                                 "Wunpol", "Wpol",
                                 "CoaddedMapsContained"]:
                        meta_info_frame[k] = mp_fr[k]
                meta_info_frame["CoaddedMapsContained"] = False
                meta_info_frames_to_return.append(meta_info_frame)
                
                self.log("# Here is what the frame looks like:")
                self.log("")
                self.log(mp_fr)
                self.log("")
                
                
                self.log("# --------------------------- #")
            self.log("\n")
            
            return meta_info_frames_to_return + \
                   list(self.coadded_map_frames.values()) + \
                   [frame]




class AppendDirectionsToMapIDs(object):
    
    def __init__(self, map_id):
        
        self.map_id = map_id.replace("Left", "")
        self.n_mp_l = 0
        self.n_mp_r = 0
    
    def __call__(self, frame):
        
        if "CoaddedObservationIDs" in frame.keys():
            n_maps = len(frame["CoaddedObservationIDs"].values()[0])
            if "Left" in frame["Id"]:
                self.n_mp_l = n_maps
            elif "Right" in frame["Id"]:
                self.n_mp_r = n_maps
        
        if (frame.type == core.G3FrameType.Map) and \
           (frame["Id"] in ["90GHz", "150GHz", "220GHz"]):
            old_id = frame.pop("Id")
            if self.n_mp_l > self.n_mp_r:
                frame["Id"] = "Right" + old_id
                self.n_mp_r += 1
            else:
                frame["Id"] = "Left" + old_id
                self.n_mp_l += 1




def print_recorded_information(
        filepath, info_tuples, print_frame=True, output_path=None):
    
    iterator = core.G3File(filepath)
    
    while True:
        new_frame = iterator.next()
        if "CoaddedObservationIDs" in new_frame.keys():
            database = new_frame
            break
    
    contents  = "\n"
    if print_frame:
        contents += "# -------------------------- #\n"
        contents += "#  The list of keys stored:  #\n"
        contents += "# -------------------------- #\n"
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
# Define functions needed to start the pipeline
# ------------------------------------------------------------------------------


def is_a_good_obs_id(g3_file, min_obs_id, max_obs_id, bad_obs_ids=[]):
    
    meaningful_part_of_path = g3_file.split('/')[-1].split('.')[0]
    try:
        obs_id = int(meaningful_part_of_path.split('_')[0])
    except ValueError:
        # Probably a g3 file that contains coadded maps is encountered
        # if the name cannot be converted to an integer.
        return True
    
    if (obs_id >= min_obs_id) and \
       (obs_id <= max_obs_id) and \
       (obs_id not in bad_obs_ids):
        return True
    else:
        return False



def run(input_files=[], output_file='./coadded_maps.g3', map_ids=["90GHz"],
        temperature_maps_only=False,
        trick_pipeline_into_getting_left_right_maps=False,
        maps_split_by_scan_direction=False,
        combine_left_right=False,
        combine_different_wafers=False,
        subtract_existing_maps=False, 
        collect_averages_from_flagging_statistics=False,
        calculate_pW_to_K_conversion_factors=False,
        calculate_map_rms_and_weight_stat=False,
        calculate_noise_from_individual_maps=False,
        calculate_noise_from_coadded_maps=False,
        calculate_cross_spectrum_with_coadded_maps=False,
        point_source_file='',
        bad_map_list_file='',
        calculate_pointing_discrepancies=False,
        sources=["ra0hdec-44.75", "ra0hdec-52.25",
                 "ra0hdec-59.75", "ra0hdec-67.25"],
        min_obs_id=0,max_obs_id=99999999,
        bad_obs_ids=[], min_file_size=0.01,
        logger_name="", log_file=None):
    
    
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
    
    
    # - Check consistencies among some options that are supposed to be used

    """if calculate_noise_from_coadded_maps:
        if combine_left_right:
            log("")
            log("If noise is to be calculated from running coadds,")
            log("then the combine_left_right option needs to be False")
            log("so that the coadded left-going and right-going maps")
            log("can be used later as well!")
            log("")
            sys.exit()"""

    if calculate_noise_from_individual_maps and \
       calculate_noise_from_coadded_maps:
        log("")
        log("Currently, calculating noise from both individual maps")
        log("and coadded maps is not supported...")
        log("")
        return 0
    
    if calculate_noise_from_coadded_maps and \
       subtract_existing_maps:
        log("")
        log("If noise is to be calculated from running coadded maps,")
        log("then subtracting some maps that were previously used in")
        log("the coadded maps is not allowed because it seems that")
        log("this would make the results of previous noise level")
        log("calculations not meaningful.")
        log("")
        return 0
    
    if calculate_cross_spectrum_with_coadded_maps and \
       maps_split_by_scan_direction:
        log("")
        log("Currently, calculating cross spectra between")
        log("individual observations' maps and running coadded maps")
        log("in the case when the maps are split by scan direction")
        log("is not supported...")
        log("")
        return 0


    all_good_g3_files = []
    for g3_file in input_files:
        if os.path.isfile(g3_file) and \
           is_a_good_obs_id(g3_file, min_obs_id, max_obs_id, bad_obs_ids) and \
           (os.path.getsize(g3_file) > min_file_size*2**30):
            all_good_g3_files.append(g3_file)


    if (len(all_good_g3_files) <= 1):
        log("")
        log("No applicable input files, so nothing to do!")
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

    log("- Observation IDs that were originally on the list")
    log("  but then excluded:")
    log("")
    log(bad_obs_ids)
    log("")


    # - Construct a pipeline that makes coadded maps
        
    pipeline = core.G3Pipeline()

    pipeline.Add(core.G3Reader,
                 filename=all_good_g3_files)
    
    if trick_pipeline_into_getting_left_right_maps:
        pipeline.Add(AppendDirectionsToMapIDs,
                     map_id=map_ids[0])
    
    pipeline.Add(lambda frame: log(frame))
    
    pipeline.Add(CoaddMapsAndDoSomeMapAnalysis,
                 map_ids=map_ids,
                 map_sources=sources,
                 temperature_only=temperature_maps_only,
                 maps_split_by_scan_direction=\
                     maps_split_by_scan_direction,
                 combine_left_right=\
                     combine_left_right,
                 combine_different_wafers=\
                     combine_different_wafers,
                 allow_subtraction=\
                     subtract_existing_maps,
                 collect_averages_from_flagging_stats=\
                     collect_averages_from_flagging_statistics,
                 calculate_pW_to_K_conversion_factors=\
                     calculate_pW_to_K_conversion_factors,
                 calculate_map_rms_and_weight_stat=\
                     calculate_map_rms_and_weight_stat,
                 calculate_noise_from_individual_maps=\
                     calculate_noise_from_individual_maps,
                 calculate_noise_from_coadded_maps=\
                     calculate_noise_from_coadded_maps,
                 calculate_xspec_with_coadded_maps=\
                     calculate_cross_spectrum_with_coadded_maps,
                 point_source_file=\
                     point_source_file,
                 bad_map_list_file=\
                     bad_map_list_file,
                 calculate_pointing_discrepancies=\
                     calculate_pointing_discrepancies,
                 logging_function=log)

    pipeline.Add(lambda frame: "CoaddedMapsContained" in frame)

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
                             "analysis on the maps such as calculating noise.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input_files",
                        type=str, action="store", nargs="+",
                        help="The paths to the g3 files "
                             "that contain map frames. ")

    parser.add_argument("-o", "--output_file",
                        type=str, action="store", default="./coadded_maps.g3",
                        help="The location and name of the output g3 file "
                             "that will contain the coadded maps.")

    parser.add_argument("-i", "--map_ids",
                        type=str, action="store", nargs="+", default=["90GHz"],
                        help="The relevant map IDs for which corresponding "
                             "map data will be added. ")

    parser.add_argument("-t", "--temperature_maps_only",
                        action="store_true", default=False,
                        help="Whether the maps stored in the data frames "
                             "contain temperature-only maps or not.")
    
    parser.add_argument("-P", "--trick_pipeline_into_getting_left_right_maps",
                        action="store_true", default=False,
                        help="Whether to trick the pipeline module that "
                             "co-adds maps and performs analyses into "
                             "receiving maps whose IDs imply that "
                             "they are made from only one scan direction. "
                             "This option is useful in the case where "
                             "we want to calculate noise levels from "
                             "running coadded maps but the input files "
                             "do not contain maps that are split by "
                             "scan direction.")

    parser.add_argument("-r", "--maps_split_by_scan_direction",
                        action="store_true", default=False,
                        help="Whether the maps stored in the data frames are "
                             "split by the scan direction.")

    parser.add_argument("-c", "--combine_left_right",
                        action="store_true", default=False,
                        help="Whether to add maps made from left-going scans "
                             "and those made from right-going scans together. "
                             "If True, a map whose 'Id' is 'Left'+map_id and "
                             "one whose 'Id' is 'Right'+map_id will be added, "
                             "where map_id is an ID specified in the argument "
                             "map_ids above.")

    parser.add_argument("-C", "--combine_different_wafers",
                        action="store_true", default=False,
                        help="Whether to add maps made from different wafers "
                             "together. If True, maps whose 'Id's look like "
                             "map_id+a_wafer_name will be added, where"
                             "map_id is an ID specified in the argument "
                             "map_ids above.")

    parser.add_argument("-u", "--subtract_existing_maps",
                        action="store_true", default=False,
                        help="Whether to subtract certain observation's data "
                             "contained in an input file from the coadded maps "
                             "if the data already exist in the coadded maps.")

    parser.add_argument("-a", "--collect_averages_from_flagging_statistics",
                        action="store_true", default=False,
                        help="Whether to calculate average number of "
                             "bolometers that were not flagged "
                             "over all scans of an observation, "
                             "average number of bolometers that "
                             "were removed, and so on.")
    
    parser.add_argument("-k", "--calculate_pW_to_K_conversion_factors",
                        action="store_true", default=False,
                        help="Whether to calculate the mean value of the "
                             "pW/K absolute calibration factors for each band "
                             "of each wafer.")
    
    parser.add_argument("-w", "--calculate_map_rms_and_weight_stat",
                        action="store_true", default=False,
                        help="Whether to calculate an RMS value of the map, "
                             "a mean weight, and a product of these for "
                             "individual maps and coadded maps.")

    parser.add_argument("-n", "--calculate_noise_from_individual_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels from "
                             "each map that goes into the coadded maps.")

    parser.add_argument("-N", "--calculate_noise_from_coadded_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels for the "
                             "running coadded maps. In other words, the noise "
                             "will be calculated from the coadded maps "
                             "every time a new map is added.")
    
    parser.add_argument("-x", "--calculate_cross_spectrum_with_coadded_maps",
                        action="store_true", default=False,
                        help="Whether to calculate the cross spectrum between "
                             "each observation's map and the running "
                             "coadded maps.")

    parser.add_argument("-M", "--point_source_file",
                        type=str, action="store", default="",
                        help="Path to a point source list, which will be used "
                             "when making a mask for calcuting Cls.")

    parser.add_argument("-p", "--calculate_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to calculate the difference between the "
                             "true positions of some point sources and the "
                             "measured positions")

    parser.add_argument("-S", "--sources",
                        type=str, action="store", nargs="+",
                        default=["ra0hdec-44.75", "ra0hdec-52.25",
                                 "ra0hdec-59.75", "ra0hdec-67.25"],
                        help="The sub-field(s) that will be included "
                             "in the coadded maps.")

    parser.add_argument("-s", "--min_obs_id",
                        type=int, action="store", default=0,
                        help="The smallest observation ID that will be "
                             "considered to be used in the coadded map. "
                             "The default is 00000000.")

    parser.add_argument("-l", "--max_obs_id",
                        type=int, action="store", default=99999999,
                        help="The largest observation ID that will be "
                             "considered to be used in the coadded map. "
                             "The default is 99999999.")

    parser.add_argument("-b", "--bad_obs_ids",
                        type=int, action="store", nargs="+", default=[],
                        help="The observation IDs that will be excluded from "
                             "making the coadded map. The script has some "
                             "simple criteria to decide whether to include "
                             "certain observations, but one can manually "
                             "specify bad observations, too.")
    
    parser.add_argument("-B", "--bad_map_list_file",
                        type=str, action="store", default="",
                        help="A text file in which maps that are found to be "
                             "bad during the coadding processes are recorded.")

    parser.add_argument("-m", "--min_file_size",
                        type=float, action="store", default=0.01,
                        help="The minimum size (in GB) a g3 file needs "
                             "to have for it to be considered as a good file. "
                             "This is to reduce the occurrence of a situation "
                             "where the script crashes due to reading a "
                             "problematic file.")
    
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

