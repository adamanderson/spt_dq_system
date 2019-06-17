## *******  Always need to import modules!  ******* ##

from    spt3g              import  core
from    spt3g              import  mapmaker
from    spt3g              import  util
from    spt3g              import  coordinateutils
from    spt3g              import  mapspectra
from    spt3g.mapspectra   import  map_analysis
from    scipy              import  signal
import  os
import  gc
import  sys
import  glob
import  numpy
import  argparse
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
    
    """
    ## Ignore this part for now
    too_big_map_vals = {"T": {}}
    for map_id in arguments.map_ids:
        if "90GHz" in map_id:
            too_big_map_vals["T"][map_id] = numpy.inf
        elif "150GHz" in map_id:
            too_big_map_vals["T"][map_id] = numpy.inf
        elif "220GHz" in map_id:
            too_big_map_vals["T"][map_id] = numpy.inf
    pctl_10 = numpy.nanpercentile(numpy.asarray(map_weight_divided), 90)
    pctl_90 = numpy.nanpercentile(numpy.asarray(map_weight_divided), 10)
    larger_abs = numpy.max([numpy.abs(pctl_10), numpy.abs(pctl_90)])
    if larger_abs/core.G3Units.mK > dict_thresh[tqu_type][id_for_thresh]:
        print()
        print("The fluctuations seem too large!")
        print("10th pctl. value in the map:", pctl_10, "mK")
        print("90th pctl. value in the map:", pctl_90, "mK")
        return False
    """
    
    return True



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

    # print(n_bolos_from_each_scan)
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
                # print(reason, number)
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
        pw_k  = calframe["CalibratorResponse"][bolo]
        pw_k *= calframe[flux_calibration_key][bolo]
        pw_k *= calframe[integral_flux_key][bolo]
        pw_k *= calframe[sky_transmission_key][band]
        pw_k /= flx[source][band]
        values_by_band_and_wafer[band][wafer].append(pw_k)
    
    for band in values_by_band_and_wafer.keys():
        for wafer in wafers:
            vals = values_by_band_and_wafer[band][wafer]
            if True not in numpy.isfinite(vals):
                values_by_band_and_wafer[band][wafer] = numpy.nan
            else:
                values_by_band_and_wafer[band][wafer] = numpy.nanmean(vals)
    
    return values_by_band_and_wafer




def calculate_pointing_discrepancies(fr_one, fr_two, sub_field, temp_only=True):

    u = core.G3Units

    if temp_only:
        if fr_two != None:
            total_weight = fr_one["Wunpol"].TT + fr_two["Wunpol"].TT
            sum_map      = fr_one["T"] + fr_two["T"]
            map_to_use   = sum_map
        else:
            map_to_use = fr_one["T"] / fr_one["Wunpol"].TT
    
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
        
        if len(numpy.where(numpy.isnan(mini_map))[0]) > 0:
            discrep_dict[point_source_number]["delta_ra"]  = numpy.nan
            discrep_dict[point_source_number]["delta_dec"] = numpy.nan
            continue
        
        fit_parameters = util.fitting.fit_gaussian2d(mini_map)
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
        
    return discrep_dict, flux_dict




def create_new_map_frame_with_smaller_map_region(
        map_frame, temperature_only,
        center_ra, center_dec,
        look_for_noise_map=False, divisive_factor=1.0):

    new_mp_fr = core.G3Frame(core.G3FrameType.Map)

    if temperature_only:
        if look_for_noise_map and \
           ("T_Noise" in map_frame.keys()):
            original_sky_map = map_frame["T_Noise"] / divisive_factor
        else:
            original_sky_map = map_frame["T"]
        original_wght_map = map_frame["Wunpol"].TT

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
        smaller_wght_map = coordinateutils.FlatSkyMap(
                x_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                y_len=int(5.0*core.G3Units.deg/original_sky_map.res),
                res=original_sky_map.res,
                proj=original_sky_map.proj,
                alpha_center=center_ra,
                delta_center=center_dec,
                pol_type=original_sky_map.pol_type,
                coord_ref=original_sky_map.coord_ref)

        coordinateutils.reproj_map(original_sky_map,  smaller_sky_map,  1)
        coordinateutils.reproj_map(original_wght_map, smaller_wght_map, 1)

        new_mp_fr["T"] = smaller_sky_map
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
                frame, temp_only, center_ra, center_dec, None)
    
    if (False in numpy.isfinite(numpy.asarray(mp_fr["T"])))         or \
       (False in numpy.isfinite(numpy.asarray(mp_fr["Wunpol"].TT))) or \
       (0.0 in numpy.asarray(mp_fr["Wunpol"].TT)):
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
        look_for_noise_map=False, divisive_factor=None):

    if fr_two != None:
        mp_fr = map_analysis.subtract_two_maps(
                    fr_one, fr_two, divide_by_two=True, in_place=False)
    else:
        mp_fr = fr_one
    
    if smaller_region:
        mp_fr = create_new_map_frame_with_smaller_map_region(
                    mp_fr, temp_only,
                    center_ra, center_dec,
                    look_for_noise_map=look_for_noise_map,
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
                    """old_dict[map_id][source] = \
                    type(new_ids)(sorted(
                    set(old_dict[map_id][source]) | set(new_ids)))"""
                    # ** The order of the IDs may be important.
    
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
                 calculate_noise_from_individual_maps=False,
                 calculate_noise_from_coadded_maps=False,
                 calculate_xspec_with_coadded_maps=False,
                 point_source_file=None,
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

        
        # - Initialize variables related to storing coadded maps
        
        self.map_sources = map_sources
        self.temp_only   = temperature_only
        
        self.coadded_map_frames = {map_id: core.G3Frame(core.G3FrameType.Map) \
                                   for map_id in self.map_ids}
        self.coadded_obs_ids    = {map_id: core.G3MapVectorInt()    \
                                   for map_id in self.map_ids}
        self.coadded_map_ids    = {map_id: core.G3MapVectorString() \
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
            self.log("")
            self.log("-----------------------------------------------------")
            self.log("Found an observation frame!")
            self.log("(Probably a new g3 file has arrived to the pipeline.)")
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
                    self.log("* because its map ID is not one of the desired IDs.")
                    self.log("\n")
                    return []
                else:
                    obs_ids_from_this_frame = \
                        {id_for_coadds: frame["CoaddedObservationIDs"]}
                    map_ids_from_this_frame = \
                        {id_for_coadds: frame["CoaddedMapIDs"]}
            else:
                if frame["Id"] not in self.id_mapping.keys():
                    self.log("")
                    self.log("* Skipping the map frame above")
                    self.log("* because its map ID is not one of the desired IDs.")
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
                    self.log("* because the observation type is not of interest.")
                    self.log("\n")
                    return []
                else:
                    id_for_coadds  = self.id_mapping[frame["Id"]]
                    oid = self.obs_info["ObservationID"]
                    src = self.obs_info["SourceName"]
                    dmid = str(oid) + frame["Id"]   # ** more detailed map id
                    obs_ids_from_this_frame = \
                        {id_for_coadds: {src: core.G3VectorInt([oid])}}
                    map_ids_from_this_frame = \
                        {id_for_coadds: {src: core.G3VectorString([dmid])}}
            
            
            # - Perform a few basic checks on the frame
            # - to furthur figure out whether it can go into the coadds
            
            if not map_seems_fine(frame, self.log):
                self.log("")
                self.log("* Well, the map doesn't look good,")
                self.log("* so, this one will be skipped...")
                self.log("\n")
                return []
            
            self.coadded_obs_ids, this_frame_already_added = \
                combine_ids(self.coadded_obs_ids,
                            obs_ids_from_this_frame)
            self.coadded_map_ids, this_frame_already_added = \
                combine_ids(self.coadded_map_ids,
                            map_ids_from_this_frame)
            
            if this_frame_already_added and \
               (not self.allow_subtraction):
                self.log("")
                self.log("* Well, it looks that this map frame")
                self.log("* has already been added (at least partially),")
                self.log("* so, it will be skipped...")
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
                    xspec_sqrt_cls_this_fr = {id_for_coadds: frame[cl_key]}
                elif len(self.coadded_obs_ids[id_for_coadds][src]) <= 1:
                    self.log("")
                    self.log("* Currently, the maps in the cache")
                    self.log("* (or maps to be added to the cache) are the same as")
                    self.log("* the maps in this frame, so the cross spectrum")
                    self.log("* calculation will be skipped for this frame.")
                    xspec_sqrt_cls = core.G3MapMapDouble()
                    xspec_sqrt_cls[src] = core.G3MapDouble()
                    xspec_sqrt_cls[src][str(oid)] = numpy.nan
                    xspec_sqrt_cls_this_fr = {id_for_coadds: xspec_sqrt_cls}
                else:
                    self.log("")
                    self.log("* Calculating avg.(sqrt(Cl)) of the cross spectrum")
                    self.log("* between the coadded maps and the map from")
                    self.log("* this observation ...")
                    center_ra  = 0.0 * core.G3Units.deg
                    center_dec = float(src[-6:]) * core.G3Units.deg
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
                        
            if len(self.coadded_map_frames[id_for_coadds].keys()) == 0:
                self.log("")
                self.log("* A map frame with ID %s is here!", frame["Id"])
                self.log("* Both the sky map and weight map are added to the ")
                self.log("* empty coadded maps as %s ...", id_for_coadds)
                if self.temp_only:
                    self.coadded_map_frames \
                        [id_for_coadds]["T"] = frame["T"]
                    self.coadded_map_frames \
                        [id_for_coadds]["Wunpol"] = frame["Wunpol"]
                self.log("* Done.")
                self.log("")
            elif not subtract_maps_in_this_frame:
                self.log("")
                self.log("* Adding another set of sky map and weight map "
                         "(%s)", frame["Id"])
                self.log("* to the coadded maps for %s ...", id_for_coadds)
                if self.temp_only:
                    existing_t_map = \
                        self.coadded_map_frames[id_for_coadds].pop("T")
                    existing_w_map = \
                        self.coadded_map_frames[id_for_coadds].pop("Wunpol")
                    self.coadded_map_frames[id_for_coadds]["T"] = \
                        existing_t_map + frame["T"]
                    self.coadded_map_frames[id_for_coadds]["Wunpol"] = \
                        existing_w_map + frame["Wunpol"]
                self.log("* Done.")
                self.log("")
            else:
                self.coadded_obs_ids = \
                    remove_ids(self.coadded_obs_ids, obs_ids_from_this_frame)
                self.coadded_map_ids = \
                    remove_ids(self.coadded_map_ids, map_ids_from_this_frame)
                self.log("")
                self.log("* Removing a set of sky map and weight map (%s)", frame["Id"])
                self.log("* from the coadded maps for %s ...", id_for_coadds)
                if self.temp_only:
                    existing_t_map = \
                        self.coadded_map_frames[id_for_coadds].pop("T")
                    self.coadded_map_frames[id_for_coadds]["T"] = \
                        existing_t_map + (frame["T"] * (-1))
                    self.coadded_map_frames[id_for_coadds]["Wunpol"].TT += \
                        (frame["Wunpol"].TT * (-1))
                self.log("* Done.")
                self.log("")
            
            
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
                
                if self.calc_noise_from_individ_maps:
                    self.noise_from_individual_maps = \
                        remove_partial_mapmapdouble(
                            self.noise_from_individual_maps,
                            obs_ids_from_this_frame)
                """elif self.calc_noise_from_coadded_maps:
                    self.noise_from_coadded_maps = \
                        remove_partial_mapmapdouble(
                            self.noise_from_coadded_maps,
                            obs_ids_from_this_frame)
                    if not self.maps_split_by_scan_direction:
                        ops_done = self.operations_done_to_maps[src][oid]
                        inverse_ops = -1.0 * ops_done
                        if self.temp_only:
                            if inverse_ops == 0.0:
                                pass
                            else:
                                map_from_this_obs  = \
                                    frame["T"] / frame["Wunpol"].TT
                                map_from_this_obs.is_weighted = False
                                existing_noise_map = \
                                    self.coadded_map_frames\
                                    [id_for_coadds].pop("T_Noise")
                                self.coadded_map_frames["T_Noise"] = \
                                    existing_noise_map + \
                                    (map_from_this_obs * inverse_ops)
                        
                        self.operations_done_to_maps = \
                            remove_partial_mapmapdouble(
                                self.operations_done_to_maps,
                                obs_ids_from_this_frame)"""
                
                time_to_calculate_pointing_and_noise = False
            
            elif self.calc_noise_from_coadded_maps or \
                 self.calc_noise_from_individ_maps or \
                 self.calc_pointing_discrepancies:
                
                if frame_has_coadds:
                    time_to_calculate_pointing_and_noise = True
                    # ** The assignment above is just a dummy one
                
                elif self.maps_split_by_scan_direction:
                    self.map_frame_arrival_counters[frame["Id"]] = 1
                    if self.calc_noise_from_individ_maps or \
                       self.calc_pointing_discrepancies:
                        self.map_frame_cache[frame["Id"]] = frame
                    
                    time_to_calculate_pointing_and_noise = False
                    for di_mid in self.direction_independent_map_ids:
                        if "SomeDirection" in di_mid:
                            id_l_ver = map_id.replace("SomeDirection", "Left")
                            id_r_ver = map_id.replace("SomeDirection", "Right")
                        else:
                            id_l_ver = "Left"  + map_id
                            id_r_ver = "Right" + map_id
                        if (self.map_fr_arrival_counters[id_l_ver] == 1) and \
                           (self.map_fr_arrival_counters[id_r_ver] == 1):
                            time_to_calculate_pointing_and_noise = True
                            break
                else:
                    time_to_calculate_pointing_and_noise = True
                
                if time_to_calculate_pointing_and_noise:
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
                    self.log("* Done.")
                    self.log("")
                else:
                    if self.maps_split_by_scan_direction:
                        self.log("")
                        self.log("* Pointing discrepancy and flux calculations")
                        self.log("* are about to start because both")
                        self.log("* %s and %s", id_l_ver, id_r_ver)
                        self.log("* of obs. %s", oid)
                        self.log("* have passed the pipeline ...")
                        
                        pt_offsets_to_be_reorganized, \
                        fluxes_to_be_reorganized =    \
                            calculate_pointing_discrepancies(
                                self.map_frame_cache[id_r_ver],
                                self.map_frame_cache[id_l_ver],
                                src, temp_only=True)
                    else:
                        self.log("")
                        self.log("* Pointing discrepancy and flux calculations")
                        self.log("* are about to start for")
                        self.log("* obs. %s %s ...", oid, id_for_coadds)
                        
                        pt_offsets_to_be_reorganized, \
                        fluxes_to_be_reorganized =    \
                            calculate_pointing_discrepancies(
                                frame, None, src, temp_only=True)
                        
                    for rank, diffs in pt_offsets_to_be_reorganized.items():
                        for delta_coord, diff in diffs.items():
                            self.log("* Source %s " + delta_coord.ljust(9) + " %s arcseconds",
                                     rank, diff/core.G3Units.arcsec)
                        flux  = fluxes_to_be_reorganized[rank]
                        usys  = core.G3Units
                        flux /= (usys.mK * usys.arcmin * usys.arcmin)
                        self.log("* Source %s " + "flux".ljust(9) + " %s mK.arcmin.arcmin",
                                 rank, flux)
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
                            
                            self.fluxes[rank] = \
                                combine_mapmapdoubles(
                                    self.fluxes[rank],
                                    fluxes_from_this_fr[rank])
            
            
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
                        if not self.maps_split_by_scan_direction:
                            if self.temp_only:
                                if "T_Noise" in self.coadded_map_frames[id_for_coadds].keys():
                                    existing_noise_map  = self.coadded_map_frames[id_for_coadds].pop("T_Noise")
                                    existing_noise_map += frame["T_Noise"]
                                    self.coadded_map_frames[id_for_coadds]["T_Noise"] = existing_noise_map
                                else:
                                    self.coadded_map_frames[id_for_coadds]["T_Noise"] = frame["T_Noise"]
                            ops_key = "NoiseCalculationsOperationsDoneToMaps"
                            self.operations_done_to_maps = \
                                combine_mapmapdoubles(
                                    self.operations_done_to_maps,
                                    {id_for_coadds: frame[ops_key]})
                    
                    self.log("* Done.")
                    self.log("")
                else:
                    center_ra  = 0.0 * core.G3Units.deg
                    center_dec = float(src[-6:]) * core.G3Units.deg
                    if self.maps_split_by_scan_direction:
                        self.log("")
                        self.log("* Noise calculation is about to start")
                        self.log("* b/c both %s and %s", id_l_ver, id_r_ver)
                        self.log("* of obs. %s", oid)
                        self.log("* have passed the pipeline...")
                        
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
                            
                            operation = decide_operation_to_do_with_new_map(
                                            frame,
                                            self.operations_done_to_maps\
                                            [id_for_coadds][src],
                                            self.coadded_obs_ids\
                                            [id_for_coadds][src],
                                            self.temp_only,
                                            center_ra, center_dec)
                            op_dict = {"Copy"    :  1,
                                       "Add"     :  1,
                                       "Subtract": -1,
                                       "Ignore"  :  0}
                            
                            self.operations_done_to_maps\
                                [id_for_coadds][src][str(oid)] = op_dict[operation]                            
                            self.log("* (The operation to be performed on the map ")
                            self.log("*  from this frame is %s.)", operation)
                            
                            if operation == "Copy":
                                if self.temp_only:
                                    self.coadded_map_frames\
                                    [id_for_coadds]["T_Noise"] = \
                                        frame["T"] / frame["Wunpol"].TT
                                    self.coadded_map_frames\
                                    [id_for_coadds]["T_Noise"].\
                                        is_weighted = False
                            elif operation == "Ignore":
                                pass
                            else:
                                if operation == "Add":
                                    multiplicative_factor = 1.0
                                elif operation == "Subtract":
                                    multiplicative_factor = -1.0
                                
                                if self.temp_only:
                                    new_mp = frame["T"] / frame["Wunpol"].TT
                                    new_mp.is_weighted = False
                                    existing_map = \
                                        self.coadded_map_frames\
                                        [id_for_coadds].pop("T_Noise")
                                    self.coadded_map_frames\
                                    [id_for_coadds]["T_Noise"] = \
                                        existing_map + \
                                        (multiplicative_factor * new_mp)
                            
                            n_added = list(self.operations_done_to_maps\
                                         [id_for_coadds][src].values()).\
                                      count(1.0)
                            n_subed = list(self.operations_done_to_maps\
                                         [id_for_coadds][src].values()).\
                                      count(-1.0)
                            
                            if operation == "Ignore":
                                self.log("* (Since this map was bad in some way,")
                                self.log("*  and it was not added nor subtracted,")
                                self.log("*  the noise calculation will not occur.)")
                                noise = numpy.nan
                                pass
                            elif (n_added == n_subed) and (n_added > 0):
                                self.log("* (Since %s pairs of " +\
                                         "difference map have been combined", n_added)
                                self.log("*  into the coadded noise maps,")
                                self.log("*  now is a good time to calculate " +\
                                         "the noise level.)")
                                divide_map_by = n_added * 2
                                noise = calculate_noise_level(
                                            self.coadded_map_frames\
                                            [id_for_coadds],
                                            None,
                                            self.point_source_file,
                                            temp_only=self.temp_only,
                                            smaller_region=True,
                                            center_ra=center_ra,
                                            center_dec=center_dec,
                                            look_for_noise_map=True,
                                            divisive_factor=divide_map_by)
                            else:
                                self.log("* (Since the current " +\
                                         "coadded 'noise' maps don't")
                                self.log("*  contain a non-zero integer number " +\
                                         "of pair(s) of difference map,")
                                self.log("*  now is not a good time to " +\
                                         "calculate the noise level.)")
                                noise = numpy.nan
                                                            
                    n = noise/(core.G3Units.uK*core.G3Units.arcmin)
                    self.log("* ... the noise level was calculated to be")
                    self.log("* %s uK.arcmin.", n)
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
                    
                    self.map_frame_arrival_counters[id_r_ver] = 0
                    self.map_frame_arrival_counters[id_l_ver] = 0
            
            
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
                            self.log(" "*6 + "%s", data.keys())
                            self.log(" "*6 + "%s", [int(avg) for avg in data.values()])
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
                            self.log(" "*6 + "%s", data.keys())
                            self.log(" "*6 + "%s", [cf/uni for cf in data.values()])
                            self.log("")
                    self.log("\n")
                
                
                if self.calc_pointing_discrepancies:
                    mp_fr["DeltaRasAndDecsContained"]  = True
                    mp_fr["IntegratedFluxesContained"] = True
                    self.log("# Pointing discrepancies and fluxes:")
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
                            
                            self.log(" "*3+" - %s", d_typ)
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
                            sqrt_cl  = data[str(obs_id)]
                            sqrt_cl /= (core.G3Units.uK*core.G3Units.arcmin)
                            self.log(" "*3 + "%s %s uK.arcmin", obs_id, sqrt_cl)
                        self.log("")
                    self.log("\n")
                
                
                if self.calc_noise_from_individ_maps or \
                   self.calc_noise_from_coadded_maps:
                    mp_fr["NoiseLevelsContained"] = True
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
                            noise = data[str(obs_id)]
                            noise /= (core.G3Units.uK*core.G3Units.arcmin)
                            self.log(" "*3 +" %s, %s uK.arcmin", obs_id, noise)
                        self.log("")
                        if (not self.calc_noise_from_individ_maps) and \
                           (not self.maps_split_by_scan_direction):
                            self.log(" "*3 + "Operations done to maps" +\
                                "during noise calculations: %s",
                                [mp_fr[ops_k][sub_field][str(oid)] for oid in \
                                 mp_fr["CoaddedObservationIDs"][sub_field]])
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
        temperature_maps_only=False, maps_split_by_scan_direction=False,
        combine_left_right=False, combine_different_wafers=False,
        subtract_existing_maps=False, 
        collect_averages_from_flagging_statistics=False,
        calculate_pW_to_K_conversion_factors=False,
        calculate_noise_from_individual_maps=False,
        calculate_noise_from_coadded_maps=False,
        calculate_cross_spectrum_with_coadded_maps=False,
        point_source_file='',
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

    if calculate_noise_from_coadded_maps:
        """if len(sources) > 1:
            print()
            print("If noise is to be calculated from running coadds,")
            print("then only one observation type can be specified!")
            print()
            sys.exit()"""
        if combine_left_right:
            log("")
            log("If noise is to be calculated from running coadds,")
            log("then the combine_left_right option needs to be False")
            log("so that the coadded left-going and right-going maps")
            log("can be used later as well!")
            log("")
            sys.exit()

    if calculate_noise_from_individual_maps and \
       calculate_noise_from_coadded_maps:
        log("")
        log("Currently, calculating noise from both individual maps")
        log("and coadded maps is not supported...")
        log("")
        sys.exit()
    
    if calculate_noise_from_coadded_maps and \
       subtract_existing_maps:
        log("")
        log("If noise is to be calculated from running coadded maps,")
        log("then subtracting some maps that were previously used in")
        log("the coadded maps is not allowed because it seems that")
        log("this would make the results of previous noise level")
        log("calculations not meaningful.")
        log("")
        sys.exit()
    
    if calculate_cross_spectrum_with_coadded_maps and \
       maps_split_by_scan_direction:
        log("")
        log("Currently, calculating cross spectra between")
        log("individual observations' maps and running coadded maps")
        log("in the case when the maps are split by scan direction")
        log("is not supported...")
        log("")


    all_good_g3_files = []
    for g3_file in input_files:
        if os.path.isfile(g3_file) and \
           is_a_good_obs_id(g3_file, min_obs_id, max_obs_id, bad_obs_ids) and \
           (os.path.getsize(g3_file) > min_file_size*2**30):
            all_good_g3_files.append(g3_file)


    if (len(all_good_g3_files) <= 1):
        log("")
        log("No applicable input files, so nothing to do!")
        log("\n\n")
        return


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
                 calculate_noise_from_individual_maps=\
                     calculate_noise_from_individual_maps,
                 calculate_noise_from_coadded_maps=\
                     calculate_noise_from_coadded_maps,
                 calculate_xspec_with_coadded_maps=\
                     calculate_cross_spectrum_with_coadded_maps,
                 point_source_file=\
                     point_source_file,
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
    
    log("\n")
    log("# ============================ #")
    log("\n\n\n")

# ==============================================================================





# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                 description="This script makes coadded maps of the "+\
                             "1500 sq. deg. field and can also perform some "+\
                             "analysis on the maps such as calculating noise.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input_files",
                        type=str, action="store", nargs="+",
                        help="The paths to the g3 files "+\
                             "that contain map frames. ")

    parser.add_argument("-o", "--output_file",
                        type=str, action="store", default="./coadded_maps.g3",
                        help="The location and name of the output g3 file "+\
                             "that will contain the coadded maps.")

    parser.add_argument("-i", "--map_ids",
                        type=str, action="store", nargs="+", default=["90GHz"],
                        help="The relevant map IDs for which corresponding "+\
                             "map data will be added. ")

    parser.add_argument("-t", "--temperature_maps_only",
                        action="store_true", default=False,
                        help="Whether the maps stored in the data frames contain "+\
                             "temperature-only maps or not.")

    parser.add_argument("-r", "--maps_split_by_scan_direction",
                        action="store_true", default=False,
                        help="Whether the maps stored in the data frames are "+\
                             "split by the scan direction.")

    parser.add_argument("-c", "--combine_left_right",
                        action="store_true", default=False,
                        help="Whether to add maps made from left-going scans and "+\
                             "those made from right-going scans together. "+\
                             "If True, a map whose 'Id' is 'Left'+map_id and "+\
                             "one whose 'Id' is 'Right'+map_id will be added, "+\
                             "where map_id is an ID specified in the argument "+\
                             "map_ids above.")

    parser.add_argument("-C", "--combine_different_wafers",
                        action="store_true", default=False,
                        help="Whether to add maps made from different wafers "+\
                             "together. If True, maps whose 'Id's look like "+\
                             "map_id+a_wafer_name will be added, "+\
                             "where map_id is an ID specified in the argument "+\
                             "map_ids above.")

    parser.add_argument("-u", "--subtract_existing_maps",
                        action="store_true", default=False,
                        help="Whether to subtract certain observation's data "+\
                             "contained in an input file from the coadded maps "+\
                             "if the data already exist in the coadded maps.")

    parser.add_argument("-a", "--collect_averages_from_flagging_statistics",
                        action="store_true", default=False,
                        help="Whether to calculate average number of bolometers "+\
                             "that were not flagged over all scans of "+\
                             "an observation, average number of bolometers that "+\
                             "were removed, and so on.")
    
    parser.add_argument("-k", "--calculate_pW_to_K_conversion_factors",
                        action="store_true", default=False,
                        help="Whether to calculate the mean value of the "+\
                             "pW/K absolute calibration factors for each band "+\
                             "of each wafer.")

    parser.add_argument("-n", "--calculate_noise_from_individual_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels from "+\
                             "each map that goes into the coadded maps.")

    parser.add_argument("-N", "--calculate_noise_from_coadded_maps",
                        action="store_true", default=False,
                        help="Whether to calculate map noise levels for the "+\
                             "running coadded maps. In other words, the noise "+\
                             "will be calculated from the coadded maps "+\
                             "every time a new map is added.")
    
    parser.add_argument("-x", "--calculate_cross_spectrum_with_coadded_maps",
                        action="store_true", default=False,
                        help="Whether to calculate the cross spectrum between "+\
                             "each observation's map and the running "+\
                             "coadded maps.")

    parser.add_argument("-M", "--point_source_file",
                        type=str, action="store", default="",
                        help="Path to a point source list, which will be used "+\
                             "when making a mask for calcuting Cls.")

    parser.add_argument("-p", "--calculate_pointing_discrepancies",
                        action="store_true", default=False,
                        help="Whether to calculate the difference between the "+\
                             "true positions of some point sources and the "+\
                             "measured positions")

    parser.add_argument("-S", "--sources",
                        type=str, action="store", nargs="+",
                        default=["ra0hdec-44.75", "ra0hdec-52.25",
                                 "ra0hdec-59.75", "ra0hdec-67.25"],
                        help="The sub-field(s) that will be included "+\
                             "in the coadded maps.")

    parser.add_argument("-s", "--min_obs_id",
                        type=int, action="store", default=0,
                        help="The smallest observation ID that will be "+\
                             "considered to be used in the coadded map. "+\
                             "The default is 00000000.")

    parser.add_argument("-l", "--max_obs_id",
                        type=int, action="store", default=99999999,
                        help="The largest observation ID that will be "+\
                             "considered to be used in the coadded map. "+\
                             "The default is 99999999.")

    parser.add_argument("-b", "--bad_obs_ids",
                        type=int, action="store", nargs="+", default=[],
                        help="The observation IDs that will be excluded from "+\
                             "making the coadded map. The script has some "+\
                             "simple criteria to decide whether to include "+\
                             "certain observations, but one can manually specify "+\
                             "bad observations, too.")

    parser.add_argument("-m", "--min_file_size",
                        type=float, action="store", default=0.01,
                        help="The minimum size (in GB) a g3 file needs "+\
                             "to have for it to be considered as a good file. "+\
                             "This is to reduce the occurrence of a situation "+\
                             "where the script crashes due to reading a "+\
                             "problematic file.")
    
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


