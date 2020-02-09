// check if cookies are defined, otherwise set them to defaults
var default_cookies = {'wafer': 'all',
                       'weekdir': 'plots/last_n/last_3/',
                       'mapweekdir': 'maps/figures/last_n/last_7/',
                       'mapweekdirsummer': 'maps/figures_summer/last_n/last_7/'};
for (var cookie_name in default_cookies) {
    if (Cookies.get(cookie_name) == undefined) 
        Cookies.set(cookie_name, default_cookies[cookie_name], {expires: 1});
}

// get cookie values
var wafer = Cookies.get('wafer');
var weekdir = Cookies.get('weekdir');
var mapweekdir = Cookies.get('mapweekdir');
var mapweekdirsummer = Cookies.get('mapweekdirsummer');


/**
 * Updates the source of all images defined on the summary page.
 */
function update_figs() {
    document["ts_median_cal_sn_highel"].src        = 'staticimg/'+weekdir+'/median_cal_sn_4Hz_highel_'+wafer+'.png';
    document["ts_median_cal_response_highel"].src  = 'staticimg/'+weekdir+'/median_cal_response_4Hz_highel_'+wafer+'.png';
    document["ts_alive_bolos_cal_highel"].src      = 'staticimg/'+weekdir+'/alive_bolos_cal_4Hz_highel_'+wafer+'.png';
      
    document["ts_median_cal_sn_lowel"].src         = 'staticimg/'+weekdir+'/median_cal_sn_4Hz_lowel_'+wafer+'.png';
    document["ts_median_cal_response_lowel"].src   = 'staticimg/'+weekdir+'/median_cal_response_4Hz_lowel_'+wafer+'.png';
    document["ts_alive_bolos_cal_lowel"].src       = 'staticimg/'+weekdir+'/alive_bolos_cal_4Hz_lowel_'+wafer+'.png';
    
    document["ts_median_elnod_slope"].src          = 'staticimg/'+weekdir+'/median_elnod_response_'+wafer+'.png';
    document["ts_median_elnod_sn"].src             = 'staticimg/'+weekdir+'/median_elnod_sn_'+wafer+'.png';
    document["ts_median_elnod_opacity"].src        = 'staticimg/'+weekdir+'/median_elnod_opacity_'+wafer+'.png';
    document["ts_median_elnod_iq"].src             = 'staticimg/'+weekdir+'/median_elnod_iq_phase_'+wafer+'.png';
    document["ts_alive_bolos_elnod"].src           = 'staticimg/'+weekdir+'/alive_bolos_elnod_'+wafer+'.png';
    
    document["ts_rcw38_sky_transmission"].src      = 'staticimg/'+weekdir+'/rcw38_sky_transmission_'+wafer+'.png';
    
    document["ts_median_rcw38_fluxcal"].src        = 'staticimg/'+weekdir+'/median_rcw38_fluxcal_'+wafer+'.png';
    document["ts_median_rcw38_intflux"].src        = 'staticimg/'+weekdir+'/median_rcw38_intflux_'+wafer+'.png';
    document["ts_alive_bolos_rcw38"].src           = 'staticimg/'+weekdir+'/alive_bolos_rcw38_'+wafer+'.png';
    
    document["ts_mat5a_sky_transmission"].src      = 'staticimg/'+weekdir+'/mat5a_sky_transmission_'+wafer+'.png';
    
    document["ts_median_mat5a_fluxcal"].src        = 'staticimg/'+weekdir+'/median_mat5a_fluxcal_'+wafer+'.png';
    document["ts_median_mat5a_intflux"].src        = 'staticimg/'+weekdir+'/median_mat5a_intflux_'+wafer+'.png';
    document["ts_alive_bolos_mat5a"].src           = 'staticimg/'+weekdir+'/alive_bolos_mat5a_'+wafer+'.png';
    
    document["ts_median_net_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NET_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_net_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NET_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_net_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NET_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_net_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NET_10.0Hz_to_15.0Hz_'+wafer+'.png';
    
    document["ts_median_nep_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NEP_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_nep_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NEP_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_nep_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NEP_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_nep_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NEP_10.0Hz_to_15.0Hz_'+wafer+'.png';
    
    document["ts_median_nei_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NEI_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_nei_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NEI_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_nei_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NEI_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_nei_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NEI_10.0Hz_to_15.0Hz_'+wafer+'.png';
    
    document["map_t_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map.png'
    document["map_t_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map.png'
    document["map_t_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map.png'
    
    document["map_t_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map.png'
    document["map_t_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map.png'
    document["map_t_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map.png'
    
    document["map_tt_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-TT_weight_map.png'
    document["map_tt_150"].src = 'staticimg/'+mapweekdir+'/150GHz-TT_weight_map.png'
    document["map_tt_220"].src = 'staticimg/'+mapweekdir+'/220GHz-TT_weight_map.png'
    
    document["map_tt_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-TT_weight_map.png'
    document["map_tt_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-TT_weight_map.png'
    document["map_tt_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-TT_weight_map.png'

    document["xsec_tt_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-TT_weight_map_cross_sectional_view.png'
    document["xsec_tt_150"].src = 'staticimg/'+mapweekdir+'/150GHz-TT_weight_map_cross_sectional_view.png'
    document["xsec_tt_220"].src = 'staticimg/'+mapweekdir+'/220GHz-TT_weight_map_cross_sectional_view.png'

    document["xsec_tt_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-TT_weight_map_cross_sectional_view.png'
    document["xsec_tt_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-TT_weight_map_cross_sectional_view.png'
    document["xsec_tt_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-TT_weight_map_cross_sectional_view.png'
    
    document["individ_noise_t_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_noise_levels_from_individual_maps.png'
    document["individ_noise_t_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_noise_levels_from_individual_maps.png'
    document["individ_noise_t_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_noise_levels_from_individual_maps.png'
    
    document["individ_noise_t_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_noise_levels_from_individual_maps.png'
    document["individ_noise_t_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_noise_levels_from_individual_maps.png'
    document["individ_noise_t_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_noise_levels_from_individual_maps.png'

    document["mean_tt_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_mean_of_tt_weight_map_values.png'
    document["mean_tt_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_mean_of_tt_weight_map_values.png'
    document["mean_tt_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_mean_of_tt_weight_map_values.png'

    document["mean_tt_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_mean_of_tt_weight_map_values.png'
    document["mean_tt_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_mean_of_tt_weight_map_values.png'
    document["mean_tt_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_mean_of_tt_weight_map_values.png'
    
    document["running_noise_t_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_noise_levels_from_running_coadds.png'
    document["running_noise_t_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_noise_levels_from_running_coadds.png'
    document["running_noise_t_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_noise_levels_from_running_coadds.png'

    document["running_noise_t_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_noise_levels_from_running_coadds.png'
    document["running_noise_t_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_noise_levels_from_running_coadds.png'
    document["running_noise_t_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_noise_levels_from_running_coadds.png'
    
    document["distro_noise_t_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_distributions_of_individual_noise_levels.png'
    document["distro_noise_t_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_distributions_of_individual_noise_levels.png'
    document["distro_noise_t_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_distributions_of_individual_noise_levels.png'

    document["distro_noise_t_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distributions_of_individual_noise_levels.png'
    document["distro_noise_t_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_distributions_of_individual_noise_levels.png'
    document["distro_noise_t_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_distributions_of_individual_noise_levels.png'
    
    document["distro_mean_tt_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_distributions_of_individual_mean_tt_weights.png'
    document["distro_mean_tt_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_distributions_of_individual_mean_tt_weights.png'
    document["distro_mean_tt_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_distributions_of_individual_mean_tt_weights.png'

    document["distro_mean_tt_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distributions_of_individual_mean_tt_weights.png'
    document["distro_mean_tt_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_distributions_of_individual_mean_tt_weights.png'
    document["distro_mean_tt_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_distributions_of_individual_mean_tt_weights.png'
    
    document["ra_offsets_el0"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_delta_Ras_from_point_sources_in_ra0hdec-44.75.png'
    document["ra_offsets_el1"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_delta_Ras_from_point_sources_in_ra0hdec-52.25.png'
    document["ra_offsets_el2"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_delta_Ras_from_point_sources_in_ra0hdec-59.75.png'
    document["ra_offsets_el3"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_delta_Ras_from_point_sources_in_ra0hdec-67.25.png'

    document["ra_offsets_el0_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-24.5.png'
    document["ra_offsets_el1_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-31.5.png'
    document["ra_offsets_el2_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-38.5.png'
    document["ra_offsets_el3_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-45.5.png'
    document["ra_offsets_el4_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-52.5.png'
    document["ra_offsets_el5_summer"].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_delta_Ras_from_point_sources_in_ra5hdec-59.5.png'
    
    document["dec_offsets_el0"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_delta_Decs_from_point_sources_in_ra0hdec-44.75.png'
    document["dec_offsets_el1"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_delta_Decs_from_point_sources_in_ra0hdec-52.25.png'
    document["dec_offsets_el2"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_delta_Decs_from_point_sources_in_ra0hdec-59.75.png'
    document["dec_offsets_el3"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_delta_Decs_from_point_sources_in_ra0hdec-67.25.png'

    document["dec_offsets_el0_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-24.5.png'
    document["dec_offsets_el1_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-31.5.png'
    document["dec_offsets_el2_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-38.5.png'
    document["dec_offsets_el3_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-45.5.png'
    document["dec_offsets_el4_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-52.5.png'
    document["dec_offsets_el5_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_delta_Decs_from_point_sources_in_ra5hdec-59.5.png'
    
    document['distro_ra_offsets'].src  = 'staticimg/'+mapweekdir+'/150GHz-T_map_distribution_of_average_Ra_offsets.png'
    document['distro_dec_offsets'].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_distribution_of_average_Dec_offsets.png'

    document['distro_ra_offsets_summer'].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distribution_of_average_Ra_offsets.png'
    document['distro_dec_offsets_summer'].src = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distribution_of_average_Dec_offsets.png'
    
    document["flagging_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_average_numbers_of_flagged_detectors.png'
    document["flagging_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_average_numbers_of_flagged_detectors.png'
    document["flagging_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_average_numbers_of_flagged_detectors.png'

    document["flagging_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_average_numbers_of_flagged_detectors.png'
    document["flagging_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_average_numbers_of_flagged_detectors.png'
    document["flagging_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_average_numbers_of_flagged_detectors.png'
    
    document["fractional_coverage_90"].src  = 'staticimg/'+mapweekdir+
        '/90GHz-T_map_fractional_coverage_in_tt_weight_maps.png'
    document["fractional_coverage_150"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_fractional_coverage_in_tt_weight_maps.png'
    document["fractional_coverage_220"].src = 'staticimg/'+mapweekdir+
        '/220GHz-T_map_fractional_coverage_in_tt_weight_maps.png'

    document["fractional_coverage_90_summer"].src  = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_fractional_coverage_in_tt_weight_maps.png'
    document["fractional_coverage_150_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/150GHz-T_map_fractional_coverage_in_tt_weight_maps.png'
    document["fractional_coverage_220_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/220GHz-T_map_fractional_coverage_in_tt_weight_maps.png'
    
    document["xspec_ratio_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_averages_of_power_spectra_ratios.png'
    document["xspec_ratio_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_averages_of_power_spectra_ratios.png'
    document["xspec_ratio_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_averages_of_power_spectra_ratios.png'

    document["xspec_ratio_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_averages_of_power_spectra_ratios.png'
    document["xspec_ratio_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_averages_of_power_spectra_ratios.png'
    document["xspec_ratio_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_averages_of_power_spectra_ratios.png'
    
    document["distro_avg_xsp_ratio_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_distributions_of_crass_spectra_ratios.png'
    document["distro_avg_xsp_ratio_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_distributions_of_crass_spectra_ratios.png'
    document["distro_avg_xsp_ratio_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_distributions_of_crass_spectra_ratios.png'

    document["distro_avg_xsp_ratio_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distributions_of_crass_spectra_ratios.png'
    document["distro_avg_xsp_ratio_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_distributions_of_crass_spectra_ratios.png'
    document["distro_avg_xsp_ratio_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_distributions_of_crass_spectra_ratios.png'
    
    document["observation_durations"].src = 'staticimg/'+mapweekdir+'/observation_durations.png'

    document["observation_durations_summer"].src = 'staticimg/'+mapweekdirsummer+'/observation_durations.png'
    
    document["cal_vs_el_by_field_90"].src  = 'staticimg/'+mapweekdir+
        '/90GHz-T_map_median_cal_resp_percentage_changes_by_field.png'
    document["cal_vs_el_by_field_150"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_median_cal_resp_percentage_changes_by_field.png'
    document["cal_vs_el_by_field_220"].src = 'staticimg/'+mapweekdir+
        '/220GHz-T_map_median_cal_resp_percentage_changes_by_field.png'

    document["cal_vs_el_by_field_90_summer"].src  = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_median_cal_resp_percentage_changes_by_field.png'
    document["cal_vs_el_by_field_150_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/150GHz-T_map_median_cal_resp_percentage_changes_by_field.png'
    document["cal_vs_el_by_field_220_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/220GHz-T_map_median_cal_resp_percentage_changes_by_field.png'
    
    document["cal_vs_el_by_wafer_90"].src  = 'staticimg/'+mapweekdir+
        '/90GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'
    document["cal_vs_el_by_wafer_150"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'
    document["cal_vs_el_by_wafer_220"].src = 'staticimg/'+mapweekdir+
        '/220GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'

    document["cal_vs_el_by_wafer_90_summer"].src  = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'
    document["cal_vs_el_by_wafer_150_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/150GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'
    document["cal_vs_el_by_wafer_220_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/220GHz-T_map_median_cal_resp_percentage_changes_by_wafer.png'
    
    document["distro_cal_vs_el_90"].src  = 'staticimg/'+mapweekdir+'/90GHz-T_map_distributions_of_cal_response_changes.png'
    document["distro_cal_vs_el_150"].src = 'staticimg/'+mapweekdir+'/150GHz-T_map_distributions_of_cal_response_changes.png'
    document["distro_cal_vs_el_220"].src = 'staticimg/'+mapweekdir+'/220GHz-T_map_distributions_of_cal_response_changes.png'

    document["distro_cal_vs_el_90_summer"].src  = 'staticimg/'+mapweekdirsummer+'/90GHz-T_map_distributions_of_cal_response_changes.png'
    document["distro_cal_vs_el_150_summer"].src = 'staticimg/'+mapweekdirsummer+'/150GHz-T_map_distributions_of_cal_response_changes.png'
    document["distro_cal_vs_el_220_summer"].src = 'staticimg/'+mapweekdirsummer+'/220GHz-T_map_distributions_of_cal_response_changes.png'
    
    document["pw_per_k_by_wafer_90"].src  = 'staticimg/'+mapweekdir+
        '/90GHz-T_map_median_temperature_calibration_factors_by_wafer.png'
    document["pw_per_k_by_wafer_150"].src = 'staticimg/'+mapweekdir+
        '/150GHz-T_map_median_temperature_calibration_factors_by_wafer.png'
    document["pw_per_k_by_wafer_220"].src = 'staticimg/'+mapweekdir+
        '/220GHz-T_map_median_temperature_calibration_factors_by_wafer.png'

    document["pw_per_k_by_wafer_90_summer"].src  = 'staticimg/'+mapweekdirsummer+
        '/90GHz-T_map_median_temperature_calibration_factors_by_wafer.png'
    document["pw_per_k_by_wafer_150_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/150GHz-T_map_median_temperature_calibration_factors_by_wafer.png'
    document["pw_per_k_by_wafer_220_summer"].src = 'staticimg/'+mapweekdirsummer+
        '/220GHz-T_map_median_temperature_calibration_factors_by_wafer.png'
}


/**
 * Show/hide plots when switching between yearly and non-yearly view for maps
 */
function update_visibility() {
    if(!window["mapweekdir"].includes('yearly')) {
        $(".maps_yearly").hide();
        $(".maps_non_yearly").show();
    }
    else {
        $(".maps_yearly").show();
        $(".maps_non_yearly").hide();
    }
    if(!window["mapweekdirsummer"].includes('yearly')) {
        $(".maps_summer_yearly").hide();
        $(".maps_summer_non_yearly").show();
    }
    else {
        $(".maps_summer_yearly").show();
        $(".maps_summer_non_yearly").hide();
    }
}

function update_lastmodified() {
	$.get('/lastmodified_calibration', function(data) {
		$("#lastmodified_calibration").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
	});
	$.get('/lastmodified_maps_winter', function(data) {
		$("#lastmodified_maps_winter").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
	});
	$.get('/lastmodified_maps_summer', function(data) {
		$("#lastmodified_maps_summer").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
	});
}

/**
 * Sets a global variable to records its value as a cookie for retrieval later.
 */
function set_variable(variable, newVal)
{
  window[variable] = newVal;
  Cookies.set(variable, newVal, { expires: 1 });
  update_figs();
  update_visibility();
}

/**
 * Builds the buttons on the summary page that select different time intervals
 * of data to display. These are constructed by appending the DOM directly.
 * Note also that this function also initializes the jQuery UI "controlgroup"
 * and binds the click event to it after appending all the buttons.
 * 
 * Arguments:
 *  interval : time interval to traverse, 'weekly' or 'monthly'
 *  subdirectory : subdirectory that we should traverse to get dated data
 *  tab : tab in which to add buttons, 'summary' or 'maps'
 */
function add_date_buttons(interval, subdirectory, tab, boundVar)
{
    // now rebuild the div
	$.get('/staticdirs', {subdirectory:subdirectory, interval:interval},
		  function(data, status) {
			  div_id = '#datalist_'+tab+'_'+interval;

			  data.reverse();
			  for (jdir=0; jdir<data.length; jdir++)
			  {
				  if (tab === 'summary')
					  datestring = data[jdir].split('/')[2];
				  else if (tab == 'maps')
					  datestring = data[jdir].split('/')[3];
                                  else if (tab == 'mapssummer')
                                          datestring = data[jdir].split('/')[3];

				  $(div_id).append("<input type='radio' id='dates-"+tab+"-"+datestring+"' name='dates' value='" +
				  				   data[jdir] + "'>\n" + 
				  				   "<label for='dates-"+tab+"-"+datestring+"'>"+datestring+"</label>");

				  // can use jquery <select ...></select> with this line instead of radio buttons
				  // $(div_id).append("<option value='"+data[jdir]+"'>"+datestring+"</option>");
			  }
			  $(div_id).controlgroup();
			  $("[id^=dates-"+tab+"-]").click(function(event) {
				  set_variable(boundVar, event.target.value);
			  });
		  });
}

// Page initialization
$( document ).ready(function()
{
	// Initialize jQuery UI elements and make dynamic modifications to the DOM
	$("#tabs").tabs();
	$("#waferlist").controlgroup();
	add_date_buttons('last_n', 'plots', 'summary', 'weekdir');
	add_date_buttons('monthly', 'plots', 'summary', 'weekdir');
	add_date_buttons('weekly', 'plots', 'summary', 'weekdir');

	add_date_buttons('last_n', 'maps/figures', 'maps', 'mapweekdir');
	add_date_buttons('yearly', 'maps/figures', 'maps', 'mapweekdir');
	add_date_buttons('monthly', 'maps/figures', 'maps', 'mapweekdir');
	add_date_buttons('weekly', 'maps/figures', 'maps', 'mapweekdir');

	add_date_buttons('last_n', 'maps/figures_summer', 'mapssummer', 'mapweekdirsummer');
	add_date_buttons('yearly', 'maps/figures_summer', 'mapssummer', 'mapweekdirsummer');
	add_date_buttons('monthly', 'maps/figures_summer', 'mapssummer', 'mapweekdirsummer');
	add_date_buttons('weekly', 'maps/figures_summer', 'mapssummer', 'mapweekdirsummer');

	
	// Bind the click event to the wafer buttons
	$("[id^=wafers-]").click(function(event) {
		set_variable("wafer", event.target.value);
	});

	// Bind the click event to the tabs so that we can keep track of the active tab
	$("#tabs").click(function(event) {
		Cookies.set('activeTab', $("#tabs").tabs("option", "active"), {expires: 1});
	});

	if(Cookies.get('activeTab') !== undefined)
		$("#tabs").tabs("option", "active", Cookies.get('activeTab'))
	else {
		Cookies.set('activeTab', $("#tabs").tabs("option", "active"), {expires: 1});
	}


	// Update all the figures. Need to call this on load because the values of
	// wafer and weekdir might be pulled from a cookie, which probably differs
	// from the default values.
	update_figs();
	update_visibility();
	update_lastmodified();
});
