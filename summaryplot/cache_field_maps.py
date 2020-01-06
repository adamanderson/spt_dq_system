# ============================================================================ #
#  This script is intended to be called by spt_dq_system/update_summary.py.    #
#  When update_summary.py calls this script, the former specifies              #
#  whether we want to process individual maps and then co-add them or          #
#  we want to make figures showing maps and their related quantities,          #
#  what time intervals we are interested in using maps from, and so on.        #
#                                                                              #
#  Then, Based on the specifications, this script in turn calls another        #
#  script, which is either fields_coadding.py or fields_plotting.py,           #
#  with approprite arguments.                                                  #
#                                                                              #
# ============================================================================ #


import argparse
import datetime
import os
import glob
import shutil
import json
import sys
import logging
import matplotlib
matplotlib.use("Agg")
import numpy

from spt3g import core
from spt3g import std_processing
from summaryplot import fields_coadding
from summaryplot import fields_plotting



# ==============================================================================
# Define the function that calls the map coadding and plotting scripts
# with appropriate arguments
# ------------------------------------------------------------------------------


def update(season, mode, action,
           oldest_time_to_consider=None, current_time=None,
           time_interval=None, last_how_many_days=0,
           original_maps_dir='.',
           calibration_data_dir='.', bolo_timestreams_dir='.',
           coadds_dir='.', figs_dir='.',
           bands=['90GHz', '150GHz', '220GHz'],
           sub_fields=['ra0hdec-44.75', 'ra0hdec-52.25',
                       'ra0hdec-59.75', 'ra0hdec-67.25'],
           logger_name='', log_file=None,
           just_see_commands=False):
    
    # - Define global variables
    
    # -- Variables related to logging
    
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
    
    log = logger.info
    
    log("")
    log("# --------------------------------------------- #")
    log("#  The script cache_field_maps.py was invoked!  #")
    log("# --------------------------------------------- #")
    log("")
    
    # -- Variables related to file i/o
    
    if current_time is None:
        current_time = datetime.datetime.utcnow()
    else:
        current_time = datetime.datetime.strptime(current_time, '20%y%m%d')
        current_time = current_time.replace(hour=23, minute=59, second=59)
    
    beginning_of_the_record = \
        datetime.datetime.strptime(oldest_time_to_consider, '20%y%m%d')
    
    desired_obs_id_ranges = []
    desired_time_ranges   = []
    desired_dir_names     = []
    
    script_coadding_maps = 'summaryplot/fields_coadding.py'
    script_plotting_data = 'summaryplot/fields_plotting.py'
    
    aux_files_directory  = 'summaryplot/fields_aux_files/'
    bad_map_list_file    = \
        os.path.join(aux_files_directory, 'bad_map_list.txt')
    planck_map_fits_file = \
        os.path.join(aux_files_directory,
                    'HFI_SkyMap_BAND_2048_R3.01_MISSION_cut_C_G3Units.fits')
    point_source_list_file = \
        'spt3g_software/sources/1500d_ptsrc_3band_50mJy.txt'
    dummy_input_file = \
        os.path.join(aux_files_directory, 'dummy.g3')
    
    # -- Variables that depend on the season
    
    if season == 'summer':
        all_sub_fields = ['ra5hdec-24.5', 'ra5hdec-31.5', 'ra5hdec-38.5',
                          'ra5hdec-45.5', 'ra5hdec-52.5', 'ra5hdec-59.5']
        if (len(sub_fields) != 0) and (set(sub_fields) < set(all_sub_fields)):
            # Probably a debugging mode.
            pass
        else:
            sub_fields = all_sub_fields
        planck_map_fits_file = \
            os.path.join(
                aux_files_directory,
                'HFI_SkyMap_BAND_2048_R3.01_MISSION_cut_summer_C_G3Units.fits')
        point_source_list_file = \
            os.path.join(
                aux_files_directory,
                '1500d_summer_point_source_list_from_at20g.txt')
    
    
    # - Figure out what appropriate time intervals are and
    #   make sure an appropriate directory structure exists
    
    # -- Figure out appropriate time intervals
    
    def convert_to_obs_id(datetime_object):
        
        return std_processing.time_to_obsid(
                   datetime_object.strftime("20%y%m%d_%H%M%S"))
    
    
    if time_interval == "last_n":
        delta_t    = datetime.timedelta(days=-1*last_how_many_days)
        end_time   = current_time
        start_time = current_time + delta_t
        if start_time < beginning_of_the_record:
            if mode == "coadding":
                start_time = beginning_of_the_record
        start_time = start_time.replace(minute=0, second=0, microsecond=0)
        
        end_obs_id   = convert_to_obs_id(end_time)
        start_obs_id = convert_to_obs_id(start_time)
        
        desired_obs_id_ranges.append((start_obs_id, end_obs_id))
        desired_time_ranges.append((start_time, end_time))
    
    else:        
        def get_the_beginning_of_the_interval(interval, end_time):
            
            if interval == "yearly":
                start_of_interval = end_time.replace(
                                        month=1, day=1, hour=0,
                                        minute=0, second=0, microsecond=0)
            if interval == "monthly":
                start_of_interval = end_time.replace(
                                        day=1, hour=0,
                                        minute=0, second=0, microsecond=0)
            if interval == "weekly":
                current_weekday   = end_time.weekday()
                delta_t           = datetime.timedelta(days=-1*current_weekday)
                start_of_interval = (end_time + delta_t).replace(
                                        hour=0,
                                        minute=0, second=0, microsecond=0)
            return start_of_interval
        
        
        if mode == "plotting":
            if time_interval == "yearly":
                end_time_at_start_of_loop = current_time
            if time_interval == "monthly":
                if current_time.month != 12:
                    beginning_of_next_month = \
                        current_time.replace(
                            month=current_time.month+1, day=1,
                            hour=0, minute=0, second=0, microsecond=0)
                else:
                    beginning_of_next_month = \
                        current_time.replace(year=current_time.year+1,
                            month=1, day=1,
                            hour=0, minute=0, second=0, microsecond=0)
                end_of_this_month = \
                    beginning_of_next_month + \
                    datetime.timedelta(seconds=-1)
                end_time_at_start_of_loop = end_of_this_month
            if time_interval == "weekly":
                end_of_this_week = \
                    current_time + \
                    datetime.timedelta(days=(6-current_time.weekday()))
                end_of_this_week = \
                    end_of_this_week.replace(
                        hour=23, minute=59, second=59, microsecond=0)
                end_time_at_start_of_loop = end_of_this_week
        if mode == "coadding":
            end_time_at_start_of_loop = current_time
        
        while end_time_at_start_of_loop > beginning_of_the_record:
            beginning_of_the_interval = \
                get_the_beginning_of_the_interval(
                    time_interval, end_time_at_start_of_loop)
            if (mode == "coadding") and \
               (beginning_of_the_interval < beginning_of_the_record):
                time_to_use_as_start = beginning_of_the_record
            else:
                time_to_use_as_start = beginning_of_the_interval
            
            start_obs_id_of_the_interval = \
                convert_to_obs_id(time_to_use_as_start)
            end_obs_id_of_the_interval = \
                convert_to_obs_id(end_time_at_start_of_loop)
            
            desired_obs_id_ranges.append(
                (start_obs_id_of_the_interval, end_obs_id_of_the_interval))
            desired_time_ranges.append(
                (beginning_of_the_interval, end_time_at_start_of_loop))
            
            end_time_at_start_of_loop = \
                beginning_of_the_interval + datetime.timedelta(seconds=-1)
    
    # -- Check the input and output directory structure
    
    log("")
    log("--------------------------------------------------")
    log(" Making sure the directory structure is valid ... ")
    log("--------------------------------------------------")
    log("")
    
    log("* Checking whether relevant directories exist or not ...")
    log("  (If not, they will be created.)")
    
    def convert_time_intervals_to_dir_names(interval_type, datetime_obj):
        
        if interval_type == "weekly":
            str_fmt = "20%y%m%d"
        elif interval_type == "monthly":
            str_fmt = "20%y%m"
        elif interval_type == "yearly":
            str_fmt = "20%y"
        return datetime_obj.strftime(str_fmt)
    
    
    for root_directory in [coadds_dir, figs_dir]:
        if not os.path.isdir(root_directory):
            os.mkdir(root_directory)
        if time_interval not in os.listdir(root_directory):
            os.mkdir(os.path.join(root_directory, time_interval))
        for time_range in desired_time_ranges:
            if time_interval == "last_n":
                sub_dir = time_interval.replace("n", str(last_how_many_days))
            else:
                sub_dir = convert_time_intervals_to_dir_names(
                              time_interval, time_range[0])
            full_path = os.path.join(root_directory, time_interval, sub_dir)
            desired_dir_names.append(full_path)
    
    for dir_name in desired_dir_names:
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
    
    log("")
    
    
    # - Call the relevant scripts with appropriate arguments
    
    log("")
    log("----------------------------------------------------------")
    log(" The mode, action, and time interval that were specified: ")
    log("----------------------------------------------------------")
    log("")
    log("   %s, %s, %s", mode, action, time_interval)
    log("")
    
    log("")
    log("----------------------------------------------")
    log(" Relevant obs_id range(s) to take actions on: ")
    log("----------------------------------------------")
    log("")
    for counter, interval in enumerate(desired_obs_id_ranges, 1):
        log("Interval %s :", counter)
        log("  from %s ("+\
            str(std_processing.obsid_to_g3time(interval[0])).split(".")[0]+")",
            interval[0])
        log("  to   %s ("+\
            str(std_processing.obsid_to_g3time(interval[1])).split(".")[0]+")",
            interval[1])
    log("")
    
    log("")
    log("---------------------------")
    log(" Relevant commands to run: ")
    log("---------------------------")
    log("")
    
    # -- First, define some functions
    
    def decide_what_ids_to_use_and_not(
            existing_ids, desired_id_range, file_name):
        
        if len(existing_ids) == 0:
            raise RuntimeError("It appears that coadded maps stored in "
                               + file_name + "are empty. Probably something "
                               "went wrong last time!")
        
        max_eid = numpy.max(existing_ids)
        min_eid = numpy.min(existing_ids)
        max_gid = desired_id_range[1]  # gid = good ids!
        min_gid = desired_id_range[0]
        
        if max_eid<min_gid:
            ignore_coadds  = True
            ids_to_exclude = []
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid
        
        elif (min_eid<min_gid) and (max_eid>min_gid) and (max_eid<max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids if obs_id >= min_gid]
            subtract_maps  = True
            id_lower_bound = min_eid
            id_upper_bound = max_gid
        
        elif (min_eid<min_gid) and (max_eid>max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids \
                 if (obs_id >= min_gid) and (obs_id <= max_gid)]
            subtract_maps  = True
            id_lower_bound = min_eid
            id_upper_bound = max_eid
        
        elif (min_eid>=min_gid) and (max_eid<=max_gid):
            ignore_coadds  = False
            ids_to_exclude = existing_ids
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid
        
        elif (min_eid>=min_gid) and (min_eid<max_gid) and (max_eid>max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids if obs_id <= max_gid]
            subtract_maps  = True
            id_lower_bound = min_gid
            id_upper_bound = max_eid
        
        elif min_eid>max_gid:
            ignore_coadds  = True
            ids_to_exclude = []
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid
        
        else:
            raise RuntimeError("There seems to be a forgotten scenario about "
                               "the relation between existing observation IDs "
                               "and desired ones to use!")
        
        return ignore_coadds,  subtract_maps, \
               ids_to_exclude, int(id_lower_bound), int(id_upper_bound)
    
    
    
    def get_bad_obs_ids_from_the_list(
            text_file, band_to_compare, sub_field_to_compare, old_bad_ids):
        
        ids_to_exclude = []
        try:
            bands, sub_fields, obs_ids = \
                numpy.loadtxt(text_file, dtype=str, delimiter="|",
                              usecols=(0, 1, 2), unpack=True)
        except:
            bands      = []
            sub_fields = []
            obs_ids    = []
        
        for index, band in enumerate(bands):
            band = band.replace(" ", "")
            if band == band_to_compare:
                sub_field = sub_fields[index].replace(" ", "")
                if sub_field == sub_field_to_compare:
                    ids_to_exclude.append(int(obs_ids[index].replace(" ", "")))
        
        rd = {'bad_obs_ids': old_bad_ids + ids_to_exclude}
        
        return rd
    
    
    
    def figure_out_io_arguments_for_coadding(
            band, sub_field, rg_idx, n_iter=1):
        
        arguments = {}
        
        arguments['sub_fields'] = [sub_field]
        arguments['temperature_maps_only'] = True
        
        g3_files_for_individual_observations = \
            sorted(glob.glob(os.path.join(original_maps_dir, sub_field,
                             '*{}_tonly.g3.gz'.format(band))))
        g3_file_for_coadded_maps = \
            os.path.join(desired_dir_names[rg_idx],
                         'coadded_maps_from_{}_{}.g3.gz'.\
                         format(sub_field, band))
        
        arguments['output_file'] = g3_file_for_coadded_maps.replace('g3', 'g4')
        
        if (action == 'update') and \
           (os.path.isfile(g3_file_for_coadded_maps)):
            iterator = core.G3File(g3_file_for_coadded_maps)
            existing_ids = []
            
            for iteration_time in range(n_iter):
                frame = iterator.next()
                for obs_id_types \
                in  ["CoaddedObservationIDs", "IgnoredObservationIDs"]:
                    for sf, oids in frame[obs_id_types].items():
                        for oid in oids:
                            if oid not in existing_ids:
                                existing_ids.append(oid)
            
            ignore_coadds,  subtract_maps, \
            ids_to_exclude, id_lower_bound, id_upper_bound = \
                decide_what_ids_to_use_and_not(
                    existing_ids, desired_obs_id_ranges[rg_idx],
                    g3_file_for_coadded_maps)
            
            if ignore_coadds:
                arguments['input_files'] = g3_files_for_individual_observations
            else:
                arguments['input_files'] = [g3_file_for_coadded_maps] + \
                                           g3_files_for_individual_observations
            if subtract_maps:
                arguments['subtract_existing_maps'] = True
            else:
                arguments['subtract_existing_maps'] = False
            arguments['min_obs_id']  = id_lower_bound
            arguments['max_obs_id']  = id_upper_bound
            arguments['bad_obs_ids'] = ids_to_exclude
        
        else:
            arguments['input_files'] = g3_files_for_individual_observations
            arguments['min_obs_id']  = desired_obs_id_ranges[rg_idx][0]
            arguments['max_obs_id']  = desired_obs_id_ranges[rg_idx][1]
        
        arguments['bad_map_list_file'] = bad_map_list_file
        arguments['log_file']    = log_file
        arguments['logger_name'] = \
            '{}_{}_{}'.format(
            sub_logger_name, band, sub_field.replace('.', ''))
        arguments['less_verbose'] = False
        
        arguments['calibration_data_dir'] = calibration_data_dir
        arguments['bolo_timestreams_dir'] = bolo_timestreams_dir
        arguments['auxiliary_files_directory'] = aux_files_directory
        arguments['point_source_list_file']    = point_source_list_file
        arguments['planck_map_fits_file']      = planck_map_fits_file
        
        return arguments
    
    
    
    def run_coadding_or_plotting_function(
            function, arguments, logger, log_file, just_see_args):
        
        log('Calling the function ...')
        log('Here are the arguments to be supplied:')
        arguments_to_show = {k: v for k, v in arguments.items()}
        for key in ['input_files', 'bad_obs_ids']:
            if key in arguments_to_show.keys():
                if len(arguments_to_show[key]) > 10:
                    arguments_to_show.pop(key)
        log(json.dumps(arguments_to_show, indent=1))
        log('')
        
        if not just_see_args:
            if log_file is None:
                rval = function(**arguments)
            else:
                try:
                    rval = function(**arguments)
                except Exception:
                    logger.exception('Something did not go well '
                                     'while running coadding/plotting!')
                    raise RuntimeError('An error occurred! '
                                       'Please check where it occurred '
                                       'in the log file {}!'.format(log_file))
        else:
            rval = 1
        log('')
        
        return rval
    
    
    
    def analyze_and_coadd_maps(
            function, arguments, logger, log_file, just_see_commands):
        
        output_file_temp = arguments['output_file']
        output_file_perm = output_file_temp.replace('g4', 'g3')
        
        rval = run_coadding_or_plotting_function(
                   function, arguments,
                   logger, log_file, just_see_commands)
        
        log('Changing the name of the temporary output file ...')
        if os.path.isfile(output_file_temp):
            log('  (Original name: %s' , output_file_temp)
            log('   New name     : %s)', output_file_perm)
            if not just_see_commands:
                if os.path.getsize(output_file_temp) < 100 * 2**20:
                    logger.exception(
                        'The size of the output file %s seems too small!',
                        output_file_temp)
                    raise RuntimeError(
                        'An error occurred! '
                        'Please check where it occurred '
                        'in the log file {}!'.format(log_file))
                shutil.move(output_file_temp, output_file_perm)
            log('Done.')
        else:
            log('Actually, %s does not exist.', output_file_temp)
        log('')
        
        return rval
    
    
    
    def save_analysis_results_in_separate_file(original_file):
        
        def load_analysis_results_only(frame, original_file):
            if frame.type != core.G3FrameType.EndProcessing:
                return []
            else:
                final_frame = core.G3File(original_file).next()
                return [final_frame, frame]
        
        output_file = original_file.\
                          replace('coadded_maps', 'some_analysis_results').\
                          replace('.gz', '')
        
        pipeline = core.G3Pipeline()
        pipeline.Add(core.G3Reader,
                     filename=dummy_input_file)
        pipeline.Add(load_analysis_results_only,
                     original_file=original_file)
        pipeline.Add(core.G3Writer,
                     filename=output_file)
        pipeline.Run()
    
    
    
    def combine_all_analysis_results(
            analysis_results_g3_files, output_file):
        
        def combine_analysis_results(frame, g3_files=[]):
            if frame.type != core.G3FrameType.EndProcessing:
                return []
            else:
                final_frame = core.G3Frame()
                for g3_file in g3_files:
                    iterator  = core.G3File(g3_file)
                    new_frame = iterator.next()
                    for k, v in new_frame.iteritems():
                        if k == "Id":
                            if "Id" not in final_frame.keys():
                                final_frame["Id"] = new_frame["Id"]
                        if not isinstance(v, core.G3MapMapDouble):
                            continue
                        if k not in final_frame.keys():
                            final_frame[k] = v
                        else:
                            old_v = {"some_id": final_frame[k]}
                            new_v = {"some_id": new_frame[k]}
                            combined_kv = \
                                fields_coadding.    \
                                AnalyzeAndCoaddMaps.\
                                combine_mapmapdoubles(
                                    "dummy", old_v, new_v)
                            final_frame.pop(k)
                            final_frame[k] = combined_kv["some_id"]
                return [final_frame, frame]
        
        pipeline = core.G3Pipeline()
        pipeline.Add(core.G3Reader,
                     filename=dummy_input_file)
        pipeline.Add(combine_analysis_results,
                     g3_files=analysis_results_g3_files)
        pipeline.Add(core.G3Writer,
                     filename=output_file)
        pipeline.Run()
    
    
    
    # -- Then, call the coadding/plotting scripts by using those functions
    
    n_time_ranges = len(desired_obs_id_ranges)
    for i in range(n_time_ranges):
        
        log("# Command set %s to run:", "{:03}".format(i+1))
        log("# ---------------------------")
        log("")
        
        sub_logger_name = \
            '{}_{}_{}'.format(
            mode, time_interval, desired_dir_names[i].split('/')[-1])
        
        
        if mode == 'coadding':
            
            for band in bands:
                rvals_sum = 0
                
                anal_yearly_args = \
                    {'map_ids': ['Left'+band, 'Right'+band],
                     'trick_pipeline_to_receive_left_right_maps': True,
                     'ids_to_append_directions_to'              : [band],
                     'maps_split_by_scan_direction'             : True,
                     'combine_left_right'                       : False,
                     'calculate_map_rmss_and_weight_stats'      : True,
                     'rmss_and_wgts_from_coadds_or_individuals' : 'c',
                     'rmss_and_wgts_from_signals_or_noises'     : 'sn',
                     'calculate_noise_from_coadded_maps'        : True}
                
                anal_mothly_args = \
                    {'map_ids': [band],
                     'calculate_map_rmss_and_weight_stats'      : True,
                     'rmss_and_wgts_from_coadds_or_individuals' : 'i',
                     'rmss_and_wgts_from_signals_or_noises'     : 's',
                     'collect_averages_from_flagging_statistics': True,
                     'calculate_pW_to_K_conversion_factors'     : True,
                     'calculate_calibrator_response_vs_el'      : True,
                     'calculate_cross_spectra_with_planck_map'  : True,
                     'calculate_pointing_discrepancies'         : True,
                     'calculate_noise_from_individual_maps'     : True}
                
                anal_simple_args = \
                    {'map_ids': [band]}
                
                
                for sub_field in sub_fields:
                    log("--- %s %s ---", band, sub_field)
                    log("")
                    
                    n_iter = 2 if time_interval == 'yearly' else 1
                    args_coadding = \
                        figure_out_io_arguments_for_coadding(
                            band, sub_field, i, n_iter)
                    
                    if time_interval == "yearly":
                        args_coadding.update(anal_yearly_args)
                    elif time_interval == "monthly":
                        args_coadding.update(anal_mothly_args)
                    else:
                        args_coadding.update(anal_simple_args)
                    
                    if time_interval != 'last_n':
                        # * If existing maps are not to be subtracted,
                        #   then there is no need to modify the id range
                        args_coadding.update(
                            {'subtract_existing_maps': False})
                        args_coadding.update(
                            {'min_obs_id': desired_obs_id_ranges[i][0]})
                        args_coadding.update(
                            {'max_obs_id': desired_obs_id_ranges[i][1]})
                    
                    rval = analyze_and_coadd_maps(
                               fields_coadding.run, args_coadding,
                               logger, log_file, just_see_commands)
                    rvals_sum += rval
                    
                    log('\n\n')
                
                log('--- %s Full field ---', band)
                log('')
                log('\n')
                
                """log('Saving the analysis results in a separate file ...')
                if not just_see_commands:
                    save_analysis_results_in_separate_file(
                        os.path.join(desired_dir_names[i],
                                     'coadded_maps_{}.g3.gz'.format(band)))
                log('')"""
                
                if rvals_sum == 0:
                    log('There was no update to any of the sub-field,')
                    log('so, there is no need to combine the 4 g3 files again!')
                    log('\n\n')
                    continue
                
                args_coadding = \
                    {'map_ids'              : [band],
                     'sub_fields'           : sub_fields,
                     'temperature_maps_only': True,
                     'input_files': [os.path.join(desired_dir_names[i],
                                     'coadded_maps_from_{}_{}.g3.gz'.\
                                      format(sf, band)) for sf in sub_fields],
                     'output_file':  os.path.join(desired_dir_names[i],
                                     'coadded_maps_{}.g3.gz'.format(band)),
                     'logger_name': '{}_{}_full_field'.\
                                     format(sub_logger_name, band),
                     'log_file'   : log_file,
                     'auxiliary_files_directory': aux_files_directory,
                     'calibration_data_dir'     : calibration_data_dir,
                     'bolo_timestreams_dir'     : bolo_timestreams_dir,
                     'bad_map_list_file'        : bad_map_list_file}
                
                if time_interval == 'yearly':
                    args_coadding.update(anal_yearly_args)
                elif time_interval == 'monthly':
                    args_coadding.update(anal_mothly_args)
                else:
                    args_coadding.update(anal_simple_args)
                
                if time_interval == 'yearly':
                    args_coadding.update({'combine_left_right': True,
                                          'map_ids'           : [band]})
                
                run_coadding_or_plotting_function(
                    fields_coadding.run, args_coadding,
                    logger, log_file, just_see_commands)
                
                log('\n')
                log('Saving the analysis results in a separate file ...')
                if not just_see_commands:
                    if os.path.isfile(args_coadding['output_file']):
                        save_analysis_results_in_separate_file(
                            args_coadding['output_file'])
                log('')
                
                log('\n\n')
        
        
        elif mode == 'plotting':
            
            fig_dir = desired_dir_names[i+n_time_ranges]
            
            for band in bands:
                log('--- %s ---', band)
                log('')
                
                obs_info_etc_file = \
                    os.path.join(desired_dir_names[i],
                                 'some_analysis_results_{}.g3'.format(band))
                analysis_results_etc_file = \
                    os.path.join(coadds_dir, 'monthly', 'all_months',
                                 'all_analysis_results_{}.g3'.format(band))
                map_data_etc_file = \
                    os.path.join(desired_dir_names[i],
                                 'coadded_maps_{}.g3.gz'.format(band))
                
                ignore_map_data = False   # * Just for temporary testing
                if ignore_map_data:
                    input_files = [obs_info_etc_file,
                                   analysis_results_etc_file,
                                   obs_info_etc_file]
                else:
                    input_files = [obs_info_etc_file,
                                   analysis_results_etc_file,
                                   map_data_etc_file]
                    if not os.path.isfile(map_data_etc_file):
                        log('The file that is supposed to contain maps '
                            'does not exist yet, so nothing to do!')
                        log('\n\n')
                        continue
                
                args_plotting = \
                    {'input_files': input_files,
                     'directory_to_save_figures': fig_dir,
                     'simpler_file_names'       : True,
                     'figure_title_font_size'   : 22,
                     'map_id'                   : band,
                     'map_type'                 : 'T',
                     'coadded_data'             : True,
                     'make_figure_for_field_map': True,
                     'smooth_map_with_gaussian' : False,
                     'make_figure_for_entire_weight_map'       : True,
                     'make_figure_for_weight_map_cross_section': True,
                     'logger_name': '{}_{}'.format(sub_logger_name, band),
                     'log_file'   : log_file}
                
                if time_interval in ["last_n", "weekly"]:
                    args_plotting.update(
                        {'smooth_map_with_gaussian': True,
                         'gaussian_fwhm'           : 1.0})
                
                if time_interval == 'yearly':
                    args_plotting.update(
                        {'make_figures_showing_time_evolution'     : True,
                         'make_figures_showing_distributions'      : True,
                         'make_figures_for_responsivity_changes'   : True,
                         'make_figures_for_fluctuation_metrics'    : True,
                         'make_figures_for_pointing_discrepancies' : True,
                         'make_figure_for_ratios_of_power_spectra': True,
                         'make_figure_for_noise_levels'           : True})
                else:
                    args_plotting.update(
                        {'make_figures_showing_time_variations'   : True,
                         'make_figure_for_flagging_statistics'    : True,
                         'make_figure_for_pW_to_K_factors'        : True,
                         'make_figures_for_responsivity_changes'  : True,
                         'make_figures_for_fluctuation_metrics'   : True,
                         'make_figures_for_pointing_discrepancies': True,
                         'make_figure_for_ratios_of_power_spectra': True,
                         'make_figure_for_noise_levels'           : True,
                         'left_xlimit_for_time_variations' : \
                              desired_obs_id_ranges[i][0],
                         'right_xlimit_for_time_variations': \
                              desired_obs_id_ranges[i][1]})
                
                if action == 'update':
                    args_plotting.update(
                        {'decide_whether_to_make_figures_at_all': True})
                
                run_coadding_or_plotting_function(
                    fields_plotting.run, args_plotting,
                    logger, log_file, just_see_commands)
                
                log('\n\n')
            
            if len(os.listdir(fig_dir)) == 0:
                os.rmdir(fig_dir)
            
        log('\n\n\n')
    
    
    if (mode == 'coadding') and (time_interval == 'monthly'):
    
        for band in bands:
            
            all_months_dir = \
                os.path.join(coadds_dir, 'monthly', 'all_months')
            if not os.path.isdir(all_months_dir):
                os.mkdir(all_months_dir)
            
            all_months_g3_files = \
                sorted(glob.glob(os.path.join(
                           coadds_dir, 'monthly', '*',
                           'coadded_maps_{}.g3.gz'.format(band))))
            
            output_file = \
                os.path.join(
                    all_months_dir,
                    'all_analysis_results_{}.g3'.format(band))
            
            log('\n')
            log('Finally, the analysis results from all months '
                'will be combined into one file ...')
            log('')
            log(' * Input files:')
            for inf in all_months_g3_files:
                log('     %s', os.path.relpath(inf, coadds_dir))
            log(' * Output file:')
            log('     %s', os.path.relpath(output_file, coadds_dir))
            
            if not just_see_commands:
                combine_all_analysis_results(
                    all_months_g3_files, output_file)
            
            log('')
            log('Done.')
            log('\n\n')
                

# ==============================================================================




# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                 description="This script can co-add field maps corresponding "
                             "to different time intervals (a week, a month, "
                             "last n days, etc.) and make figures for the "
                             "resultant co-added maps.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                 epilog="Hopefully this works well! (Wei)")
    
    parser.add_argument("-s", "--season",
                        type=str, action="store",
                        choices=["winter", "summer"], default="winter",
                        help="If we want to analyze the maps from "
                             "the 1500 square degrees field, "
                             "we should choose 'winter'. If we want to anayze "
                             "the summer fields, we should choose 'summer'."
                             "This choice affects some of the "
                             "hardcoded variables used in the script.")
    
    parser.add_argument("-m", "--mode",
                        type=str, action="store",
                        choices=["coadding", "plotting"], default="plotting",
                        help="For each of the two seasons, "
                             "two modes are available: coadding and plotting. "
                             "In the former mode, this script will run another "
                             "script that coadds maps. In the latter mode, "
                             "a different script that makes figures for maps "
                             "and related quantities will be run.")
    
    parser.add_argument("-a", "--action",
                        type=str, action="store",
                        choices=["rebuild", "update"], default="update",
                        help="For each of the two modes mentioned above, "
                             "two actions are available: update and rebuild. "
                             "In the former case, the script will "
                             "make co-adds or generate figures "
                             "all over again for all relevant time intervals "
                             "specified by the argument time_interval below. "
                             "In the latter case, the script will "
                             "just update the co-adds/figures for those "
                             "time intervals that contain new information.")
    
    parser.add_argument("-o", "--oldest-time-to-consider",
                        type=str, action="store", default="20190321",
                        help="The oldest day to consider "
                             "when making coadds/figures. In other words, "
                             "data taken before this date will be ignored. "
                             "The format needs to match that of the default.")
    
    parser.add_argument("-c", "--current-time",
                        type=str, action="store", default=None,
                        help="The 'current' time when the script is run. "
                             "Having the ability to arbitrarily define "
                             "a 'current' time may be useful in testing. "
                             "The format of the string needs to match that of "
                             "the default of the previous argument.")
    
    parser.add_argument("-t", "--time-interval",
                        type=str, action="store",
                        choices=["last_n", "weekly", "monthly", "yearly"],
                        default="last_n",
                        help="Time intervals for which co-added maps "
                             "will be made or figures will be generated. "
                             "last_n refers to the last N days, "
                             "where N will be specified by the argument "
                             "last_how_many_days defined below.")
    
    parser.add_argument("-n", "--last-how-many-days",
                        type=int, action="store", default=1,
                        help="Number that specifies how long the time interval "
                             "last_n mentioned above corresponds to.")
    
    parser.add_argument("-d", "--original-maps-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory that contains "
                             "field maps from individual observations. It is "
                             "expected that this directory contains "
                             "four directories named "
                             "ra0hdec-44.75, ra0hdec-52.25, "
                             "ra0hdec-59.75, and ra0hdec-67.25, and "
                             "each sub-directory in turn contains g3 files "
                             "named like 12345678_abc.g3, which stores "
                             "temperature maps and weight maps.")
    
    parser.add_argument("-C", "--calibration-data-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory that contains "
                             "auto-processed calibration results. The results "
                             "from calibrator observations will be used when "
                             "calculating how much detectors' response changed "
                             "at the top of a field compared to the bottom.")
    
    parser.add_argument("-T", "--bolo-timestreams-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory that contains "
                             "raw timestreams. The elevation timestreams of "
                             "calibrator observations will be used to make sure "
                             "the observations taken at the correct elevations "
                             "are used.")
    
    parser.add_argument("-D", "--coadds-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory where co-added maps "
                             "are to be stored. This directory will in turn "
                             "contain four directories named last_n, weekly, "
                             "monthly, and yearly, each of which then contain "
                             "sub-directories corresponding to different "
                             "weeks, months, and so on. In each of these, "
                             "g3 file storing co-added maps and weight maps "
                             "will be saved.")
    
    parser.add_argument("-F", "--figs-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory where figures related "
                             "to co-added maps are to be stored. The directory "
                             "structure is the same as the one used for the "
                             "coadds-dir.")
    
    parser.add_argument("-i", "--bands",
                        type=str, action="store",
                        nargs="+", default=["90GHz", "150GHz", "220GHz"],
                        help="Relevant bands to take actions on. This option "
                             "is here mainly for debugging purposes.")
    
    parser.add_argument("-f", "--sub_fields",
                        type=str, action="store",
                        nargs="+", default=["ra0hdec-44.75", "ra0hdec-52.25",
                                            "ra0hdec-59.75", "ra0hdec-67.25"],
                        help="Relevant sub-fields to take actions on. This "
                             "option is here mainly for debugging purposes.")
    
    parser.add_argument("-j", "--just-see-commands",
                        action="store_true", default=False,
                        help="Whether the commands to be run that add maps or "
                             "make figures will just be shown instead of "
                             "actually being run. This is also for debugging.")
    
    parser.add_argument("-L", "--logger-name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "
                             "record log messages.")
    
    parser.add_argument("-l", "--log-file",
                        type=str, action="store", default=None,
                        help="The file to which the logs will be recorded.")
    
    arguments = parser.parse_args()
    
    update(**vars(arguments))


# ==============================================================================

