# ******* Always need to import modules! ******* #

import os
import sys
import shutil
import json
import argparse
import datetime
import logging

import matplotlib
matplotlib.use("Agg")
import numpy

from glob  import glob
from spt3g import core
from spt3g import std_processing

from summaryplot import fields_coadding
from summaryplot import fields_plotting



# ==============================================================================
# Define a function that calls the map coadding and plotting scripts
# with appropriate arguments
# ------------------------------------------------------------------------------


def update(mode, action, oldest_time_to_consider=None, current_time=None,
           time_interval=None, last_how_many_days=0, original_maps_dir='.',
           coadds_dir='.', figs_dir='.',
           map_ids=["90GHz", "150GHz", "220GHz"],
           sub_fields=["ra0hdec-44.75", "ra0hdec-52.25",
                       "ra0hdec-59.75", "ra0hdec-67.25"],
           just_see_commands=False, logger_name='', log_file=None):
    
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
    
    
    # -- Other variables.
    
    if current_time is None:
        current_time = datetime.datetime.utcnow()
    else:
        current_time = datetime.datetime.strptime(
                           current_time, "20%y%m%d")
    beginning_of_the_record = \
        datetime.datetime.strptime(oldest_time_to_consider, "20%y%m%d")

    desired_obs_id_ranges   = []
    desired_time_ranges     = []
    desired_dir_names       = []

    script_coadding_maps = "summaryplot/fields_coadding.py"
    script_plotting_data = "summaryplot/fields_plotting.py"
    point_source_file    = "spt3g_software/sources/1500d_ptsrc_3band_50mJy.txt"
    bad_map_list_file    = "summaryplot/fields_bad_map_list.txt"

    if original_maps_dir[-1] != "/":
        original_maps_dir += "/"
    if coadds_dir[-1] != "/":
        coadds_dir += "/"
    if figs_dir[-1] != "/":
        figs_dir += "/"
    
    
    
    # - Figure out what appropriate time intervals are and
    #   make sure an appropriate directory structure exists.
    
    # -- Figure out appropriate time intervals
    
    def convert_to_obs_id(datetime_object):
        return std_processing.time_to_obsid(
                   datetime_object.strftime("20%y%m%d_%H%M%S"))

    if time_interval == "last_n":

        delta_t      = datetime.timedelta(days=-1*last_how_many_days)
        end_time     = current_time
        start_time   = current_time + delta_t
        if start_time < beginning_of_the_record:
            start_time = beginning_of_the_record
        """start_time   = start_time.replace(hour=0, minute=0,
                                          second=0, microsecond=0)"""
        end_obs_id   = convert_to_obs_id(end_time)
        start_obs_id = convert_to_obs_id(start_time)

        desired_obs_id_ranges.append((start_obs_id, end_obs_id))
        desired_time_ranges.append((start_time, end_time))

    else:
        def get_the_beginning_of_the_interval(interval, end_time):
            if interval == "yearly":
                return end_time.replace(month=1, day=1, hour=0,
                                        minute=0, second=0, microsecond=0)
            elif interval == "monthly":
                return end_time.replace(day=1, hour=0,
                                        minute=0, second=0, microsecond=0)
            elif interval == "weekly":
                current_weekday = end_time.weekday()
                delta_t         = datetime.timedelta(days=-1*current_weekday)
                return (end_time + delta_t).replace(hour=0, minute=0,
                                                    second=0, microsecond=0)

        if mode == "plotting":
            if time_interval == "yearly":
                end_time_at_start_of_loop = current_time.replace(microsecond=0)
            elif time_interval == "monthly":
                beginning_of_next_month = \
                    current_time.replace(
                        month=current_time.month+1, day=1,
                        hour=0, minute=0, second=0, microsecond=0)
                end_of_this_month = \
                    beginning_of_next_month + \
                    datetime.timedelta(seconds=-1)
                end_time_at_start_of_loop = \
                    end_of_this_month
            elif time_interval == "weekly":
                end_day_of_this_week = \
                    current_time + \
                    datetime.timedelta(days=(6-current_time.weekday()))
                end_time_this_week = \
                    end_day_of_this_week.replace(
                        hour=23, minute=59, second=59, microsecond=0)
                end_time_at_start_of_loop = \
                    end_time_this_week
        elif mode == "coadding":
            end_time_at_start_of_loop = current_time

        while end_time_at_start_of_loop > beginning_of_the_record:
            beginning_of_the_interval = \
                get_the_beginning_of_the_interval(time_interval,
                                                  end_time_at_start_of_loop)
            if (mode == "coadding") and \
               (beginning_of_the_interval < beginning_of_the_record):
                start_obs_id_of_the_interval = \
                    convert_to_obs_id(beginning_of_the_record)
            else:
                start_obs_id_of_the_interval = \
                    convert_to_obs_id(beginning_of_the_interval)
            end_obs_id_of_the_interval = \
                convert_to_obs_id(end_time_at_start_of_loop)

            desired_obs_id_ranges.append((start_obs_id_of_the_interval,
                                          end_obs_id_of_the_interval))
            desired_time_ranges.append((beginning_of_the_interval,
                                        end_time_at_start_of_loop))

            end_time_at_start_of_loop = beginning_of_the_interval + \
                                            datetime.timedelta(seconds=-1)


    # -- Check the input and output directory structure
    
    log("")
    log("-------------------------------------------------")
    log(" Making sure the directory structure is valid... ")
    log("-------------------------------------------------")
    log("")
    
    log("* Checking whether relevant directories exist or not.")
    log("  If not, they will be created.")
    
    
    def convert_time_intervals_to_dir_names(interval_type, datetime_obj):
        if interval_type == "weekly":
            str_fmt = "20%y%m%d"
        elif interval_type == "monthly":
            str_fmt = "20%y%m"
        elif interval_type == "yearly":
            str_fmt = "20%y"
        return datetime_obj.strftime(str_fmt)


    for root_directory in [coadds_dir, figs_dir]:
        if time_interval not in os.listdir(root_directory):
            os.mkdir(root_directory+time_interval)

        for time_range in desired_time_ranges:
            if time_interval == "last_n":
                subdir_name = time_interval[:-1] + \
                              str(last_how_many_days) + "/"
            else:
                subdir_name = convert_time_intervals_to_dir_names(
                                  time_interval, time_range[0]) + "/"
            fulldir_name = root_directory + \
                           time_interval + "/" + \
                           subdir_name
            desired_dir_names.append(fulldir_name)


    for dir_name in desired_dir_names:
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)

    log("")
    

    
    # - Call the relevant functions with appropriate arguments
    
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
    for counter, interval in enumerate(desired_obs_id_ranges, 0):
        log("Interval %s :", counter+1)
        log("  from %s ("+\
            str(std_processing.obsid_to_g3time(interval[0])).split(".")[0]+")"+\
            "("+str(desired_time_ranges[counter][0])+")", interval[0])
        log("  to   %s ("+\
            str(std_processing.obsid_to_g3time(interval[1])).split(".")[0]+")"+\
            "("+str(desired_time_ranges[counter][1])+")", interval[1])
    log("")

    log("")
    log("--------------------------")
    log("Relevant commands to run: ")
    log("--------------------------")
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

        if max_eid < min_gid:
            ignore_coadds  = True
            ids_to_exclude = []
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid

        elif (min_eid<min_gid) and (max_eid>min_gid) and (max_eid<max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids \
                 if obs_id >= min_gid]
            subtract_maps  = True
            id_lower_bound = min_eid
            id_upper_bound = max_gid

        elif (min_eid < min_gid) and (max_eid > max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids \
                 if (obs_id >= min_gid) and (obs_id <= max_gid)]
            subtract_maps  = True
            id_lower_bound = min_eid
            id_upper_bound = max_eid

        elif (min_eid >= min_gid) and (max_eid <= max_gid):
            ignore_coadds  = False
            ids_to_exclude = existing_ids
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid

        elif (min_eid>=min_gid) and (min_eid<max_gid) and (max_eid>max_gid):
            ignore_coadds  = False
            ids_to_exclude = \
                [obs_id for obs_id in existing_ids \
                 if obs_id <= max_gid]
            subtract_maps  = True
            id_lower_bound = min_gid
            id_upper_bound = max_eid

        elif min_eid > max_gid:
            ignore_coadds  = True
            ids_to_exclude = []
            subtract_maps  = False
            id_lower_bound = min_gid
            id_upper_bound = max_gid

        else:
            raise RuntimeError("There seems to be a forgotten scenario!")

        return ignore_coadds,  subtract_maps, \
               ids_to_exclude, int(id_lower_bound), int(id_upper_bound)
    
    
    
    def figure_out_arguments_to_use_for_coadding(
            map_id, sub_field, time_range_id):
        
        arguments = {}
        
        g3_files_for_individual_observations = \
            sorted(glob(os.path.join(original_maps_dir,
                                     sub_field,
                                     '*{}_tonly.g3.gz'.format(map_id))))
        if sub_field == '*':
            g3_file_for_coadded_maps = \
                os.path.join(desired_dir_names[time_range_id],
                             'coadded_maps_{}.g3'.format(map_id))
        else:
            g3_file_for_coadded_maps = \
                os.path.join(desired_dir_names[time_range_id],
                             'coadded_maps_from_{}_{}.g3'.format(sub_field,
                                                                 map_id))
        
        arguments['output_file'] = '{}.g4'.format(g3_file_for_coadded_maps[:-3])
        arguments['map_ids'] = [map_id]
        
        if (action == 'update') and \
           (os.path.isfile(g3_file_for_coadded_maps)):
            iterator = core.G3File(g3_file_for_coadded_maps)
            existing_ids = []
            while True:
                try:
                    frame = iterator.next()
                    if "CoaddedObservationIDs" in frame.keys():
                        for sub_field, obs_ids \
                        in  frame["CoaddedObservationIDs"].items():
                            for obs_id in obs_ids:
                                if obs_id not in existing_ids:
                                    existing_ids.append(obs_id)
                        break
                except StopIteration:
                    break
        
            ignore_coadds,  subtract_maps, \
            ids_to_exclude, id_lower_bound, id_upper_bound = \
                decide_what_ids_to_use_and_not(
                    existing_ids,
                    desired_obs_id_ranges[time_range_id],
                    g3_file_for_coadded_maps)
            
            if ignore_coadds:
                arguments['input_files'] = g3_files_for_individual_observations
            else:
                arguments['input_files'] = [g3_file_for_coadded_maps] + \
                                           g3_files_for_individual_observations
            
            arguments['min_obs_id']  = id_lower_bound
            arguments['max_obs_id']  = id_upper_bound
            arguments['bad_obs_ids'] = ids_to_exclude
            
            if subtract_maps:
                arguments['subtract_existing_maps'] = True
                
        else:
            arguments['input_files'] = g3_files_for_individual_observations
            arguments['min_obs_id'] = desired_obs_id_ranges[i][0]
            arguments['max_obs_id'] = desired_obs_id_ranges[i][1]
                    
        arguments['temperature_maps_only'] = True
        arguments['point_source_file'] = point_source_file
        arguments['bad_map_list_file'] = bad_map_list_file
        arguments['log_file'] = log_file
        
        return arguments

    
    
    def run_command(function, arguments, logger, log_file, just_see_args):
        
        log('Calling the main function ...')
        log('Here are the arguments to be supplied:')
        arguments_to_show = {k: v for k, v in arguments.items()}
        for key in ['input_files', 'bad_obs_ids']:
            if key in arguments_to_show.keys():
                if len(arguments_to_show[key]) > 5:
                    arguments_to_show.pop(key)
        log(json.dumps(arguments_to_show, indent=1))
        log('')

        if not just_see_args:
            if log_file is None:
                function(**arguments)
            else:
                try:
                    function(**arguments)
                except Exception:
                    logger.exception('Something did not go well '
                                     'while running coadding/plotting!')
                    raise RuntimeError('An error occurred! '
                                       'Please check where it occurred '
                                       'in the log file {}!'.format(log_file))
        log('')
        
    
    
    def generate_new_coadded_maps(
            function, arguments, logger, log_file,
            just_see_commands, back_up_good_coadds):
        
        output_file_temp = arguments['output_file']
        output_file_perm = '{}.g3'.format(output_file_temp[:-3])
        
        run_command(function, arguments, logger, log_file, just_see_commands)
        
        log('Changing the name of the temporary output file ...')
        if os.path.isfile(output_file_temp):
            log('  (Original name: %s' , output_file_temp)
            log('   New name     : %s)', output_file_perm)
            if not just_see_commands:
                shutil.move(output_file_temp, output_file_perm)
            log('Done.')
        else:
            log('Actually, %s does not exist.', output_file_temp)
        log('')
        
        if back_up_good_coadds and os.path.isfile(output_file_perm):
            make_backup = "yes"
            backup      = '{}_backup.g3'.format(output_file_perm[:-3])
            iterator    = core.G3File(output_file_perm)
            new_frame   = iterator.next()
            obs_info    = new_frame["CoaddedObservationIDs"]
            noise_info  = new_frame["NoiseFromCoaddedMaps"]
            map_id      = new_frame["Id"]
            all_bad_ids = {}
            for sub_field, obs_ids in obs_info.items():
                n_over_time = \
                    numpy.array([noise_info[sub_field][str(obs_id)] \
                                 for obs_id in obs_ids])
                if not numpy.isfinite(n_over_time[-1]):
                    make_backup = "unsure"
                    break
                else:
                    valid_indices = numpy.where(numpy.isfinite(n_over_time))[0]
                    valid_data    = n_over_time[valid_indices]
                    if valid_data[-1] > valid_data[-2]:
                        bad_obs_ids = obs_ids[int(valid_indices[-2])+1:]
                        for boid in bad_obs_ids:
                            all_bad_ids[boid] = sub_field
                        make_backup = "no"
                        break
            
            if make_backup == "unsure":
                pass
            elif make_backup == "yes":
                log('Creating a backup for %s ...', output_file_perm)
                log('(b/c the data seem fine and may be useful in the future)')
                log('  (Original file: %s' , output_file_perm)
                log('   Backup file  : %s)', backup)
                if not just_see_commands:
                    shutil.copy(output_file_perm, backup)
                log('Done.')
                log('')
            else:
                log('The output file seems to contain problematic data')
                log('(noise in the running noise map '
                    'does not decrease monotonically).')
                log('The possibly bad observations are being recorded ...')
                for obs_id, sub_field in all_bad_ids.items():
                    fields_coadding.record_bad_obs_id(
                        bad_ids_file, map_id, sub_field, obs_id, "Noisy")
                log('Done.')
                log('')
                
                if os.path.isfile(backup):
                    log('Since fortunately a backup file exists,')
                    log('the ouput file will be replaced with this one,')
                    log('and coadding will be redone without the')
                    log('bad observations.')
                    log('Replacing the file ...')
                    if not just_see_commands:
                        os.remove(output_file_perm)
                        shutil.copy(backup, output_file_perm)
                else:
                    log('Since unfortunately there is no backup file,')
                    log('the output file will just be deleted,')
                    log('and the coadding will be redone,')
                    log('but without the bad observations.')
                    log('Deleting the file ...')
                    if not just_see_commands:
                        os.remove(output_file_perm)
                log('Done.')
                log('')
    
    
    # -- Then, execute commands by utilizing those functions
    
    n_time_ranges = len(desired_obs_id_ranges)
    
    for i in range(n_time_ranges):
        
        n_cmd_set = "{:03}".format(i+1)

        log("# Command set %s to run:", n_cmd_set)
        log("# ---------------------------")
        log("")
        
        sub_logger_name = '{}_{}_{}'.format(mode, time_interval, desired_dir_names[i].split('/')[-2])
        
        if mode == "coadding":
            if time_interval != "yearly":
                for map_id in map_ids:
                    log('--- %s ---', map_id)
                    log('')
                    
                    coadd_args = figure_out_arguments_to_use_for_coadding(map_id, '*', i)
                    coadd_args['collect_averages_from_flagging_statistics'] = True
                    coadd_args['calculate_pW_to_K_conversion_factors'] = True
                    coadd_args['calculate_pointing_discrepancies'] = True
                    coadd_args['calculate_noise_from_individual_maps'] = True
                    coadd_args['logger_name'] = '{}_{}'.format(sub_logger_name, map_id)

                    generate_new_coadded_maps(fields_coadding.run, coadd_args, logger, log_file, just_see_commands, False)
                    log('\n\n')
            else:
                for map_id in map_ids:
                    coadd_all_fields_args = {}
                    for sub_field in sub_fields:
                        log('--- %s %s ---', map_id, sub_field)
                        log('')
                        
                        coadd_args = figure_out_arguments_to_use_for_coadding(map_id, sub_field, i)
                        coadd_args.pop('subtract_existing_maps', None)
                        coadd_args['sources'] = [sub_field]
                        coadd_args['calculate_noise_from_coadded_maps'] = True
                        coadd_args['calculate_cross_spectrum_with_coadded_maps'] = True
                        coadd_args['logger_name'] = '{}_{}_{}'.format(sub_logger_name, map_id, sub_field.replace('.', ''))
                        """coadd_args = fields_coadding.gather_bad_obs_ids_from_list(
                                         bad_ids_file, map_id, sub_field, coadd_args)"""
                        
                        generate_new_coadded_maps(fields_coadding.run, coadd_args, logger, log_file, just_see_commands, False)
                        log('\n\n')
                    
                    log('--- %s Full field ---', map_id)
                    log('')
                    
                    coadd_all_fields_args = {}
                    coadd_all_fields_args['input_files'] = \
                        [os.path.join(desired_dir_names[i],
                                      'coadded_maps_from_{}_{}.g3'.format(sf, map_id)) for sf in sub_fields]
                    coadd_all_fields_args['output_file'] = \
                        os.path.join(desired_dir_names[i],
                                     'coadded_maps_{}.g3'.format(map_id))
                    coadd_all_fields_args['map_ids'] = [map_id]
                    coadd_all_fields_args['temperature_maps_only'] = True
                    coadd_all_fields_args['calculate_noise_from_coadded_maps'] = True
                    coadd_all_fields_args['calculate_cross_spectrum_with_coadded_maps'] = True
                    coadd_all_fields_args['logger_name'] = '{}_{}_full_field'.format(sub_logger_name, map_id)
                    coadd_all_fields_args['log_file'] = log_file
                    
                    run_command(fields_coadding.run, coadd_all_fields_args, logger, log_file, just_see_commands)
                    log('\n\n')

        elif mode == 'plotting':    
            for map_id in map_ids:
                log('--- %s ---', map_id)
                log('')
                
                plotting_args = {'input_files': [os.path.join(desired_dir_names[i],
                                                              'coadded_maps_{}.g3'.format(map_id))],
                                 'directory_to_save_figures': desired_dir_names[i+n_time_ranges],
                                 'simpler_file_names': True,
                                 'figure_title_font_size': 15,
                                 'map_id': map_id,
                                 'map_type': 'T',
                                 'coadded_data': True,
                                 'make_figures_for_field_maps': True,
                                 'smooth_map_with_gaussian': True,
                                 'gaussian_fwhm': 1,
                                 'make_figures_for_noise_levels': True,
                                 'log_file': log_file,
                                 'logger_name': '{}_{}'.format(sub_logger_name, map_id)}

                if time_interval != 'yearly':
                    plotting_args['make_figures_for_flagging_statistics'] = True
                    plotting_args['make_figures_for_pointing_discrepancies'] = True
                    plotting_args['left_xlimit_for_time_variations'] = desired_obs_id_ranges[i][0]
                    plotting_args['right_xlimit_for_time_variations'] = desired_obs_id_ranges[i][1]
                else:
                    plotting_args['make_figures_for_entire_weight_maps'] = True
                    plotting_args['make_figures_for_weight_maps_cross_section'] = True

                if action == 'update':
                    plotting_args['decide_whether_to_make_figures_at_all'] = True

                run_command(fields_plotting.run, plotting_args, logger, log_file, just_see_commands)
                log('\n\n')
        
        log('\n\n\n')


# ==============================================================================




# ==============================================================================
# Run the script from command line if desired
# ------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                 description="This script can co-add field maps corresponding to "
                             "different time intervals (a week, a month, etc.) "
                             "and make figures for the resultant co-added maps.",
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                 epilog="Hopefully this works well! (Wei)")

    parser.add_argument("-m", "--mode",
                        type=str, action="store",
                        choices=["coadding", "plotting"], default="plotting",
                        help="Two modes are available: coadding and plotting. "
                             "In the former mode, this script will run another "
                             "script that coadds maps. In the latter mode, "
                             "a different script that makes figures will be run.")

    parser.add_argument("-a", "--action",
                        type=str, action="store",
                        choices=["rebuild", "update"], default="update",
                        help="For each of the two modes mentioned above, "
                             "two actions are available: update and rebuild. "
                             "In the former case, the script will "
                             "make co-adds and/or generate figures "
                             "all over again for all relevant time intervals "
                             "specified by the argument time-interval below. "
                             "In the latter case, the script will "
                             "just update the co-adds/figures corresponding to "
                             "the most recent time interval.")

    parser.add_argument("-o", "--oldest-time-to-consider",
                        type=str, action="store", default="20190218",
                        help="The oldest time to consider "
                             "when making coadds/figures. In other words, "
                             "data taken before this date will be ignored. "
                             "The format needs to match that of the default.")

    parser.add_argument("-c", "--current-time",
                        type=str, action="store", default=None,
                        help="The 'current' time when the script it run. "
                             "Having the ability to arbitrarily define "
                             "a 'current' time may be useful in testing. "
                             "The format of the string needs to match that of "
                             "the default of the previous argument.")

    parser.add_argument("-t", "--time-interval",
                        type=str, action="store",
                        choices=["last_n", "weekly", "monthly", "yearly"],
                        default="last_n",
                        help="Time intervals for which co-added maps "
                             "will be made and/or figures will be generated. "
                             "last_n refers to the last N days, where N will be "
                             "specified by the argument last_how_many_days below.")

    parser.add_argument("-n", "--last-how-many-days",
                        type=int, action="store", default=3,
                        help="Number that specifies how long the time interval "
                             "last_n mentioned above corresponds to.")

    parser.add_argument("-d", "--original-maps-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory that contains field maps "
                             "from individual observations. It is expected that "
                             "this directory contains four directories named "
                             "ra0hdec-44.75, ra0hdec-52.25, "
                             "ra0hdec-59.75, and ra0hdec-67.25, and "
                             "each sub-directory in turn contains g3 files "
                             "named like 12345678.g3, which stores "
                             "temperature maps and weight maps.")

    parser.add_argument("-D", "--coadds-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory where co-added maps "
                             "are to be stored. This directory may already/will "
                             "contain four directories named last_n, weekly, "
                             "monthly, and yearly. The weekly and monthly ones "
                             "in turn contain multiple sub-directories "
                             "corresponding to different weeks and months. "
                             "Then, in each directory, a g3 file storing co-added "
                             "temperature maps and weight maps will be saved.")

    parser.add_argument("-F", "--figs-dir",
                        type=str, action="store", default=".",
                        help="The path to the directory where figures related to "
                             "co-added maps are to be stored. This directory "
                             "may already/will contain four directories named "
                             "last_n, weekly, monthly, and yearly. "
                             "The weekly and monthly ones "
                             "in turn contain multiple sub-directories "
                             "corresponding to different weeks and months. "
                             "Then, in each directory, a g3 file storing co-added "
                             "temperature maps and weight maps will be saved.")
    
    parser.add_argument("-i", "--map_ids",
                        type=str, action="store",
                        nargs="+", default=["90GHz", "150GHz", "220GHz"],
                        help="Relevant map IDs to take actions on. This option is "
                             "here mainly for debugging purposes.")
    
    parser.add_argument("-f", "--sub_fields",
                        type=str, action="store",
                        nargs="+", default=["ra0hdec-44.75", "ra0hdec-52.25",
                                            "ra0hdec-59.75", "ra0hdec-67.25"],
                        help="Relevant sub-fields to take actions on. This option is "
                             "here mainly for debugging purposes.")

    parser.add_argument("-j", "--just-see-commands",
                        action="store_true", default=False,
                        help="Whether the commands to be run will just be shown "
                             "instead of actually being run.")
    
    parser.add_argument("-L", "--logger-name",
                        type=str, action="store", default="",
                        help="The name of the logger that will be used to "
                             "record log messages.")
    
    parser.add_argument("-l", "--log-file",
                        type=str, action="store", default=None,
                        help="The file to which the logger will send messages.")

    arguments = parser.parse_args()
    
        
    update(mode=arguments.mode,
           action=arguments.action,
           oldest_time_to_consider=arguments.oldest_time_to_consider,
           current_time=arguments.current_time,
           time_interval=arguments.time_interval,
           last_how_many_days=arguments.last_how_many_days,
           original_maps_dir=arguments.original_maps_dir,
           coadds_dir=arguments.coadds_dir,
           figs_dir=arguments.figs_dir,
           map_ids=arguments.map_ids,
           sub_fields=arguments.sub_fields,
           just_see_commands=arguments.just_see_commands,
           logger_name=arguments.logger_name,
           log_file=arguments.log_file)


# ==============================================================================

