# ============================================================================ #
#  This script is intended to be called by spt3g_dq_system/update_summary.py.  #
#  When update_summary.py calls this script, options such as                   #
#  whether we want to process the data saved in the ARC files or               #
#  want to make plots from the processed data are specified.                   #
#                                                                              #
#  Then, based on the specifications, this script in turn calls                #
#  another script, which is either arc_reducing.py or arc_plotting.py, with    #
#  the appropriate arguments.                                                  #
#                                                                              #
# ============================================================================ #


import logging
import datetime
import os
import glob
import pickle

from spt3g       import std_processing
from summaryplot import arc_reducing
from summaryplot import arc_plotting



def run(mode=None,
        action=None,
        interval=None,
        oldest_time_to_consider=None,
        newest_time_to_consider=None,
        arc_files_dir=None,
        pickles_dir=None,
        figures_dir=None,
        just_see_commands=True,
        logger_name=None,
        log_file=None):    


    # - Set up logging

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

    logger.info("")
    logger.info("# ----------------------------------------------- #")
    logger.info("#  The script cache_archive_data.py was invoked!  #")
    logger.info("# ----------------------------------------------- #")
    logger.info("")


    # - Decide on the time intervals we care about

    if newest_time_to_consider is None:
        newest_time_to_consider = datetime.datetime.utcnow()
    else:
        newest_time_to_consider = datetime.datetime.strptime(
                                      newest_time_to_consider, '20%y%m%d')
    newest_time_to_consider = newest_time_to_consider.replace(
                                  hour=23, minute=59,
                                  second=59, microsecond=0)
    oldest_time_to_consider = datetime.datetime.strptime(
                                  oldest_time_to_consider, '20%y%m%d')
    oldest_time_to_consider = oldest_time_to_consider.replace(
                                  hour=0, minute=0,
                                  second=0, microsecond=0)

    def convert_datetime_obj_to_obs_id(datetime_obj):
        return std_processing.time_to_obsid(
                   datetime_obj.strftime('20%y%m%d_%H%M%S'))

    if mode == 'reducing':
        newest_obs_id = convert_datetime_obj_to_obs_id(newest_time_to_consider)
        oldest_obs_id = convert_datetime_obj_to_obs_id(oldest_time_to_consider)
        desired_obs_id_ranges = [(oldest_obs_id,
                                  newest_obs_id)]
        desired_time_ranges   = [(oldest_time_to_consider,
                                  newest_time_to_consider)]
    elif mode == 'plotting':
        desired_obs_id_ranges = []
        desired_time_ranges   = []
        desired_outdir_paths  = []
        one_sec_ago = datetime.timedelta(seconds=-1)

        if 'last' in interval:
            n = int(interval.split('_')[1])
            n_days_ago   = datetime.timedelta(days=-1*n)
            end_time     = newest_time_to_consider
            start_time   = newest_time_to_consider + n_days_ago
            start_time   = start_time.replace(hour=0, minute=0, second=0)
            end_obs_id   = convert_datetime_obj_to_obs_id(end_time)
            start_obs_id = convert_datetime_obj_to_obs_id(start_time)
            desired_obs_id_ranges.append([start_obs_id, end_obs_id])
            desired_time_ranges.append((start_time, end_time))
        else:
            if interval == 'weekly':
                current_week_day = newest_time_to_consider.weekday()
                n_days_later     = datetime.timedelta(days=6-current_week_day)
                end_of_this_week = newest_time_to_consider + n_days_later
                end_time_at_start_of_loop = end_of_this_week
            elif interval == 'monthly':
                this_year  = newest_time_to_consider.year
                this_month = newest_time_to_consider.month
                delta_year = 1 if this_month == 12 else 0
                new_month  = 1 if this_month == 12 else this_month + 1
                start_of_next_month = newest_time_to_consider.replace(
                                          year=this_year+delta_year,
                                          month=new_month, day=1,
                                          hour=0, minute=0, second=0)
                end_time_at_start_of_loop = start_of_next_month + one_sec_ago

            while end_time_at_start_of_loop > oldest_time_to_consider:
                if interval == 'weekly':
                    six_days_ago    = datetime.timedelta(days=-6)
                    start_of_a_week = end_time_at_start_of_loop + six_days_ago
                    start_of_this_interval = start_of_a_week.replace(
                                                 hour=0, minute=0, second=0)
                elif interval == 'monthly':
                    start_of_this_interval = end_time_at_start_of_loop.replace(
                                                 day=1, hour=0,
                                                 minute=0, second=0)
                end_obs_id   = convert_datetime_obj_to_obs_id(
                                   end_time_at_start_of_loop)
                start_obs_id = convert_datetime_obj_to_obs_id(
                                   start_of_this_interval)
                desired_obs_id_ranges.append(
                    (start_obs_id, end_obs_id))
                desired_time_ranges.append(
                    (start_of_this_interval, end_time_at_start_of_loop))

                end_time_at_start_of_loop = start_of_this_interval + one_sec_ago


    # - Check the existence of needed output directories

    if mode == 'reducing':
        if not os.path.isdir(pickles_dir):
            os.mkdir(pickles_dir)

    elif mode == 'plotting':
        if not os.path.isdir(figures_dir):
            os.mkdir(figures_dir) 
        if 'last' in interval:
            interval_root_dir = 'last_n'
        else:
            interval_root_dir = interval
        if interval_root_dir not in os.listdir(figures_dir):
            os.mkdir(os.path.join(figures_dir, interval_root_dir))

        for time_range in desired_time_ranges:
            if 'last' in interval:
                sub_dir = interval
            elif interval == 'weekly':
                sub_dir = time_range[0].strftime('20%y%m%d')
            elif interval == 'monthly':
                sub_dir = time_range[0].strftime('20%y%m')
            full_path = os.path.join(figures_dir, interval_root_dir, sub_dir)
            desired_outdir_paths.append(full_path)
        for path in desired_outdir_paths:
            if not os.path.isdir(path):
                os.mkdir(path)


    # - Make sure things were set up correctly

    logger.info('')
    logger.info('----------------------------------------------------------')
    logger.info(' The mode, action, and time interval that were specified: ')
    logger.info('----------------------------------------------------------')
    logger.info('')
    logger.info('   %s, %s, %s', mode, action, interval)
    logger.info('\n')
    logger.info('----------------------------------------------')
    logger.info(' Relevant obs_id range(s) to take actions on: ')
    logger.info('----------------------------------------------')
    logger.info('')
    for counter, interval in enumerate(desired_obs_id_ranges, 1):
        logger.info('Interval %s :', counter)
        logger.info(
            '  from %s (%s)',
            str(std_processing.obsid_to_g3time(interval[0])).split('.')[0],
            interval[0])
        logger.info(
            '  to   %s (%s)',
            str(std_processing.obsid_to_g3time(interval[1])).split(".")[0],
            interval[1])
    logger.info("")


    # - Actually run the data reduction or plotting script

    if mode == 'reducing':

        # -- Gather the appropriate input files
        
        def find_arc_files_from_some_time_range(arc_files_dir, t_start, t_end):
            all_files = glob.glob(os.path.join(arc_files_dir, '20*_*.dat'))
            all_files = sorted(all_files)
            for i, f in enumerate(all_files, 0):
                time_stamp   = os.path.basename(all_files[i]).split('.')[0]
                datetime_obj = datetime.datetime.strptime(
                                   time_stamp, '20%y%m%d_%H%M%S')
                if datetime_obj > t_start:
                    idx_start = i - 1
                    if idx_start == -1:
                        idx_start = 0
                    break
            for i, f in enumerate(all_files, 1):
                time_stamp   = os.path.basename(all_files[-i]).split('.')[0]
                datetime_obj = datetime.datetime.strptime(
                                   time_stamp, '20%y%m%d_%H%M%S')
                if datetime_obj < t_end:
                    idx_end = -i + 1
                    if idx_end >= -1:
                        idx_end = len(all_files) - 1
                    break
            return set(all_files[idx_start:idx_end+1])

        start_time = desired_time_ranges[0][0]
        end_time   = desired_time_ranges[0][1]
        all_arcs_from_this_period = find_arc_files_from_some_time_range(
                                        arc_files_dir, start_time, end_time)
        path_to_reduction_record = os.path.join(
                                       pickles_dir,
                                       'files_processed_so_far.pickle')
        try:
            with open(path_to_reduction_record, 'rb') as fobj:
                files_already_used = pickle.load(fobj)
        except:
            files_already_used = set()

        if (action == 'rebuild'):
            input_files = all_arcs_from_this_period
        else:
            input_files = all_arcs_from_this_period - files_already_used
        input_files = sorted(list(input_files))


        # -- Run the data reduction process

        if just_see_commands:
            logger.info('')
            logger.info('----------------------------')
            logger.info(' These are the input files: ')
            logger.info('----------------------------')
            logger.info('')
            for counter, f in enumerate(input_files, 1):
                logger.info(' file %06d: %s', counter, f)
            logger.info('')
            logger.info('That is it because this is a dry run!')
        else:
            try:
                logger.info('')
                logger.info('Starting the data reduction process...')
                logger.info('')
                arc_reducing.run(input_files=input_files,
                                 pickles_dir=pickles_dir,
                                 logger=logger)
                logger.info('')
                logger.info('Recording the files processed this time...')
                logger.info('')
                with open(path_to_reduction_record, 'wb') as fobj:
                    files_already_used = files_already_used | set(input_files)
                    pickle.dump(set(files_already_used), fobj)
                logger.info('')
                logger.info('All done!')
            except:
                logger.info('')
                logger.exception('An error occurred:\n------------------')
                raise RuntimeError('An error occurred! '
                                   'Check {}!'.format(log_file))


    elif mode == 'plotting':

        def find_pickle_files_from_some_time_range(
                pickles_dir, start_obs_id, end_obs_id):
            relevant_pickles = []
            for some_pickle in sorted(os.listdir(pickles_dir)):
                try:
                    this_obs_id = int(some_pickle.split('.')[0])
                    if start_obs_id <= this_obs_id <= end_obs_id:
                        relevant_pickles.append(
                            os.path.join(pickles_dir, some_pickle))
                except:
                    ### Probably a file whose name is not 12345678.pickle
                    pass
            return set(relevant_pickles)

        for counter_r, obs_id_range in enumerate(desired_obs_id_ranges, 1):

            # -- Gather the appropriate input files

            start_obs_id = obs_id_range[0]
            end_obs_id   = obs_id_range[1]
            output_dir   = desired_outdir_paths[counter_r-1]
            input_files  = find_pickle_files_from_some_time_range(
                               pickles_dir, start_obs_id, end_obs_id)

            path_to_plotting_record = os.path.join(
                                          output_dir,
                                          'files_used_last_time.pickle')
            try:
                with open(path_to_plotting_record, 'rb') as fobj:
                    files_processed_last_time = pickles.load(fobj)
            except:
                files_processed_last_time = set()

            if action != 'rebuild':
                input_files = input_files - files_processed_last_time
            input_files = sorted(list(input_files))

            logger.info('')
            logger.info('# Starting the plotting process '
                        'for time interval %d!', counter_r)
            logger.info('')

            if just_see_commands:
                logger.info('------------------------------------')
                logger.info('Relevant arguments for interval %03d', counter_r)
                logger.info('------------------------------------')
                logger.info('')
                logger.info('Start time as in observation ID: %d', start_obs_id)
                logger.info('End   time as in observation ID: %d', end_obs_id)
                logger.info('Output directory:')
                logger.info(' %s', output_dir)
                logger.info('Input files:')
                if len(input_files) == 0:
                    logger.info(' Actually, all the files have'
                                ' already been used for plotting!')
                else:
                    for counter_f, f in enumerate(input_files, 1):
                        logger.info(' file %04d: %s', counter_f, f)
                logger.info('')
                logger.info('That is it because this is a dry run!')
                logger.info('')
            else:
                try:
                    if len(input_files) == 0:
                        logger.info('Actually, all the files have '
                                    'already been used for plotting')
                        logger.info('(or there are no files yet), '
                                    'so, there is nothing to do!')
                    else:
                        arc_plotting.run(input_files=input_files,
                                         output_dir=output_dir,
                                         start_obs_id=start_obs_id,
                                         end_obs_id=end_obs_id,
                                         logger=logger)
                        logger.info('')
                        logger.info('Recording the files used this time...')
                        with open(path_to_plotting_record, 'wb') as fobj:
                            pickle.dump(set(input_files), fobj)
                    logger.info('')
                except:
                    logger.info('')
                    logger.info('An error occurred:\n-----------------')
                    raise RuntimeError('An error occurred! '
                                       'Check {}'.format(log_file))

            logger.info('Finished the process '
                        'for time interval %d!', counter_r)
            logger.info('\n')
