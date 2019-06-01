# ******* Always need to import modules! ******* #

import os
import sys
import argparse
import datetime
import numpy

from spt3g import core
from spt3g import std_processing



# ==============================================================================
# Understand what to do and define global variables for later use
# ------------------------------------------------------------------------------


# Parse arguments
# ------------------------------------------------------------------------------

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

parser.add_argument("-j", "--just-see-commands",
                    action="store_true", default=False,
                    help="Whether the commands to be run will just be shown "
                         "instead of actually being run.")

arguments = parser.parse_args()


# ------------------------------------------------------------------------------



# Define global variables
# ------------------------------------------------------------------------------

print()
print("# --------------------------------------------- #")
print("#  The script cache_field_maps.py was invoked!  #")
print("# --------------------------------------------- #")
print()

if arguments.current_time is None:
    current_time = datetime.datetime.utcnow()
else:
    current_time = datetime.datetime.strptime(
                       arguments.current_time, "20%y%m%d")
beginning_of_the_record = \
    datetime.datetime.strptime(arguments.oldest_time_to_consider, "20%y%m%d")

desired_obs_id_ranges   = []
desired_time_ranges     = []
desired_dir_names       = []

script_coadding_maps = "summaryplot/fields_coadding.py"
script_plotting_data = "summaryplot/fields_plotting.py"

sub_fields = ["ra0hdec-44.75", "ra0hdec-52.25",
              "ra0hdec-59.75", "ra0hdec-67.25"] 
map_ids    = ["90GHz", "150GHz", "220GHz"]
# map_ids    = ["90GHz"]

if arguments.original_maps_dir[-1] != "/":
    arguments.original_maps_dir += "/"
if arguments.coadds_dir[-1] != "/":
    arguments.coadds_dir += "/"
if arguments.figs_dir[-1] != "/":
    arguments.figs_dir += "/"


# ==============================================================================





# ==============================================================================
# Figure out what appropriate time intervals are and
# make sure an appropriate directory structure exists.
# ------------------------------------------------------------------------------


# List the time intervals of interest and the corresponding observation IDs
# ------------------------------------------------------------------------------

def convert_to_obs_id(datetime_object):
    return std_processing.time_to_obsid(
               datetime_object.strftime("20%y%m%d_%H%M%S"))


if arguments.time_interval == "last_n":

    delta_t      = datetime.timedelta(days=-1*arguments.last_how_many_days)
    end_time     = current_time
    start_time   = current_time + delta_t
    if start_time < beginning_of_the_record:
        start_time = beginning_of_the_record
    start_time   = start_time.replace(hour=0, minute=0,
                                      second=0, microsecond=0)
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
    
    if arguments.mode == "plotting":
        if arguments.time_interval == "yearly":
            end_time_at_start_of_loop = current_time.replace(microsecond=0)
        elif arguments.time_interval == "monthly":
            beginning_of_next_month = \
                current_time.replace(
                    month=current_time.month+1, day=1,
                    hour=0, minute=0, second=0, microsecond=0)
            end_of_this_month = \
                beginning_of_next_month + \
                datetime.timedelta(seconds=-1)
            end_time_at_start_of_loop = \
                end_of_this_month
        elif arguments.time_interval == "weekly":
            end_day_of_this_week = \
                current_time + \
                datetime.timedelta(days=(6-current_time.weekday()))
            end_time_this_week = \
                end_day_of_this_week.replace(
                    hour=23, minute=59, second=59, microsecond=0)
            end_time_at_start_of_loop = \
                end_time_this_week
    elif arguments.mode == "coadding":
        end_time_at_start_of_loop = current_time
    
    while end_time_at_start_of_loop > beginning_of_the_record:
        beginning_of_the_interval = \
            get_the_beginning_of_the_interval(arguments.time_interval,
                                              end_time_at_start_of_loop)
        if arguments.mode == "coadding":
            if beginning_of_the_interval < beginning_of_the_record:
                beginning_of_the_interval = beginning_of_the_record
        
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


"""if arguments.action == "update":
    desired_obs_id_ranges = [desired_obs_id_ranges[0]]
    desired_time_ranges   = [desired_time_ranges[0]]
elif arguments.action == "rebuild":
    desired_obs_id_ranges = desired_obs_id_ranges
    desired_time_ranges   = desired_time_ranges"""
# ** The update mode ahould check all time intervals, too,
# ** because maps may not appear on disk chronologically.


# ------------------------------------------------------------------------------



# Check the input and output directory structure
# ------------------------------------------------------------------------------

print()
print("-------------------------------------------------")
print(" Making sure the directory structure is valid... ")
print("-------------------------------------------------")
print()

"""
Deleting files seems unnecessary after all...

if (arguments.action == "rebuild") or (arguments.time_interval == "last_n"):
    if arguments.mode == "coadding":
        file_extension = "g3"
    elif arguments.mode == "plotting":
        file_extension = "png"
    
    print("* Since the mode is 'rebuild', or the time interval is 'last_n',")
    print("relevant", file_extension, "files will be deleted!", "\n")
    
    find_cmd = "find " + arguments.coadds_and_figs_dir + \
               arguments.time_interval + " -name '*'"  + \
               file_extension + " -type f"
    del_cmd  = find_cmd + " -delete"
    
    os.system(find_cmd)
    os.system(del_cmd)
"""

print("* Checking whether relevant directories exist or not.")
print("  If not, they will be created.")


def convert_time_intervals_to_dir_names(interval_type, datetime_obj):
    if interval_type == "weekly":
        str_fmt = "20%y%m%d"
    elif interval_type == "monthly":
        str_fmt = "20%y%m"
    elif interval_type == "yearly":
        str_fmt = "20%y"
    return datetime_obj.strftime(str_fmt)


for root_directory in [arguments.coadds_dir, arguments.figs_dir]:
    if arguments.time_interval not in os.listdir(root_directory):
        os.mkdir(root_directory+arguments.time_interval)
    
    for time_range in desired_time_ranges:
        if arguments.time_interval == "last_n":
            subdir_name = arguments.time_interval[:-1] + \
                          str(arguments.last_how_many_days) + "/"
        else:
            subdir_name = convert_time_intervals_to_dir_names(
                              arguments.time_interval, time_range[0]) + "/"
        fulldir_name = root_directory + \
                       arguments.time_interval + "/" + \
                       subdir_name
        desired_dir_names.append(fulldir_name)
        
        
for dir_name in desired_dir_names:
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


print()


# ==============================================================================





# ==============================================================================
# Run the relevant script with appropriate arguments
# ------------------------------------------------------------------------------

print()
print("----------------------------------------------------------")
print(" The mode, action, and time interval that were specified: ")
print("----------------------------------------------------------")
print()
print("  "+arguments.mode+", "+arguments.action+", "+arguments.time_interval)
print()


print()
print("----------------------------------------------")
print(" Relevant obs_id range(s) to take actions on :")
print("----------------------------------------------")
print()

for counter, interval in enumerate(desired_obs_id_ranges, 0):
    print("Interval", counter+1, ":")
    print("  from", interval[0], "("+\
          str(std_processing.obsid_to_g3time(interval[0])).split(".")[0]+")",
          "("+str(desired_time_ranges[counter][0])+")")
    print("  to  ", interval[1], "("+\
          str(std_processing.obsid_to_g3time(interval[1])).split(".")[0]+")",
          "("+str(desired_time_ranges[counter][1])+")")
print()


print("\n")
print("-------------------------")
print("Relevant commands to run:")
print("-------------------------")
print()


def decide_what_ids_to_use_and_not(existing_ids, desired_id_range):

    if len(existing_ids) == 0:
        raise RuntimeError("It appears that coadded maps were not made "
                           "successfully last time!")
    
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
    
    elif (min_eid < min_gid) and (max_eid > min_gid) and (max_eid < max_gid):
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
    
    elif (min_eid >= min_gid) and (min_eid < max_gid) and (max_eid > max_gid):
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
           ids_to_exclude, id_lower_bound, id_upper_bound

    
n_ranges = len(desired_obs_id_ranges)
for i in range(n_ranges):

    if arguments.mode == "coadding":
        script = script_coadding_maps
        
        if arguments.time_interval != "yearly":
            
            cmds = []
            
            for map_id in map_ids:
                g3_files_for_individual_observations = \
                    arguments.original_maps_dir + "ra0hdec*/" + \
                    "*" + map_id + "_tonly.g3.gz"
                g3_file_for_coadded_maps = \
                    desired_dir_names[i] + "coadded_maps_"+map_id+".g3"
                
                output_file = g3_file_for_coadded_maps[:-3] + ".g4"
                
                flags = []
                flags.append("-o " + output_file)
                flags.append("-i " + map_id)

                if (arguments.action == "update") and \
                   (os.path.isfile(g3_file_for_coadded_maps)):
                    map_frames   = core.G3File(g3_file_for_coadded_maps)
                    existing_ids = [] 
                    while True:
                        try:
                            map_frame = map_frames.next()
                            if "CoaddedObservationIDs" in map_frame.keys():
                                for sub_field, obs_ids \
                                in  map_frame["CoaddedObservationIDs"].items():
                                    for obs_id in obs_ids:
                                        if obs_id not in existing_ids:
                                            existing_ids.append(obs_id)
                                break
                        except StopIteration:
                            break
                    
                    ignore_coadds,  subtract_maps,
                    ids_to_exclude, id_lower_bound, id_upper_bound = \
                        decide_what_ids_to_use_and_not(
                            existing_ids, desired_obs_id_ranges[i])
                    
                    if ignore_coadds:
                        input_files = g3_files_for_individual_observations
                    else:
                        input_files = g3_file_for_coadded_maps + " " + \
                                      g3_files_for_individual_observations
                    
                    flags.append("-s " + str(id_lower_bound))
                    flags.append("-l " + str(id_upper_bound))
                    
                    if len(ids_to_exclude) > 0:
                        ids_to_exclude = \
                            " ".join([str(oid) \
                                      for oid in ids_to_exclude])
                        flags.append("-b " + ids_to_exclude)
                    
                    if subtract_maps:
                        flags.append("-u")
                                
                else:
                    input_files = g3_files_for_individual_observations
                    flags.append("-s " + str(desired_obs_id_ranges[i][0]))
                    flags.append("-l " + str(desired_obs_id_ranges[i][1]))
                
                flags.append("-t -a -p -n")
                flags = " ".join(flags)
                
                cmds.append(
                    "python" + " " + script + " " + input_files + " " + flags)
                if os.path.isfile(g3_file_for_coadded_maps):
                    cmds.append("cp" + " " + g3_file_for_coadded_maps + " " + \
                                g3_file_for_coadded_maps[:-3] + "_backup.g3")
                cmds.append("mv" + " " + output_file + " " + \
                            g3_file_for_coadded_maps)
        
        else:
            cmds = []
            
            for map_id in map_ids:        
                for sub_field in sub_fields:
                    g3_files_for_individual_observations = \
                        arguments.original_maps_dir + sub_field+ \
                        "/" + "*" + map_id + "_tonly.g3.gz"
                    g3_file_for_coadded_maps = \
                        desired_dir_names[i] + "coadded_maps_from_" + \
                        sub_field + "_" + map_id + ".g3"
                    
                    output_file = g3_file_for_coadded_maps[:-3] + ".g4"
                    
                    flags = []
                    flags.append("-o " + output_file)
                    flags.append("-i " + map_id)
                    flags.append("-S " + sub_field)

                    if (arguments.action == "update") and \
                       (os.path.isfile(g3_file_for_coadded_maps)):
                        map_frames   = core.G3File(g3_file_for_coadded_maps)
                        existing_ids = [] 
                        while True:
                            try:
                                mp_fr = map_frames.next()
                                if "CoaddedObservationIDs" in mp_fr.keys():
                                    for source, obs_ids \
                                    in  mp_fr["CoaddedObservationIDs"].items():
                                        if source != sub_field:
                                            continue
                                        for obs_id in obs_ids:
                                            if obs_id not in existing_ids:
                                                existing_ids.append(obs_id)
                                    break
                            except StopIteration:
                                break
                        
                        ignore_coadds,  subtract_maps,
                        ids_to_exclude, id_lower_bound, id_upper_bound = \
                            decide_what_ids_to_use_and_not(
                                existing_ids, desired_obs_id_ranges[i])
                        
                        if ignore_coadds:
                            input_files = g3_files_for_individual_observations
                        else:
                            input_files = g3_file_for_coadded_maps + " " + \
                                          g3_files_for_individual_observations
                        
                        flags.append("-s " + str(id_lower_bound))
                        flags.append("-l " + str(id_upper_bound))
                        
                        if len(ids_to_exclude) > 0:
                            ids_to_exclude = \
                                " ".join([str(oid) \
                                          for oid in ids_to_exclude])
                            flags.append("-b " + ids_to_exclude)
                        
                        if subtract_maps:
                            flags.append("-u")
                    
                    else:
                        input_files = g3_files_for_individual_observations
                        flags.append("-s " + str(desired_obs_id_ranges[i][0]))
                        flags.append("-l " + str(desired_obs_id_ranges[i][1]))
                    
                    
                    flags.append("-t -N")
                    flags = " ".join(flags)
                    
                    cmds.append(
                        "python"+" "+script + " " + input_files+" "+flags)
                    if os.path.isfile(g3_file_for_coadded_maps):
                        cmds.append("cp"+" "+ g3_file_for_coadded_maps + " "+\
                                    g3_file_for_coadded_maps[:-3]+"_backup.g3")
                    cmds.append(
                        "mv" + " " + output_file + " " + \
                         g3_file_for_coadded_maps)
                
                
                four_g3_files = \
                    " ".join([desired_dir_names[i] + "coadded_maps_from_" + \
                              sf + "_" + map_id + ".g3" for sf in sub_fields])
                final_out_file = \
                    desired_dir_names[i] + "coadded_maps_"+map_id+".g3"
                
                flags = []
                flags.append("-o " + final_out_file)
                flags.append("-i " + map_id)
                flags.append("-t -N")
                flags = " ".join(flags)
                cmd_final = \
                    "python" + " " + script + " " + four_g3_files + " " + flags
                cmds.append(cmd_final)
             
    
    else:
        script = script_plotting_data
        cmds   = []
        
        for map_id in map_ids:
            input_file = desired_dir_names[i] + "coadded_maps_" + map_id + ".g3"
            output_dir = desired_dir_names[i+n_ranges]
            
            flags = []
            flags.append("-d " + output_dir)
            flags.append("-i " + map_id)
            flags.append("-f -c -w -n")
            if arguments.time_interval != "yearly":
                flags.append("-F -p")
                flags.append("-l " + str(desired_obs_id_ranges[i][0]))
                flags.append("-r " + str(desired_obs_id_ranges[i][1]))
            else:
                flags.append("-N")
            flags.append("-z 15")
            flags = " ".join(flags)
            
            cmd = "python" + " " + script + " " + input_file + " " + flags
            cmds.append(cmd)
    
    
    if i+1 < 10:
        n_cmd_set = "00" + str(i+1)
    elif i+1 < 100:
        n_cmd_set = "0"  + str(i+1)
    else:
        n_cmd_set = str(i+1)
    
    print("# Command set", n_cmd_set, "to run:")
    print("# ---------------------------")
    print()
    for cmd in cmds:
        print(cmd, "\n")
        if not arguments.just_see_commands:
            os.system(cmd)
    print("\n")

print()


# ==============================================================================


