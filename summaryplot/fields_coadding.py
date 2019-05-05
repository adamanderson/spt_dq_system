# ----  Always need to import modules!  ---- #

from    spt3g              import  core
from    spt3g              import  mapmaker
from    spt3g.mapspectra   import  map_analysis
import  os
import  gc
import  sys
import  glob
import  numpy
import  argparse
import  spt3g.util.filtered_gauss_fit  as  twod_gaussian_fit




# ==============================================================================
# Define global variables and gather input files
# ------------------------------------------------------------------------------


# Parse arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
             description="This script makes a coadded map of the "+\
                         "1500 sq. deg. field.")

parser.add_argument("input_files",
                    type=str, action="store", nargs="+",
                    help="The paths to the g3 files "+\
                         "that contain map frames. "+\
                         "Every file should be named like 12345678.g3, "
                         "except for a file that contains coadded maps.")

parser.add_argument("-o", "--output_file",
                    type=str, action="store", default="./coadded_map.g3",
                    help="The location and name of the output g3 file "+\
                         "that will contain the coadded map. The default "+\
                         "is ./coadded_map.g3.")

parser.add_argument("-i", "--map_ids",
                    type=str, action="store", nargs="+", default=["90GHz"],
                    help="The relevant map IDs for which data will be added. "+\
                         "The default is 90GHz.")

parser.add_argument("-c", "--combine_left_right",
                    action="store_true", default=False,
                    help="Whether to add maps made from left-going scans and "+\
                         "those made from right-going scans together. "+\
                         "If True, a map whose 'Id' is 'Left'+map_id and "+\
                         "one whose 'Id' is 'Right'+map_id will be added, "+\
                         "where map_id is an ID specified in the argument "+\
                         "map_ids above. The default is False.")

parser.add_argument("-C", "--combine_different_wafers",
                    action="store_true", default=False,
                    help="Whether to add maps made from different wafers "+\
                         "together. If True, maps whose 'Id's look like "+\
                         "map_id+a_wafer_name will be added, "+\
                         "where map_id is an ID specified in the argument "+\
                         "map_ids above. The default is False.")

parser.add_argument("-n", "--calculate_noise_from_individual_maps",
                    action="store_true", default=False,
                    help="Whether to calculate map noise from the difference "+\
                         "between a map made with only right-going scans and "+\
                         "a map made with only left-going scans for "+\
                         "each observation that will go into the coadded maps.")

parser.add_argument("-N", "--calculate_noise_from_coadded_maps",
                    action="store_true", default=False,
                    help="Whether to calculate map noise from the difference "+\
                         "between a map made with only right-going scans and "+\
                         "a map made with only left-going scans for "+\
                         "running coadded maps. In other words, the noise "+\
                         "will be calculated from coadded right-going maps "+\
                         "and left-going ones every time a new map is added.")

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
                         "in the coadded map. The default includes "+\
                         "all four fields.")

parser.add_argument("-s", "--min_obs_id",
                    type=int, action="store", default="00000000",
                    help="The smallest observation ID that will be "+\
                         "considered to be used in the coadded map. "+\
                         "The default is 00000000.")

parser.add_argument("-l", "--max_obs_id",
                    type=int, action="store", default="99999999",
                    help="The largest observation ID that will be "+\
                         "considered to be used in the coadded map. "+\
                         "The default is 99999999.")

parser.add_argument("-b", "--bad_obs_ids",
                    type=int, action="store", nargs="+", default=[],
                    help="The observation IDs that will be excluded from "+\
                         "making the coadded map. The script has some "+\
                         "simple criteria to decide whether to include "+\
                         "certain observations, but one can manually specify "+\
                         "bad observations, too. The default is None.")

parser.add_argument("-m", "--min_file_size",
                    type=float, action="store", default=0.01,
                    help="The minimum size (in GB) a g3 file needs "+\
                         "to have for it to be considered as a good file. "+\
                         "This is to reduce the occurrence of a situation "+\
                         "where the script crashes due to reading a "+\
                         "problematic file. The default is 0.01 GB (10 MB).")

arguments = parser.parse_args()



# -----------------------------------------------------------------------------
# Define some global variables
# -----------------------------------------------------------------------------

tqu_type = "T"

too_big_map_vals = {"T": {}}
for map_id in arguments.map_ids:
    ## Don't set thresholds for now
    if "90GHz" in map_id:
        too_big_map_vals["T"][map_id] = numpy.inf
    elif "150GHz" in map_id:
        too_big_map_vals["T"][map_id] = numpy.inf
    elif "220GHz" in map_id:
        too_big_map_vals["T"][map_id] = numpy.inf


if (arguments.calculate_noise_from_individual_maps or \
    arguments.calculate_noise_from_coadded_maps    or \
    arguments.calculate_pointing_discrepancies)   and \
   (not arguments.combine_left_right):
    for map_id in arguments.map_ids:
        if ("Left"  not in map_id) and \
           ("Right" not in map_id):
            print()
            print(map_id, "is an invalid ID because it doesn't contain")
            print("information on scan direction and because")
            print("noise is to be calculated from 'left-right' maps.")
            print("Also, currently, if pointing discrepancies")
            print("are to be calculated, then the pipeline needs to")
            print("receive a map made from left-going scans and")
            print("one made from right-going scans for a given obs..")
            print()
            sys.exit()

if arguments.calculate_noise_from_coadded_maps:
    if len(arguments.sources) > 1:
        print()
        print("If noise is to be calculated from running coadds,")
        print("then only one observation type can be specified!")
        print()
        sys.exit()
    if arguments.combine_left_right:
        print()
        print("If noise is to be calculated from running coadds,")
        print("then the combine_left_right option needs to be False")
        print("so that the coadded left-going and right-going maps")
        print("can be used later as well!")
        print()
        sys.exit()

if arguments.calculate_noise_from_individual_maps and \
   arguments.calculate_noise_from_coadded_maps:
    print()
    print("Currently, calculating noise from both individual maps")
    print("and coadded maps does not work...")
    print()
    sys.exit()



# -----------------------------------------------------------------------------
# Gather input files
# -----------------------------------------------------------------------------

def is_a_good_obs_id(g3_file):
    try:
        obs_id = int(g3_file.split("/")[-1].split(".")[0])
    except ValueError:
        # Probably a file that has coadded maps is encountered
        # if the name cannot be converted to an integer.
        return True
    if (obs_id >= arguments.min_obs_id) and \
       (obs_id <= arguments.max_obs_id) and \
       (obs_id not in arguments.bad_obs_ids):
        return True
    else:
        return False


all_good_g3_files = []
for g3_file in arguments.input_files:
    if os.path.isfile(g3_file) and \
       is_a_good_obs_id(g3_file) and \
       (os.path.getsize(g3_file) > arguments.min_file_size*2**30):
        all_good_g3_files.append(g3_file)


if (len(all_good_g3_files) == 0):
    print()
    print("No applicable input files, so nothing to do!")
    print("\n\n")
    sys.exit()



# -----------------------------------------------------------------------------
# Show settings related to I/O
# -----------------------------------------------------------------------------

print()
print("# =========================== #")
print("# Start making coadded maps! #")
print("# =========================== #")
print()

print("- Input files to check:")
for input_file in all_good_g3_files:
    print(input_file)
print()
print("- File to which the coadded map will be saved:")
print(arguments.output_file)
print()
print("- The range of observation ID of interest:")
print("from", arguments.min_obs_id, "to", arguments.max_obs_id)
print()
print("- Observation IDs to definitely exclude:")
print(arguments.bad_obs_ids)
print()


# ==============================================================================





# ==============================================================================
# Define modules and functions needed for the pipeline
# ------------------------------------------------------------------------------

def map_seems_fine(map_frame, tqu_type, dict_thresh, id_for_thresh):
    
    if map_frame[tqu_type].units != core.G3TimestreamUnits.Tcmb:
        print("The units of this map are not in Tcmb!")
        return False
    
    if tqu_type == "T":
        map_values = map_frame["T"]
    
    if not numpy.isfinite(numpy.nanmean(numpy.asarray(map_values))):
        print()
        print("There seem to be only NaNs in the map!")
        return False
    
    """
    ## Ignore this part for now
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



def combine_info_on_ids(old_dict, new_dict):
    
    this_frame_already_added = False
    for map_id in new_dict.keys():
        for source, new_ids in new_dict[map_id].items():
            if source not in old_dict[map_id].keys():
                old_dict[map_id][source] = new_ids
            else:
                common_ids = set(new_ids) & \
                             set(old_dict[map_id][source])
                if len(common_ids) > 0:
                    this_frame_already_added = True
                else:
                    old_dict[map_id][source] = \
                        type(new_ids)(sorted(
                            set(old_dict[map_id][source]) | set(new_ids)))
    
    return old_dict, this_frame_already_added



def calculate_noise_from_difference_map(rg_mp_fr, lg_mp_fr, tqu_type):
    
    df_mp  = map_analysis.subtract_two_maps(
                 rg_mp_fr, lg_mp_fr, divide_by_two=True)
    df_cls = map_analysis.calculateCls(
                 df_mp, apod_mask="default", qu=False,
                 kspace_filt=None, t_only=True,
                 ell_min=300, ell_max=5000)
    
    idx = numpy.where((df_cls["ell"]>3000) & (df_cls["ell"]<5000))[0]
    
    if tqu_type == "T":
        noise = numpy.sqrt(numpy.mean(df_cls["TT"][idx]))
    
    return noise



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



def calculate_pointing_discrepancies(rg_mp_fr, lg_mp_fr, tqu_type, sub_field):
    if tqu_type == "T":
        total_weight = rg_mp_fr["Wunpol"].TT + lg_mp_fr["Wunpol"].TT
        sum_map      = rg_mp_fr["T"] + lg_mp_fr["T"]
        map_to_use   = (sum_map / total_weight) / core.G3Units.mK

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
    discrep_dict = {"1": {"delta_ra": numpy.inf, "delta_dec": numpy.inf},
                    "2": {"delta_ra": numpy.inf, "delta_dec": numpy.inf},
                    "3": {"delta_ra": numpy.inf, "delta_dec": numpy.inf}}    
    
    info_on_point_sources = some_brightest_sources[sub_field]
    
    for point_source_number, info in info_on_point_sources.items():
        true_right_ascension = info[0]
        true_declination     = info[1]
        
        source_nominal_pixel_number = \
            map_to_use.angle_to_pixel(
                true_right_ascension*core.G3Units.deg,
                true_declination*core.G3Units.deg)
        source_y_index, source_x_index = \
            numpy.unravel_index(
                source_nominal_pixel_number,
                map_to_use.shape)
        
        mini_map_width    = 5 * core.G3Units.arcmin
        mini_map_height   = 5 * core.G3Units.arcmin
        mini_map_x_length = int(mini_map_width  / map_to_use.x_res)
        mini_map_y_length = int(mini_map_height / map_to_use.y_res)
        half_x_length     = mini_map_x_length // 2
        half_y_length     = mini_map_y_length // 2
        mini_map_x_left   = source_x_index - half_x_length
        mini_map_x_right  = source_x_index + half_x_length
        mini_map_y_bottom = source_y_index - half_y_length
        mini_map_y_top    = source_y_index + half_y_length
        mini_map = numpy.asarray(map_to_use)\
                       [mini_map_y_bottom:mini_map_y_top+1,
                        mini_map_x_left:mini_map_x_right+1]
        
        if False in numpy.isfinite(mini_map):
            discrep_dict[point_source_number]["delta_ra"]  = numpy.nan
            discrep_dict[point_source_number]["delta_dec"] = numpy.nan
            continue
        
        fit_parameters = twod_gaussian_fit.fitFiltered2Dgaussian(mini_map)
        gaussian_center_x_index = int(numpy.round(fit_parameters[2]))
        gaussian_center_y_index = int(numpy.round(fit_parameters[3]))
        
        measured_source_x_index = \
            gaussian_center_x_index - half_x_length + source_x_index
        measured_source_y_index = \
            gaussian_center_y_index - half_y_length + source_y_index
        
        measured_right_ascension, measured_declination = \
            map_to_use.pixel_to_angle(
                int(measured_source_x_index), int(measured_source_y_index))
        
        measured_right_ascension /= core.G3Units.deg
        measured_declination     /= core.G3Units.deg
        if measured_right_ascension < 0.0:
                measured_right_ascension += 360.0
        
        discrepancy_in_right_ascension = \
            (measured_right_ascension - true_right_ascension) * 3600
        discrepancy_in_right_ascension *= \
            numpy.cos(true_declination * numpy.pi / 180)
        
        discrepancy_in_declination = \
            (measured_declination - true_declination) * 3600
        
        discrep_dict[point_source_number]["delta_ra"]  = \
            discrepancy_in_right_ascension
        discrep_dict[point_source_number]["delta_dec"] = \
            discrepancy_in_declination
            
    return discrep_dict



class CoaddMapsAndCalculateNoises(object):
    
    def __init__(self, map_sources, map_ids, tqu_type,
                 too_big_map_vals,
                 combine_left_right=False,
                 combine_different_wafers=False,
                 calculate_noise_from_individual_maps=False,
                 calculate_noise_from_coadded_maps=False,
                 calculate_pointing_discrepancies=False):
        
        
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
        
        self.id_mapping = {}
        for map_id in map_ids:
            for direction in directions:
                for wafer in wafers:
                    self.id_mapping[direction+map_id+wafer] = map_id
        
        if combine_left_right or combine_different_wafers:
            print("- An ID stored in a map frame will be treated")
            print("  as a different ID according to the relations below:")
            for id_from in sorted(self.id_mapping.keys()):
                print("     ", id_from, "==>", self.id_mapping[id_from])
        
        self.map_ids     = map_ids
        self.map_sources = map_sources
        self.tqu_type    = tqu_type
        self.dict_thresh = too_big_map_vals
        
        
        # - Initialize variables storing maps and relevant info

        self.coadded_obs_ids = {map_id: core.G3MapVectorInt()    \
                                for map_id in self.map_ids}
        self.coadded_map_ids = {map_id: core.G3MapVectorString() \
                                for map_id in self.map_ids}
        
        self.coadded_map_frames = {map_id: None for map_id in self.map_ids}
        
        self.obs_frame = None
        
        
        # - Initialize variables related to noise and pointing calculations
        
        self.calculate_noise_from_individual_maps = \
            calculate_noise_from_individual_maps
        self.calculate_noise_from_coadded_maps    = \
            calculate_noise_from_coadded_maps
        self.calculate_pointing_discrepancies     = \
            calculate_pointing_discrepancies

        
        if self.calculate_noise_from_individual_maps or \
           self.calculate_noise_from_coadded_maps    or \
           self.calculate_pointing_discrepancies:
            self.direction_independent_map_ids = []
            self.map_frame_arrival_counters    = {}
            for map_id in self.map_ids:
                if "Left" in map_id:
                    self.direction_independent_map_ids.append(
                        map_id.replace("Left", "SomeDirection"))
                    self.map_frame_arrival_counters[map_id] = 0
                elif "Right" in map_id:
                    self.direction_independent_map_ids.append(
                        map_id.replace("Right", "SomeDirection"))
                    self.map_frame_arrival_counters[map_id] = 0
                else:
                    self.direction_independent_map_ids.append(
                        map_id)
                    self.map_frame_arrival_counters["Left"+map_id]  = 0
                    self.map_frame_arrival_counters["Right"+map_id] = 0
            self.direction_independent_map_ids = \
                set(self.direction_independent_map_ids)
        
        if self.calculate_noise_from_individual_maps or \
           self.calculate_pointing_discrepancies:
            self.map_frame_cache = \
                {map_id: None for map_id in self.map_ids}
        
        if self.calculate_noise_from_individual_maps:
            self.noise_from_individual_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        elif self.calculate_noise_from_coadded_maps:
            self.noise_from_coadded_maps = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
        else:
            # Even though noise is not going to be calculated,
            # the map frames that will be read may contain
            # some information on noise that I want to pass
            # to the output file.
            self.noise_data = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}

        
        
        if self.calculate_pointing_discrepancies:
            self.delta_ras_from_1st_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            self.delta_ras_from_2nd_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            self.delta_ras_from_3rd_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            self.delta_decs_from_1st_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            self.delta_decs_from_2nd_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            self.delta_decs_from_3rd_sources = \
                {map_id: core.G3MapMapDouble() for map_id in self.map_ids}
            
            self.mapping_between_dicts = \
                {"1": {"delta_ra" : self.delta_ras_from_1st_sources,
                       "delta_dec": self.delta_decs_from_1st_sources},
                 "2": {"delta_ra" : self.delta_ras_from_2nd_sources,
                       "delta_dec": self.delta_decs_from_2nd_sources},
                 "3": {"delta_ra" : self.delta_ras_from_3rd_sources,
                       "delta_dec": self.delta_decs_from_3rd_sources}}
    
    
    
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Observation:
            print()
            print("-----------------------------------------------------")
            print("Found an observation frame!")
            print("(Probably a new g3 file has arrived to the pipeline.)")
            print("-----------------------------------------------------")
            print()
            self.obs_frame = frame
        
        if frame.type == core.G3FrameType.Map:
            
            # - Load a map frame
            
            if "CoaddedMaps" in frame.keys():
                if (frame["Id"] in self.id_mapping.keys()) or \
                   (frame["Id"] in self.id_mapping.values()):
                    if frame["Id"] in self.id_mapping.values():
                        mapped_id = frame["Id"]
                    elif frame["Id"] in self.id_mapping.keys():
                        mapped_id = self.id_mapping[frame["Id"]]
                    obs_ids_from_this_frame = \
                        {mapped_id: frame["CoaddedObservationIDs"]}
                    map_ids_from_this_frame = \
                        {mapped_id: frame["CoaddedMapIDs"]}
                    if "NoiseFromIndividualMaps" in frame.keys():
                        existing_noise_data = \
                            {mapped_id: frame["NoiseFromIndividualMaps"]}
                    elif "NoiseFromCoaddedMaps" in frame.keys():
                        existing_noise_data = \
                            {mapped_id: frame["NoiseFromCoaddedMaps"]}
                    elif "NoiseData" in frame.keys():
                        existing_noise_data = \
                            {mapped_id: frame["NoiseData"]}
                    else:
                        existing_noise_data = \
                            {mapped_id: core.G3MapMapDouble()}
                    if self.calculate_noise_from_individual_maps:
                        self.noise_from_individual_maps = \
                            combine_mapmapdoubles(
                                self.noise_from_individual_maps,
                                existing_noise_data)
                    elif self.calculate_noise_from_coadded_maps:
                        self.noise_from_coadded_maps = \
                            combine_mapmapdoubles(
                                self.noise_from_coadded_maps,
                                existing_noise_data)
                    else:
                        self.noise_data = \
                            combine_mapmapdoubles(
                                self.noise_data,
                                existing_noise_data)
                    if "DeltaDecsFromSources1" in frame.keys():
                        for rank in ["1", "2", "3"]:
                            for coordinate in ["Ra", "Dec"]:
                                relevant_key = "Delta" + coordinate + \
                                               "sFromSources" + rank
                                existing_pointing_data = \
                                    {mapped_id: frame[relevant_key]}
                                relevant_dict = self.mapping_between_dicts \
                                                [rank] \
                                                ["delta_"+coordinate.lower()]
                                relevant_dict = \
                                    combine_mapmapdoubles(
                                        relevant_dict,
                                        existing_pointing_data)
                
                else:
                    print()
                    print("* Skipping the map frame above")
                    print("* because its map ID is not one of the desired IDs.")
                    print("\n")
                    return []
            
            elif frame["Id"] not in self.id_mapping.keys():
                print()
                print("* Skipping the map frame above")
                print("* because its map ID is not one of the desired IDs.")
                print("\n")
                return []
            elif self.obs_frame is None:
                print()
                print("* Skipping the map frame above")
                print("* because it is unclear what type of observation")
                print("* the map is from.")
                print("\n")
                return []
            elif self.obs_frame["SourceName"] not in self.map_sources:
                print()
                print("* Skipping the map frame above")
                print("* because the observation type is not of interest.")
                print("\n")
                return []
            else:
                mapped_id   = self.id_mapping[frame["Id"]]
                this_obs_id = self.obs_frame["ObservationID"]
                obs_ids_from_this_frame = \
                    {mapped_id: {self.obs_frame["SourceName"]: \
                                 core.G3VectorInt([this_obs_id])}}
                map_ids_from_this_frame = \
                    {mapped_id: {self.obs_frame["SourceName"]: \
                                 core.G3VectorString([str(this_obs_id)+\
                                                      frame["Id"]])}}
            
            
            # - Perform a few checks on the map
            
            if not map_seems_fine(frame, self.tqu_type,
                                  self.dict_thresh, mapped_id):
                print()
                print("* Well, the map doesn't look good,")
                print("* so, this one will be skipped...")
                print("\n")
                return []
            
            self.coadded_obs_ids, this_frame_already_added_check_one = \
                combine_info_on_ids(self.coadded_obs_ids,
                                    obs_ids_from_this_frame)
            self.coadded_map_ids, this_frame_already_added_check_two = \
                combine_info_on_ids(self.coadded_map_ids,
                                    map_ids_from_this_frame)
            
            if this_frame_already_added_check_two:
                print()
                print("* Well, it looks that this map frame has")
                print("* (at least partially) already been added,")
                print("* so, it will be skipped...")
                print("\n")
                return []
            
            
            # - Coadd maps                          
            
            if self.coadded_map_frames[mapped_id] is None:
                self.coadded_map_frames[mapped_id] = \
                    core.G3Frame(core.G3FrameType.Map)
                print()
                print("* A map frame with ID", frame["Id"], "is here!")
                print("* Both the sky map and weight map are added to the ")
                print("* empty coadded maps as", mapped_id, "...")
                print("\n")
                if self.tqu_type == "T":
                    self.coadded_map_frames[mapped_id]["T"] = frame["T"]
                    self.coadded_map_frames[mapped_id]["Wunpol"] = \
                        core.G3SkyMapWeights()
                    self.coadded_map_frames[mapped_id]["Wunpol"].TT = \
                        frame["Wunpol"].TT
            else:
                print()
                print("* Adding another set of sky map and weight map",
                      "("+frame["Id"]+")")
                print("* to the coadded maps for", mapped_id, "...")
                print("\n")
                if self.tqu_type == "T":
                    existing_t_map = \
                        self.coadded_map_frames[mapped_id].pop("T")
                    self.coadded_map_frames[mapped_id]["T"] = \
                        existing_t_map + frame["T"]
                    
                    existing_w_map = \
                        self.coadded_map_frames[mapped_id].pop("Wunpol")
                    self.coadded_map_frames[mapped_id]["Wunpol"] = \
                        existing_w_map + frame["Wunpol"]
            
            
            # - Calculate noise and/or pointing discrepancies
            
            if self.calculate_noise_from_coadded_maps    or \
               self.calculate_noise_from_individual_maps or \
               self.calculate_pointing_discrepancies:
                
                if self.obs_frame is None:
                    return []
                
                self.map_frame_arrival_counters[frame["Id"]] = 1
                if self.calculate_noise_from_individual_maps or \
                   self.calculate_pointing_discrepancies:
                    self.map_frame_cache[frame["Id"]] = frame
                
                for map_id in self.direction_independent_map_ids:
                    if "SomeDirection" in map_id:
                        id_l_ver = map_id.replace("SomeDirection", "Left")
                        id_r_ver = map_id.replace("SomeDirection", "Right")
                    else:
                        id_l_ver = "Left"  + map_id
                        id_r_ver = "Right" + map_id
                    
                    if (self.map_frame_arrival_counters[id_l_ver] == 1) and \
                       (self.map_frame_arrival_counters[id_r_ver] == 1):
                        if self.calculate_noise_from_individual_maps or \
                           self.calculate_noise_from_coadded_maps:
                            print("* Noise calculations are about to start")
                            print("* because both", id_l_ver, "and", id_r_ver)
                            print("* of obs.", self.obs_frame["ObservationID"])
                            print("* have passed the pipeline...")
                        
                        if self.calculate_noise_from_individual_maps:
                            noise = calculate_noise_from_difference_map(
                                        self.map_frame_cache[id_r_ver],
                                        self.map_frame_cache[id_l_ver],
                                        self.tqu_type)
                        if self.calculate_noise_from_coadded_maps:
                            noise = calculate_noise_from_difference_map(
                                        self.coadded_map_frames[id_r_ver],
                                        self.coadded_map_frames[id_l_ver],
                                        self.tqu_type)
                        
                        if self.calculate_noise_from_individual_maps or \
                           self.calculate_noise_from_coadded_maps:
                            print("* ...the noise level was calculated to be")
                            n = noise/(core.G3Units.uK*core.G3Units.arcmin)
                            print("*", n, "uK.arcmin.")
                            print("\n")
                        
                        
                        if self.calculate_pointing_discrepancies:
                            print()
                            print("* Pointing discrepancy calculations")
                            print("* are about to start because both")
                            print("*", id_l_ver, "and", id_r_ver)
                            print("* of obs.", self.obs_frame["ObservationID"])
                            print("* have passed the pipeline...")
                            
                            discrep_dict = calculate_pointing_discrepancies(
                                               self.map_frame_cache[id_r_ver],
                                               self.map_frame_cache[id_l_ver],
                                               self.tqu_type,
                                               self.obs_frame["SourceName"])
                            for rank, diffs in discrep_dict.items():
                                for coord, diff in diffs.items():
                                    print("* Source", rank, coord, diff, "arcsec.")
                            print("* Done.")
                            print()
                        
                        
                        if self.combine_left_right:
                            ids_for_recording = [map_id]
                        else:
                            ids_for_recording = [id_r_ver, id_l_ver]
                        
                        src = self.obs_frame["SourceName"]
                        oid = str(self.obs_frame["ObservationID"])
                        
                        for each_id in ids_for_recording:
                            if self.calculate_noise_from_individual_maps:
                                appropriate_dict = \
                                    self.noise_from_individual_maps
                            elif self.calculate_noise_from_coadded_maps:
                                appropriate_dict = \
                                    self.noise_from_coadded_maps
                            else:
                                appropriate_dict = None
                            if appropriate_dict is not None:
                                if src not in appropriate_dict[each_id].keys():
                                    appropriate_dict[each_id][src] = \
                                            core.G3MapDouble()
                                appropriate_dict[each_id][src][oid] = noise
                            
                            
                            if self.calculate_pointing_discrepancies:
                                for rank, diffs in discrep_dict.items():
                                    for coord, val in diffs.items():
                                        appropriate_dict = \
                                            self.mapping_between_dicts \
                                            [rank][coord]
                                        if src not in \
                                        appropriate_dict[each_id].keys():
                                            appropriate_dict \
                                            [each_id][src] = \
                                                    core.G3MapDouble()
                                        appropriate_dict \
                                        [each_id][src][oid] = val
                            
                            
                            if self.calculate_noise_from_individual_maps or \
                               self.calculate_pointing_discrepancies:
                                self.map_frame_cache[each_id] = None
                        
                        self.map_frame_arrival_counters[id_r_ver] = 0
                        self.map_frame_arrival_counters[id_l_ver] = 0
            
            del frame
            gc.collect()
            return [] 
        
        
        # - Compile data gathered so far
        
        if frame.type == core.G3FrameType.EndProcessing:
            print("\n")
            for map_id in self.map_ids:
                map_frame = self.coadded_map_frames[map_id]
                if map_frame is None:
                    self.coadded_map_frames.pop(map_id)
                    continue
                print()
                print("# ------------------------")
                print("#", "Map ID:", map_id )
                print("# ------------------------")
                print()
                map_frame["Id"] = map_id
                map_frame["CoaddedMaps"] = True
                map_frame["CoaddedObservationIDs"] = \
                    self.coadded_obs_ids[map_id]
                map_frame["CoaddedMapIDs"]         = \
                    self.coadded_map_ids[map_id]
                if self.calculate_noise_from_individual_maps:
                    map_frame["NoiseFromIndividualMaps"] = \
                        self.noise_from_individual_maps[map_id]
                elif self.calculate_noise_from_coadded_maps:
                    map_frame["NoiseFromCoaddedMaps"] = \
                        self.noise_from_coadded_maps[map_id]
                else:
                    map_frame["NoiseData"] = self.noise_data[map_id]
                
                if self.calculate_pointing_discrepancies:
                    map_frame["DeltaRasFromSources1"]  = \
                        self.mapping_between_dicts["1"]["delta_ra"][map_id]
                    map_frame["DeltaDecsFromSources1"] = \
                        self.mapping_between_dicts["1"]["delta_dec"][map_id]
                    map_frame["DeltaRasFromSources2"]  = \
                        self.mapping_between_dicts["2"]["delta_ra"][map_id]
                    map_frame["DeltaDecsFromSources2"] = \
                        self.mapping_between_dicts["2"]["delta_dec"][map_id]
                    map_frame["DeltaRasFromSources3"] = \
                        self.mapping_between_dicts["3"]["delta_ra"][map_id]
                    map_frame["DeltaDecsFromSources3"] = \
                        self.mapping_between_dicts["3"]["delta_dec"][map_id]
                
                
                print("These are the observations used for this map ID:")
                for source, obs_ids \
                in  map_frame["CoaddedObservationIDs"].items():
                    print("-", source)
                    print(list(obs_ids))
                for source, map_ids \
                in  map_frame["CoaddedMapIDs"].items():
                    print("-", source)
                    print(list(map_ids))
                print()
                
                if self.calculate_noise_from_individual_maps:
                    noise_dict = map_frame["NoiseFromIndividualMaps"]
                elif self.calculate_noise_from_coadded_maps:
                    noise_dict = map_frame["NoiseFromCoaddedMaps"]
                elif "NoiseData" in map_frame:
                    noise_dict = map_frame["NoiseData"]
                else:
                    noise_dict = None
                if noise_dict is not None:
                    print("These are the map noise levels (uK.arcmin):")
                    for source, sub_noise_dict \
                    in  noise_dict.items():
                        print("-", source)
                        for obs_id, noise in sub_noise_dict.items():
                            print("   ", obs_id,
                                  noise/(core.G3Units.uK*core.G3Units.arcmin))
                    print("\n")
                
                if self.calculate_pointing_discrepancies:
                    print("These are the pointing discrepancies:")
                    for sub_field in \
                    map_frame["DeltaRasFromSources1"].keys():
                        print("-", sub_field)
                        print("   - Source 1, Ra discrepancies")
                        for obs_id, delta_ra in \
                        map_frame["DeltaRasFromSources1"][sub_field].items():
                            print("     ", obs_id, delta_ra)
                        print("   - Source 1, Dec discrepancies")
                        for obs_id, delta_dec in \
                        map_frame["DeltaDecsFromSources1"][sub_field].items():
                            print("     ", obs_id, delta_dec)
                        print("   - Source 2, Ra discrepancies")
                        for obs_id, delta_ra in \
                        map_frame["DeltaRasFromSources2"][sub_field].items():
                            print("     ", obs_id, delta_ra)
                        print("   - Source 2, Dec discrepancies")
                        for obs_id, delta_dec in \
                        map_frame["DeltaDecsFromSources2"][sub_field].items():
                            print("     ", obs_id, delta_dec)
                        print("   - Source 3, Ra discrepancies")
                        for obs_id, delta_ra in \
                        map_frame["DeltaRasFromSources3"][sub_field].items():
                            print("     ", obs_id, delta_ra)
                        print("   - Source 3, Dec discrepancies")
                        for obs_id, delta_dec in \
                        map_frame["DeltaDecsFromSources3"][sub_field].items():
                            print("     ", obs_id, delta_dec)   
                    print("\n")
                
                
                print("Here is what the frame looks like:")
                print(map_frame)
                print()
            
            
            print("----------------------------------------")
            print()
            
            return list(self.coadded_map_frames.values()) + [frame]


# ==============================================================================





# ==============================================================================
# Construct a pipeline that makes a coadded map
# ------------------------------------------------------------------------------

pipeline = core.G3Pipeline()

pipeline.Add(core.G3Reader,
             filename=all_good_g3_files)

pipeline.Add(core.Dump)

pipeline.Add(CoaddMapsAndCalculateNoises,
             map_sources=arguments.sources,
             map_ids=arguments.map_ids,
             tqu_type=tqu_type,
             too_big_map_vals=too_big_map_vals,
             combine_left_right=arguments.combine_left_right,
             combine_different_wafers=arguments.combine_different_wafers,
             calculate_noise_from_individual_maps=\
                 arguments.calculate_noise_from_individual_maps,
             calculate_noise_from_coadded_maps=\
                 arguments.calculate_noise_from_coadded_maps,
             calculate_pointing_discrepancies=\
                 arguments.calculate_pointing_discrepancies)

pipeline.Add(lambda frame: "CoaddedMaps" in frame)

pipeline.Add(core.G3Writer,
             filename=arguments.output_file)


# ==============================================================================





# ==============================================================================
# Time to start!
# ------------------------------------------------------------------------------

print("\n")
print("# ======================== #")
print("# Starting the pipeline... #")
print("# ======================== #")
print()

pipeline.Run(profile=True)

print("\n\n")


# ==============================================================================
