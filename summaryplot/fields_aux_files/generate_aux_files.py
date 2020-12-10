# ============================================================================ #
#  This script is intended to be used to generate useful auxiliary data        #
#  such as numpy arrays specifying pixel numbers to be used for calculating    #
#  map fluctuation metrics, masks used to calculate power spectra, and         #
#  Planck mini maps used to calculate cross spectra with SPT maps.             #
#                                                                              #
#  If these data are available on disk, then fields_coadding.py does not need  #
#  to create these files, which can save time.                                 #
#                                                                              #
# ============================================================================ #


import argparse
import pickle
import os
import shutil

from spt3g       import core
from summaryplot import fields_coadding



# ==============================================================================
# Interpret command line arguments
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
             description="This script generates useful auxiliary files and "
                         "save them to disk so that fields_coadding.py can "
                         "simpy load them instead of generating them "
                         "every time the script is run.",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter,
             epilog="This should save time! (Wei)")

parser.add_argument("-s", "--season",
                    action="store", type=str, default=None,
                    help="The observing season for which "
                         "auxiliary files are to be made.")

parser.add_argument("-d", "--directory-for-file-saving",
                    action="store", type=str, default=None,
                    help="The directory in which the auxiliary files "
                         "will be saved.")

parser.add_argument("-i", "--identify-pixels-of-uniform-coverage-regions",
                    action="store_true", default=False,
                    help="Whether to record which pixels correspond to the "
                         "central, more-or-less-uniform-coverage region of "
                         "each sub-field.")

parser.add_argument("-l", "--generate-mini-planck-maps",
                    action="store_true", default=False,
                    help="Whether to make mini Planck maps used to "
                         "cross-correlate SPT maps with.")

parser.add_argument("-m", "--generate-masks-for-power-spectrum-calculations",
                    action="store_true", default=False,
                    help="Whether to create masks (point source masking and "
                         "border windowing) used for calculating various "
                         "power spectra.")

parser.add_argument("-y", "--make-dummy-g3-file",
                    action="store_true", default=False,
                    help="Whether to make a dummy g3 file that contains "
                         "no information. This can be useful when I want to "
                         "construct a pipeline that does not really require "
                         "any input files but still saves a useful g3 file.")

parser.add_argument("-p", "--point-source-list-file",
                    action="store", type=str, default=None,
                    help="The path to a point source list file, "
                         "which will be used for creating point source masks.")

parser.add_argument("-f", "--planck-map-fits-file",
                    action="store", type=str, default=None,
                    help="The path to the FITS files that contain Planck maps. "
                         "It actually does not point to any actual file. "
                         "Rather, the path should contain the words 'BAND' "
                         "and 'MISSION' in it so that the paths of the "
                         "actual files can be obtained by replacing 'BAND' "
                         "with '100GHz', '143GHz', or '217GHz' and 'MISSION' "
                         "with 'fullmission', 'halfmission-1', and so on.")

parser.add_argument("-t", "--spt-map-data-file",
                    action="store", type=str, default=None,
                    help="The path to any g3 file that would be supplied to "
                         "fields_coadding.py. The map parameters stored there "
                         "will be necessary to create the auxiliary files.")

parser.add_argument("-D", "--dummy-file-for-dummy-file",
                    action="store", type=str, default=None,
                    help="The path to really any g3 file that will be supplied "
                         "to a pipeline that creates the aforementioned dummy "
                         "g3 file.")

arguments = parser.parse_args()


# ==============================================================================




# ==============================================================================
# Generate and save the auxiliary files
# ------------------------------------------------------------------------------


# Define some global variables
# ------------------------------------------------------------------------------

deg = core.G3Units.degrees

if arguments.season == "winter":
    dec_centers = [-44.75 * deg, -52.25 * deg,
                   -59.75 * deg, -67.25 * deg]
    ra_center   =    0.0 * deg
    field_names = ["ra0hdec-44.75", "ra0hdec-52.25",
                   "ra0hdec-59.75", "ra0hdec-67.25"]
    minimap_xylen_deg = 5.0
elif arguments.season == "summer":
    dec_centers = [-24.50 * deg, -31.50 * deg,
                   -38.50 * deg, -45.50 * deg,
                   -52.50 * deg, -59.50 * deg]
    ra_center   =   75.0 * deg
    field_names = ["ra5hdec-24.5", "ra5hdec-31.5",
                   "ra5hdec-38.5", "ra5hdec-45.5",
                   "ra5hdec-52.5", "ra5hdec-59.5"]
    minimap_xylen_deg = 5.0
elif arguments.season == "summerb":
    dec_centers = [-29.75 * deg, -33.25 * deg,
                   -36.75 * deg, -40.25 * deg]
    ra_center   =   25.0 * deg
    field_names = ["ra1h40dec-29.75", "ra1h40dec-33.25",
                   "ra1h40dec-36.75", "ra1h40dec-40.25"]
    minimap_xylen_deg = 2.0



# Generate a text file used for recording bad maps found
# when making year-to-date coadded maps.
# ------------------------------------------------------------------------------

if os.path.isfile("summaryplot/fields_aux_files/bad_map_list.txt"):
    pass
else:
    print("")
    print("Preparing a text file to potentially record bad maps later!")
    shutil.copy("summaryplot/field_aux_files/bad_map_list_template.txt",
                "summaryplot/field_aux_files/bad_map_list.txt")



# Generate a dummy g3 file used for some silly pipelines invoked by
# cache_field_maps.py
# ------------------------------------------------------------------------------

if os.path.isfile("summaryplot/fields_aux_files/dummy.g3"):
    pass
else:
    print("")
    print("Creating a dummy g3 file for some silly pipelines used later!")
    pipeline = core.G3Pipeline()
    pipeline.Add(core.G3InfiniteSource,
                 n=1)
    pipeline.Add(core.G3Writer,
                 filename="summaryplot/fields_aux_files/dummy.g3")
    pipeline.Run()



# Generate pickle files storing which pixels should be used when calculating
# a map standard deviation, a mean weight, and so on
# ------------------------------------------------------------------------------

if arguments.identify_pixels_of_uniform_coverage_regions:
    
    print("")
    print("Pixel numbers to be used for each sub-field\n"
          "when calculating map fluctuation metrics are about to be found!")
    
    
    spt_map_frame = list(core.G3File(arguments.spt_map_data_file))[-1]
    
    for field_name in field_names:
        
        output_file = \
            os.path.join(arguments.directory_for_file_saving,
                         ("pixel_numbers_for_calculating_"+\
                          "fluctuation_metrics_of_" + \
                          "{}_sub_field.pickle").format(field_name))
        
        print("- Identifying the pixels for the %s sub field ..."
              %(field_name))
        
        pixels_to_use = \
            fields_coadding.identify_pixels_of_non_atypical_region(
                field_name, spt_map_frame, arguments.point_source_list_file)
        
        with open(output_file, "wb") as f_obj:
            pickle.dump(pixels_to_use, f_obj)
    
    
    print("All pixels have been identified!")
    print("")



# Generate Planck mini maps for assessing SPT maps' calibration
# ------------------------------------------------------------------------------

if arguments.generate_mini_planck_maps:
    
    print("")
    print("Mini Planck maps for asssessing SPT calibration "
          "are about to be made!")
    
    
    class MakeMiniPlanckMaps(object):
        
        def __init__(self, ra_center, dec_center, xylen_deg):
            
            self.spt_map_frame = None
            self.ra_center     = ra_center
            self.dec_center    = dec_center
            self.xylen_deg     = xylen_deg
        
        def __call__(self, frame):
            
            if frame.type == core.G3FrameType.Map:
                if "GHz" in frame["Id"]:
                    self.spt_map_frame = frame
                    print("%4s Found the map frame "
                          "that contains needed information." %(""))
                return []
            
            if frame.type == core.G3FrameType.EndProcessing:
                print("%4s Making a mini Planck map by using %s ..."
                      %("", planck_map_fits_file))
                mini_planck_map_frame = \
                    fields_coadding.create_mini_planck_map_frame(
                        planck_map_fits_file,
                        self.spt_map_frame,
                        self.ra_center, self.dec_center,
                        t_only=True, xylen_deg=self.xylen_deg)
                return [mini_planck_map_frame, frame]
    
    
    spt_to_planck_bands = {"90": "100", "150": "143", "220": "217"}
    
    for spt_band, planck_band in spt_to_planck_bands.items():
        for planck_mission in ["fullmission", "halfmission-1", "halfmission-2"]:
            
            planck_map_fits_file = \
                arguments.planck_map_fits_file.\
                replace("BAND", planck_band+"GHz").\
                replace("MISSION", planck_mission)
            
            for counter, field_name in enumerate(field_names):
                
                output_file = \
                    os.path.join(
                        arguments.directory_for_file_saving,
                        ("mini_{}_planck_map_for_"+\
                         "spt_{}GHz_{}_sub_field.g3").\
                        format(planck_mission, spt_band, field_name))
                
                print("- Making a mini %s Planck map for "
                      "%sGHz of %s sub field ..."
                      %(planck_mission, spt_band, field_name))
                
                pipeline = core.G3Pipeline()
                
                pipeline.Add(core.G3Reader,
                             filename=arguments.spt_map_data_file)
                
                pipeline.Add(lambda frame: frame.type == core.G3FrameType.Map)
                
                pipeline.Add(MakeMiniPlanckMaps,
                             ra_center=ra_center,
                             dec_center=dec_centers[counter],
                             xylen_deg=minimap_xylen_deg)
                
                pipeline.Add(core.G3Writer,
                             filename=output_file)
                
                pipeline.Run()
    
    
    print("All mini Planck maps have been made!")
    print("")



# Generate masks used for power spectra calculations
# ------------------------------------------------------------------------------

if arguments.generate_masks_for_power_spectrum_calculations:
    
    print("")
    print("Mini masks for power spectrum calculations are about to be made!")
    
    
    class MakeMiniMasksForPowerSpectraCalculations(object):
        
        def __init__(self, ra_center, dec_center, xylen_deg):
            
            self.spt_map_frame = None
            self.ra_center     = ra_center
            self.dec_center    = dec_center
            self.xylen_deg     = xylen_deg
        
        def __call__(self, frame):
            
            if frame.type == core.G3FrameType.Map:
                if "GHz" in frame["Id"]:
                    self.spt_map_frame = frame
                    print("%4s Found the map frame "
                          "that contains needed information." %(""))
                return []
            
            if frame.type == core.G3FrameType.EndProcessing:
                print("%4s Making a map frame that contains a small region "
                      "of the original map ..." %(""))
                mini_map_frame = \
                    fields_coadding.create_new_map_frame_with_smaller_region(
                        self.spt_map_frame,
                        self.ra_center, self.dec_center,
                        t_only=True, xylen_deg=self.xylen_deg)
                print("%4s Combining a point source mask and "
                      "a window function ..." %(""))
                mini_mask = \
                    fields_coadding.\
                    create_mask_for_powspec_calc_of_small_region(
                        mini_map_frame, arguments.point_source_list_file)
                mask_frame = core.G3Frame(core.G3FrameType.Map)
                mask_frame["Mask"] = mini_mask
                return [mask_frame, frame]
    
    
    
    for counter, field_name in enumerate(field_names):
        
        output_file = \
            os.path.join(arguments.directory_for_file_saving,
                         "mask_for_power_spectrum_calculations_for_{}_"
                         "sub_field.g3".format(field_name))
        
        print("- Making a mask for the %s sub field ..."
              %(field_name))
        
        pipeline = core.G3Pipeline()
        
        pipeline.Add(core.G3Reader,
                     filename=arguments.spt_map_data_file)
        
        pipeline.Add(lambda frame: frame.type == core.G3FrameType.Map)
        
        pipeline.Add(MakeMiniMasksForPowerSpectraCalculations,
                     ra_center=ra_center,
                     dec_center=dec_centers[counter],
                     xylen_deg=minimap_xylen_deg)
        
        pipeline.Add(core.G3Writer,
                     filename=output_file)
        
        pipeline.Run()
    
    
    print("All mini masks have been made!")
    print("")



# Create a dummy g3 file that contains no information
# ------------------------------------------------------------------------------

if arguments.make_dummy_g3_file:
    
    print("")
    print("A dummy g3 file is about to be made!")
    
    def insert_an_ampty_frame(frame):
        if frame.type != core.G3FrameType.EndProcessing:
            return []
        else:
            new_frame = core.G3Frame()
            return [new_frame, frame]
    
    pipeline = core.G3Pipeline()
    
    pipeline.Add(core.G3Reader,
                 filename=arguments.dummy_file_for_dummy_file)
    
    pipeline.Add(insert_an_ampty_frame)
    
    pipeline.Add(core.G3Writer,
                 filename=os.path.join(arguments.directory_for_file_saving,
                                       "dummy.g3"))
    
    pipeline.Run()
    
    print("The file was generated!.")
    print("")


# ==============================================================================
