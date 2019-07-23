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

parser.add_argument("-d", "--directory-for-file-saving",
                    action="store", type=str, default=None,
                    help="The directory in which the auxiliary files "
                         "will be saved.")

parser.add_argument("-p", "--point-source-list-file",
                    action="store", type=str, default=None,
                    help="The path to a point source list file, "
                         "which will be used for creating point source masks.")

parser.add_argument("-f", "--planck-map-fits-file",
                    action="store", type=str, default=None,
                    help="The path to the FITS files that contain Planck maps. "
                         "It actually does not point to any actual file. "
                         "Rather, the path should contain the word 'BAND' "
                         "in it so that the paths of the actual files can be "
                         "obtained by replacing 'BAND' with "
                         "'100GHz', '143GHz', or '217GHz'.")

parser.add_argument("-m", "--spt-map-data-file",
                    action="store", type=str, default=None,
                    help="The path to any g3 file that would be supplied to "
                         "fields_coadding.py. The map parameters stored there "
                         "will be necessary to create the auxiliary files.")

arguments = parser.parse_args()


# ==============================================================================




# ==============================================================================
# Generate and save the auxiliary files
# ------------------------------------------------------------------------------


# Define some global variables
# ------------------------------------------------------------------------------

deg = core.G3Units.degrees
dec_centers = [-44.75 * deg, -52.25 * deg,
               -59.75 * deg, -67.25 * deg]
ra_center = 0.0 * deg



# Generate pickle files storing which pixels should be used when calculating
# a map standard deviation, a mean weight, and so on
# ------------------------------------------------------------------------------

print("")
print("Pixel numbers to be used for each sub-field\n"
      "when calculating map fluctuation metrics are about to be found!")

spt_map_frame = list(core.G3File(arguments.spt_map_data_file))[-1]

for dec_center in dec_centers:
    
    output_file = \
        os.path.join(arguments.directory_for_file_saving,
                     "pixel_numbers_for_calculating_fluctuation_metrics_of_" + \
                     "ra0hdec{}_sub_field.pickle".format(dec_center/deg))
    
    print("- Identifying the pixels for the ra0hdec%s sub field ..."
          %(dec_center/deg))
    
    pixels_to_use = \
        fields_coadding.identify_pixels_of_non_atypical_region(
            spt_map_frame, dec_center, arguments.point_source_list_file)
    
    with open(output_file, "wb") as f_obj:
        pickle.dump(pixels_to_use, f_obj)


print("All pixels have been identified!")
print("")



# Generate Planck mini maps for assessing SPT maps' calibration
# ------------------------------------------------------------------------------

print("")
print("Mini Planck maps for asssessing SPT calibration are about to be made!")

class MakeMiniPlanckMaps(object):
    
    def __init__(self):
        
        self.spt_map_frame = None
    
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
                    ra_center, dec_center, t_only=True)
            return [mini_planck_map_frame, frame]


spt_to_planck_bands = {"90": "100", "150": "143", "220": "217"}

for spt_band, planck_band in spt_to_planck_bands.items():
    planck_map_fits_file = \
        arguments.planck_map_fits_file.replace("BAND", planck_band+"GHz")
    
    for dec_center in dec_centers:
        
        output_file = \
            os.path.join(arguments.directory_for_file_saving,
                         "mini_planck_map_for_spt_{}GHz_ra0hdec{}_sub_field.g3".\
                         format(spt_band, dec_center/deg))
        
        print("- Making a mini Planck map for %sGHz of ra0hdec%s sub field ..."
              %(spt_band, dec_center/deg))
        
        pipeline = core.G3Pipeline()
        
        pipeline.Add(core.G3Reader,
                     filename=arguments.spt_map_data_file)
        
        pipeline.Add(lambda frame: frame.type == core.G3FrameType.Map)
        
        pipeline.Add(MakeMiniPlanckMaps)
        
        pipeline.Add(core.G3Writer,
                     filename=output_file)
        
        pipeline.Run()


print("All mini Planck maps have been made!")
print("")



# Generate masks used for power spectra calculations
# ------------------------------------------------------------------------------

print("")
print("Mini masks for power spectrum calculations are about to be made!")

class MakeMiniMasksForPowerSpectraCalculations(object):
    
    def __init__(self):
        
        self.spt_map_frame = None
    
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
                    self.spt_map_frame, ra_center, dec_center, t_only=True)
            print("%4s Combinig a point source mask and "
                  "a window function ..." %(""))
            mini_mask = \
                fields_coadding.create_mask_for_powspec_calc_of_small_region(
                    mini_map_frame, arguments.point_source_list_file)
            mask_frame = core.G3Frame(core.G3FrameType.Map)
            mask_frame["Mask"] = mini_mask
            return [mask_frame, frame]



for dec_center in dec_centers:
    
    output_file = \
        os.path.join(arguments.directory_for_file_saving,
                     "mask_for_power_spectrum_calculations_for_ra0hdec{}_"
                     "sub_field.g3".format(str(dec_center/deg)))
    
    print("- Making a mask for the ra0hdec%s sub field ..."
          %(dec_center/deg))
    
    pipeline = core.G3Pipeline()
    
    pipeline.Add(core.G3Reader,
                 filename=arguments.spt_map_data_file)
    
    pipeline.Add(lambda frame: frame.type == core.G3FrameType.Map)
    
    pipeline.Add(MakeMiniMasksForPowerSpectraCalculations)
    
    pipeline.Add(core.G3Writer,
                 filename=output_file)
    
    pipeline.Run()


print("All mini masks have been made!")
print("")


# ==============================================================================
