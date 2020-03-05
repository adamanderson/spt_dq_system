# ============================================================================ #
#  This script is intended to be called by                                     #
#  spt_dq_system/summaryplot/cache_archive_data.py, which provides a list of   #
#  input ARC files to be processed.                                            #
#                                                                              #
#  Each ARC file contains about 16 minutes' worth of data collected by GCP,    #
#  and the data from most of the measured quantities have a sampling rate of   #
#  100 Hz. Since it seems unnecessary to plot all these data points,           #
#  some reduction is done by this script.                                      #
#                                                                              #
#  For each ARC file, the script extracts the data from                        #
#  the quanties to be plotted, calculates the average value for each quantity  #
#  during that period, and then saves all the average values in a pickle file. #
#  The size of each file is about 850 bytes.                                   #
#                                                                              #
# ============================================================================ #


import os
import numpy
import pickle

from spt3g import core
from spt3g import gcp
from spt3g import util
from spt3g import std_processing



def run(input_files=None,
        pickles_dir=None,
        logger=None):

    logger.info('')
    logger.info('----------------------------')
    logger.info(' These are the input files: ')
    logger.info('----------------------------')
    logger.info('')
    for counter, f in enumerate(input_files, 1):
        logger.info(' file %06d: %s', counter, f)
    logger.info('\n')


    key_dict = {'tipper_tau' : ['array', 'tipper', 'tau' ],
                'tipper_tatm': ['array', 'tipper', 'tatm'],
                'tipper_utc' : ['array', 'tipper', 'utc' ],
                'weather_tatm': ['array', 'weather', 'airTemperature'     ],
                'weather_tint': ['array', 'weather', 'internalTemperature'],
                'weather_humi': ['array', 'weather', 'relativeHumidity'   ],
                'weather_pres': ['array', 'weather', 'pressure'           ],
                'weather_wdir': ['array', 'weather', 'windDirection'      ],
                'weather_wmag': ['array', 'weather', 'windSpeed'          ],
                'weather_utc' : ['array', 'weather', 'utc'                ],
                'cb_rx_4k_head' : ['array', 'cryo', 'temperature', 0, 13],
                'cb_rx_4k_sqstr': ['array', 'cryo', 'temperature', 0, 14],
                'cb_rx_4k_plafa': ['array', 'cryo', 'temperature', 2,  0],
                'cb_rx_50k_head': ['array', 'cryo', 'temperature', 0, 15],
                'cb_rx_50k_harn': ['array', 'cryo', 'temperature', 2,  6],
                'cb_ox_4k_head' : ['array', 'cryo', 'temperature', 1,  8],
                'cb_ox_4k_lens' : ['array', 'cryo', 'temperature', 1, 10],
                'cb_ox_4k_lyot' : ['array', 'cryo', 'temperature', 1, 15],
                'cb_ox_50k_head': ['array', 'cryo', 'temperature', 1,  4],
                'cb_ox_wbp_far' : ['array', 'cryo', 'temperature', 1,  1],
                'cb_ox_50k_tube': ['array', 'cryo', 'temperature', 1,  7],
                'cb_utc'        : ['array', 'cryo', 'utc'],
                'scu_tcryoshell': ['antenna0', 'scu', 'temp', 20],
                'scu_ticecrate' : ['antenna0', 'scu', 'temp', 21],
                'scu_tsecondary': ['antenna0', 'scu', 'temp', 22],
                'scu_tbackerack': ['antenna0', 'scu', 'temp', 23],
                'scu_tcabin'    : ['antenna0', 'scu', 'temp', 24],
                'scu_tundercb'  : ['antenna0', 'scu', 'temp', 25],
                'scu_utc'       : ['antenna0', 'scu', 'tempSampleTime']}


    for counter_f, f in enumerate(input_files, 1):
        logger.info('Processing %s (file %06d)...', f, counter_f)
        arc_data_collector = util.extractdata.MultiAccumulator(
                                 keys=key_dict)
        pipeline = core.G3Pipeline()
        pipeline.Add(gcp.ARCFileReader,
                     filename=f)
        pipeline.Add(arc_data_collector)
        
        try:
            pipeline.Run()
        except:
            logger.info('')
            logger.exception('An error occurred:\n------------------')
            logger.info('------------------')
            logger.info('This file was skipped!')
            logger.info('')
            continue
        
        all_data = arc_data_collector.extract_values()
        
        for key, array in all_data.items():
            if '_utc' in key:
                all_data.pop(key)
                new_key   = key.replace('_utc', '_oid')
                new_array = [std_processing.time_to_obsid(t) for t in array]
                all_data[new_key] = new_array

        indices_from_okay_times = {'tipper' : None,
                                   'weather': None,
                                   'cb'     : None,
                                   'scu'    : None}
        for k, l in all_data.items():
            if '_oid' in k:
                dtype = k.split('_')[0]
                oids  = numpy.asarray(l)
                filt  = numpy.where(numpy.isfinite(oids)&(oids>0))
                ### Some timestamps seem to be corrupted.
                ### I saw one from the 10th century...
                indices_from_okay_times[dtype] = filt
        for k, ts in all_data.items():
            good_indices = indices_from_okay_times[k.split('_')[0]]
            all_data[k]  = numpy.asarray(ts)[good_indices]

        unique_oids, indices = numpy.unique(
                                   all_data['tipper_oid'], return_index=True)
        for key, value in all_data.items():
            if 'tipper_' in key:
                all_data[key] = numpy.array(value)[indices]

        for key, value in all_data.items():
            if len(value) == 0:
                average = numpy.nan
            else:
                average = numpy.mean(numpy.asarray(value))
            all_data[key] = average

        t_obsid = std_processing.time_to_obsid(f.split('.')[0].split('/')[-1])
        output_pickle = '{:09d}.pickle'.format(t_obsid)
        output_path = os.path.join(
                          pickles_dir, output_pickle)
        with open(output_path, 'wb') as fobj:
            pickle.dump(all_data, fobj)
        logger.info('Saved the output to %s.', output_path)

