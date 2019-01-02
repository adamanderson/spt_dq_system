from statistics import compute_median
from functools import reduce
import operator

def median_rcw38_fluxcal(frame, boloprops, selector_dict):
    if 'RCW38FluxCalibration' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38FluxCalibration', boloprops, selector_dict)

def median_rcw38_intflux(frame, boloprops, selector_dict):
    if 'RCW38IntegralFlux' not in frame.keys():
        return None
    return compute_median(frame, 'RCW38IntegralFlux', boloprops, selector_dict)

def rcw38_sky_transmission(frame, boloprops, selector_dict):
    if 'RCW38SkyTransmission' not in frame.keys():
        return None

    data_on_selection = {}

    for select_values, f_select_list in selector_dict.items():
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # get data that satisfies the selection and compute median
        if str(select_values[-1]) in frame['RCW38SkyTransmission'].keys():
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
                frame['RCW38SkyTransmission'][str(select_values[-1])]
        else:
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = np.nan

    return data_on_selection
