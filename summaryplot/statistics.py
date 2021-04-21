import numpy as np
from functools import reduce
import operator
from spt3g import core


def convert_selector_dict(key_list, boloprops, selector_dict):
    selections = {}

    wafers = np.array([boloprops[bolo].wafer_id for bolo in key_list])
    bands = np.array([boloprops[bolo].band / core.G3Units.GHz for bolo in key_list])

    for wafer, band in selector_dict.keys():
        selection = np.ones_like(key_list, dtype=bool)
        selection &= bands == band
        if wafer != 'all':
            selection &= wafers == wafer

        selections[(wafer, band)] = selection

    return selections


# functions that define quantities to be saved
def compute_median(frame, datakey, boloprops, selector_dict):
    if len(frame[datakey]) == 0:
        return np.array([])

    key_list = np.array(frame[datakey].keys())
    data_list = np.array(frame[datakey].values())
    data_on_selection = {}

    selections = convert_selector_dict(key_list, boloprops, selector_dict)

    for select_values, selection in selections.items():

        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # get data that satisfies the selection and compute median
        reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
            np.median(data_list[selection][np.isfinite(data_list[selection])])

    return data_on_selection


def compute_nalive(frame, datakey, boloprops, selector_dict, sn_threshold, operation):
    if len(frame[datakey]) == 0:
        return np.array([])

    key_list = np.array(frame[datakey].keys())
    data_list = np.array(frame[datakey].values())
    nalive_on_selection = {}

    selections = convert_selector_dict(key_list, boloprops, selector_dict)

    for select_values, selection in selections.items():

        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, nalive_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], nalive_on_selection)[keylist[-1]] = {}

        # get data that satisfies the selection and compute # alive bolos
        data_on_selection = data_list[selection][np.isfinite(data_list[selection])]
        reduce(operator.getitem, select_values[:-1], nalive_on_selection)[select_values[-1]] = \
            len(data_on_selection[operation(data_on_selection, sn_threshold)])

    return nalive_on_selection
