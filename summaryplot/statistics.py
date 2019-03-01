import numpy as np
from functools import reduce
import operator


# functions that define quantities to be saved
def compute_median(frame, datakey, boloprops, selector_dict):
    data_list = np.array([frame[datakey][bolo] 
                          for bolo in frame[datakey].keys()])
    data_on_selection = {}

    if len(frame[datakey]) == 0:
        return np.array([])

    for select_values, f_select_list in selector_dict.items():
        selection = np.array([True for bolo in frame[datakey].keys()])
        
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, data_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], data_on_selection)[keylist[-1]] = {}

        # compute the data selection
        for select_val, f_select in zip(select_values, f_select_list):
            selection = np.array([f_select(boloprops, bolo, select_val)
                                  for bolo in frame[datakey].keys()]) & selection
        
        # get data that satisfies the selection and compute median
        reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
            np.median(data_list[selection][np.isfinite(data_list[selection])])

    return data_on_selection


def compute_nalive(frame, datakey, boloprops, selector_dict, sn_threshold, operation):
    data_list = np.array([frame[datakey][bolo] 
                          for bolo in frame[datakey].keys()])
    nalive_on_selection = {}

    if len(frame[datakey]) == 0:
        return np.array([])

    for select_values, f_select_list in selector_dict.items():
        selection = np.array([True for bolo in frame[datakey].keys()])
        
        # prep the dictionary for output
        for keylist in [select_values[:j] for j in np.arange(len(select_values))+1]:
            try:
                reduce(operator.getitem, keylist, nalive_on_selection)
            except:
                reduce(operator.getitem, keylist[:-1], nalive_on_selection)[keylist[-1]] = {}

        # compute the data selection
        for select_val, f_select in zip(select_values, f_select_list):
            selection = np.array([f_select(boloprops, bolo, select_val)
                                  for bolo in frame[datakey].keys()]) & selection
        
        # get data that satisfies the selection and compute # alive bolos
        data_on_selection = data_list[selection][np.isfinite(data_list[selection])]
        reduce(operator.getitem, select_values[:-1], nalive_on_selection)[select_values[-1]] = \
            len(data_on_selection[operation(data_on_selection, sn_threshold)])

    return nalive_on_selection
