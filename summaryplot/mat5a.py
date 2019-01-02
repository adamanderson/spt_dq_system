from statistics import compute_median

def median_mat5a_fluxcal(frame, boloprops, selector_dict):
    if 'MAT5AFluxCalibration' not in frame.keys():
        return None
    return compute_median(frame, 'MAT5AFluxCalibration', boloprops, selector_dict)

def median_mat5a_intflux(frame, boloprops, selector_dict):
    if 'MAT5AIntegralFlux' not in frame.keys():
        return None
    return compute_median(frame, 'MAT5AIntegralFlux', boloprops, selector_dict)

def mat5a_sky_transmission(frame, boloprops, selector_dict):
    if 'MAT5ASkyTransmission' not in frame.keys():
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
        if str(select_values[-1]) in frame['MAT5ASkyTransmission'].keys():
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = \
                frame['MAT5ASkyTransmission'][str(select_values[-1])]
        else:
            reduce(operator.getitem, select_values[:-1], data_on_selection)[select_values[-1]] = np.nan

    return data_on_selection

