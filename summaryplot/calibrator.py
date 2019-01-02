from spt3g import core
from statistics import compute_median, compute_nalive

def median_cal_sn_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponseSN', boloprops, selector_dict)

def median_cal_response_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponse' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_median(frame, 'CalibratorResponse', boloprops, selector_dict)

def alive_bolos_cal_4Hz(frame, boloprops, selector_dict):
    if 'CalibratorResponseSN' not in frame.keys():
        return None
    elif round(frame["CalibratorResponseFrequency"] / core.G3Units.Hz, 2) != 4.0:
        return None
    return compute_nalive(frame, 'CalibratorResponseSN', boloprops, selector_dict, 20)

