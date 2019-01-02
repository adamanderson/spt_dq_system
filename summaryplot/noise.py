from statistics import compute_median

def median_nei_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NEI_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_0.1Hz_to_0.5Hz', boloprops, selector_dict)    

def median_nei_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NEI_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_1.0Hz_to_2.0Hz', boloprops, selector_dict)    

def median_nei_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NEI_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_3.0Hz_to_5.0Hz', boloprops, selector_dict)    

def median_nei_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NEI_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEI_10.0Hz_to_15.0Hz', boloprops, selector_dict)    

def median_nep_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NEP_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_0.1Hz_to_0.5Hz', boloprops, selector_dict)    

def median_nep_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NEP_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_1.0Hz_to_2.0Hz', boloprops, selector_dict)    

def median_nep_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NEP_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_3.0Hz_to_5.0Hz', boloprops, selector_dict)    

def median_nep_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NEP_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NEP_10.0Hz_to_15.0Hz', boloprops, selector_dict)    

def median_net_01Hz_to_05Hz(frame, boloprops, selector_dict):
    if 'NET_0.1Hz_to_0.5Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_0.1Hz_to_0.5Hz', boloprops, selector_dict)    

def median_net_1Hz_to_2Hz(frame, boloprops, selector_dict):
    if 'NET_1.0Hz_to_2.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_1.0Hz_to_2.0Hz', boloprops, selector_dict)    

def median_net_3Hz_to_5Hz(frame, boloprops, selector_dict):
    if 'NET_3.0Hz_to_5.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_3.0Hz_to_5.0Hz', boloprops, selector_dict)    

def median_net_10Hz_to_15Hz(frame, boloprops, selector_dict):
    if 'NET_10.0Hz_to_15.0Hz' not in frame.keys():
        return None
    return compute_median(frame, 'NET_10.0Hz_to_15.0Hz', boloprops, selector_dict)    
