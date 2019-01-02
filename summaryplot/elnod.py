from statistics import compute_median, compute_nalive

def median_elnod_sn_slope(frame, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_median(frame, 'ElnodSNSlopes', selector_dict)

def median_elnod_iq_phase_angle(frame, selector_dict):
    if 'ElnodEigenvalueDominantVectorQ' not in frame.keys() and \
       'ElnodEigenvalueDominantVectorI' not in frame.keys():
        return None

    newframe = core.G3Frame()
    newframe['ElnodPhaseAngle'] = core.G3MapDouble()
    for bolo in frame['ElnodEigenvalueDominantVectorQ'].keys():
        if frame['ElnodEigenvalueDominantVectorI'][bolo] != 0 and \
                frame['ElnodEigenvalueDominantVectorQ'][bolo] != 0:
            newframe['ElnodPhaseAngle'][bolo] = 180/np.pi * \
                np.arctan(frame['ElnodEigenvalueDominantVectorQ'][bolo] / \
                          frame['ElnodEigenvalueDominantVectorI'][bolo])

    return compute_median(newframe, 'ElnodPhaseAngle', selector_dict)

def alive_bolos_elnod(frame, selector_dict):
    if 'ElnodSNSlopes' not in frame.keys():
        return None
    return compute_nalive(frame, 'ElnodSNSlopes', selector_dict, 20)

