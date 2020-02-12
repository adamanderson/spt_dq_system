import numpy as np
from spt3g import core

def get_elevation(fname):
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=fname)

    elevation = None
    class CalcEl(object):
        def __init__(self):
            self.elevation = None
        def __call__(self, frame):
            if frame.type == core.G3FrameType.Scan and \
               'OnlineBoresightEl' in frame.keys():
                self.elevation = np.median([frame["OnlineBoresightEl"][i]*180/np.pi \
                                            for i in range(len(frame["OnlineBoresightEl"]))])
            print(self.elevation)
    calc_el = CalcEl()
    pipe.Add(calc_el)
    pipe.Run()

    return calc_el.elevation
