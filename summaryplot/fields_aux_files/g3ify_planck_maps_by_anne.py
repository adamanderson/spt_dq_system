import os
import numpy as np
from spt3g.util import healpix_tools as hpt
from spt3g import core

high_pass = None
high_pass_width = 150

def mask_1500d(nside, coord=None):
    return hpt.latlon_mask(
        nside,
        # latrange=(-74, -38),
        # lonrange=(-68, 68),
        # latrange=[-66, -18],
        # lonrange=[45, 105],
        # latrange=(-45, -25),
        # lonrange=(-5, 55),
        latrange=(-45, -25),
        lonrange=(145, 230),
        coord=coord
    )

def hpf(cutoff, width, nside):
    # high pass in l-space
    lmax = 3 * nside - 1
    ell = np.arange(lmax+1)
    trans = np.logical_and(ell > cutoff-width/2, ell < cutoff+width/2)
    kernel = np.zeros_like(ell, dtype=float)
    kernel[ell <= cutoff-width/2] = 0
    kernel[ell >= cutoff+width/2] = 1
    kernel[trans] = 0.5 - 0.5 * np.cos(ell[trans]*np.pi/width - 
                                       np.pi*(cutoff-width/2)/width)
    return ell, kernel
    
mask = mask_1500d(2048)
map_path = '/sptlocal/user/agambrel/planck_maps'
mname = 'HFI_SkyMap_{}_2048_R3.01_{}.fits'

for freq in [100, 143, 217]:
    for mission in ['fullmission', 'halfmission-1', 'halfmission-2']:
        if high_pass is not None:
            mname_out = 'HFI_SkyMap_{}_2048_R3.01_{}_cut_C_G3Units_hpl{}.fits'.format(freq, hm, high_pass)
        else:
            # mname_out = 'HFI_SkyMap_{}_2048_R3.01_{}_cut_C_G3Units.fits'.format(freq, mission)
            # mname_out = 'HFI_SkyMap_{}_2048_R3.01_{}_cut_summer_C_G3Units.fits'.format(freq, mission)
            # mname_out = 'HFI_SkyMap_{}GHz_2048_R3.01_{}_cut_summerb_C_G3Units.fits'.format(freq, mission)
            mname_out = 'HFI_SkyMap_{}GHz_2048_R3.01_{}_cut_summerc_C_G3Units.fits'.format(freq, mission)

        if os.path.exists(os.path.join(map_path, mname_out)):
            continue

        print('{} GHz {}'.format(freq, mission))
        if mission == 'fullmission':
            mission = 'full'
        m = hpt.read_map(os.path.join(map_path, mname.format(freq, mission)),
                         field=None)
        if high_pass is not None:
            print('Applying high pass filter')
            ell, bl_hp = hpf(high_pass, high_pass_width, 2048)
            m[:3] = hpt.smoothing(m[:3], beam=np.tile(bl_hp, 3))
        print('Rotating from G to C')
        mr = hpt.rotate_map(m[:3], ['G', 'C'])
        mr[:,mask] *= core.G3Units.K
        print('Writing output map')
        hpt.write_map(
            os.path.join(".", mname_out.format(freq, mission)),
            mr,
            coord='C',
            mask=mask,
            column_names=['T', 'Q', 'U'],
            extra_header={'WEIGHTED': False,
                          'UNITS': 'Tcmb'},
        )
