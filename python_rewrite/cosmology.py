"""
cosmology module mimicing the fortran code
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad


from user_routines import P



def pweighted_volume_simple(z_target, omega0=0.3, Pparam=1.0, zref=0.0, h=0.7):
    """
    Simplified version matching the Fortran implementation more closely
    Uses comoving volume shells
    """
    # Create cosmology object
    cosmo = FlatLambdaCDM(H0=100*h, Om0=omega0)

    # P function normalization
    P_ref = P(zref, Pparam)

    # Create redshift array from 0 to z_target
    nz = 1000
    z_array = np.linspace(0, z_target, nz)

    if z_target == 0:
        return 0.0

    pv = 0.0
    for i in range(1, len(z_array)):
        z_mid = (z_array[i] + z_array[i-1]) / 2

        # Comoving distances in Mpc
        r1 = cosmo.comoving_distance(z_array[i-1]).value
        r2 = cosmo.comoving_distance(z_array[i]).value

        # Volume of shell
        volume_shell = (4/3) * np.pi * (r2**3 - r1**3)

        # P-weighted factor
        P_ratio = P(z_mid, Pparam) / P_ref

        pv += P_ratio * volume_shell

    return pv

def pweighted_volumes(redshifts,  omega0=0.3, Pparam=1.0, zref=0.0, h=0.7):
    """
    Above function but for multiple redshifts
    """
    return np.array([pweighted_volume_simple(z, omega0, Pparam, zref, h) for z in redshifts])
