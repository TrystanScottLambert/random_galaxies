"""
Sub routines for the random catlaogues
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM


from user_routines import kcorr, ecorr, decorr_du
from cosmology import pweighted_volumes
from data_classes import SurveySpec, Galaxy

def kpluse(z, u, zref):
    """
    Combined k+e correction normalized to zero at zref.

    Parameters:
        z (float): Redshift
        u (float): Luminosity evolution parameter
        zref (float): Reference redshift

    Returns:
        float: The combined k+e correction
    """
    return kcorr(z) + ecorr(z, u) - (kcorr(zref) + ecorr(zref, u))

_first_call = True  # Global flag for dkpluse_du

def dkpluse_du(z, u):
    """
    Derivative of the combined k+e correction with respect to u.

    Parameters:
        z (float): Redshift
        u (float): Luminosity evolution parameter

    Returns:
        float: Derivative of k+e correction with respect to u
    """
    global _first_call
    eps = 1.0e-4

    if _first_call:
        import numpy as np
        z_vals = np.arange(0.0, 2.01, 0.1)
        u_vals = np.arange(-2.0, 2.01, 0.1)

        for zt in z_vals:
            for ut in u_vals:
                numerical = (ecorr(zt, ut + eps) - ecorr(zt, ut)) / eps
                analytic = decorr_du(zt, ut)
                err = abs(numerical - analytic)
                if err > 0.01:
                    print(f"Warning: supplied decorr_du not consistent with numerical derivative at z={zt:.2f}, u={ut:.2f} (err={err:.3f})")

        _first_call = False

    return decorr_du(z, u)


def tabulate_bin_pvolumes(redshift_bins, fractional_area, zref, Pparam, omega0=0.3):
    """
    Computes P-weighted and standard comoving volume for redshift bins.
    """
    cosmo = FlatLambdaCDM(H0=70, Om0=omega0)

    comoving_volumes = cosmo.comoving_volume(redshift_bins).value
    volbins = (comoving_volumes[1:] - comoving_volumes[:-1]) * fractional_area

    p_volumes = pweighted_volumes(redshift_bins, zref, Pparam, omega0)
    p_vol_bins = (p_volumes[1:] - p_volumes[:-1]) * fractional_area

    return volbins, p_vol_bins

def derived_gal_props(u, zref, survey: SurveySpec, catalog: list[Galaxy])