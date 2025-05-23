"""
Sub routines for the random catlaogues
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
from scipy.interpolate import interp1d


from user_routines import kcorr, ecorr, decorr_du
from cosmology import pweighted_volumes, pweighted_volume_simple
from data_classes import SurveySpec

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

def z_at_mag(absolute_magnitude: float, u, survey):
    """
    interpolates a solution becuaes of the kpluse correction
    """
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    redshifts = np.arange(survey.zmin, survey.zmax, 0.001)
    dl = cosmo.luminosity_distance(redshifts)
    apparent_magnitudes = absolute_magnitude + 5 * np.log10(dl) + 25 + kpluse(redshifts, u, zref=0.1)

    # Sort arrays to make sure interpolation is monotonic
    sort_idx = np.argsort(apparent_magnitudes)
    mags_sorted = apparent_magnitudes[sort_idx]
    z_sorted = redshifts[sort_idx]

    # Create interpolation function: m -> z
    mag_to_z = interp1d(mags_sorted, z_sorted, bounds_error=False, fill_value=np.nan)
    return mag_to_z(survey.magfaint)


def update_catalog(catalog: pd.DataFrame, zref: float, survey: SurveySpec, u: float, a: float) -> pd.DataFrame:
    """
    Updates the catalogs magnitudes and pv calculations for every galaxy.
    """
    cosmo = FlatLambdaCDM(H0=70, Om0 = 0.3)
    p_vol_min_total = pweighted_volume_simple(survey.zmin, omega0=0.3, Pparam=a, zref = 0.1, h=0.7)
    p_vol_min = p_vol_min_total*survey.survey_fractional_area

    luminosity_distances = cosmo.luminosity_distance(catalog['z']).value

    catalog['pv'] = survey.survey_fractional_area * pweighted_volumes(catalog['z'], omega0=0.3, Pparam=a, zref=0.1, h=0.7) - p_vol_min
    catalog['dm'] = 25 + 5 * np.log10(cosmo.luminosity_distance(catalog['z']).value) + kpluse(catalog['z'], u, zref) #TODO: Check that this doesn't have to be parsecs
    mag = catalog['absmag'] + 5 * np.log10(luminosity_distances) + 25 + kpluse(catalog['z'], u,zref)
    catalog['absmag'] = catalog['absmag'] - mag + catalog['mag']
    catalog['z_max'] = z_at_mag(catalog['absmag'], u, survey)
    catalog.loc[catalog['z_max'] > survey.zmax, 'z_max'] = survey.zmax - 1e-5
    catalog['pvmax'] = survey.survey_fractional_area * pweighted_volumes(catalog['z_max'], 0.3, a, 0.1, 0.7)

    return catalog