"""
Sub routines for the random catlaogues
"""

from dataclasses import dataclass

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
from scipy.interpolate import interp1d


from user_routines import kcorr, ecorr, decorr_du
from cosmology import pweighted_volumes, pweighted_volume_simple
from data_classes import SurveySpec


def initilize_bin_dataframe(
        redshift_bins: np.ndarray, cosmology: FlatLambdaCDM, survey: SurveySpec) -> pd.DataFrame:
    """
    Creates the redshift bins that we will continously updated throughout the program. 
    redshift bins is the same bin object we would pass to make a histogram. I.e, n+1 for n bins.
    """
    upper_bins = redshift_bins[1:]
    lower_bins = redshift_bins[:-1]
    mid_points = 0.5 * (upper_bins + lower_bins)

    # Get the total volume of each bin
    comoving_volumes = cosmology.comoving_volume(redshift_bins).value
    volbins = (comoving_volumes[1:] - comoving_volumes[:-1]) * survey.survey_fractional_area

    randoms = np.ones(len(mid_points))
    n = np.ones(len(mid_points))
    deltas = np.ones(len(mid_points))
    initial_data = [lower_bins, mid_points, upper_bins, volbins, randoms, n, deltas]
    df = pd.DataFrame(np.array(initial_data).T,
                      columns=['lower', 'mid', 'upper', 'volume', 'randoms', 'n', 'delta'])
    return df

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


def dkpluse_du(z, u):
    """
    Derivative of the combined k+e correction with respect to u.

    Parameters:
        z (float): Redshift
        u (float): Luminosity evolution parameter

    Returns:
        float: Derivative of k+e correction with respect to u
    """
    return decorr_du(z, u)


def check_decorr_du() -> bool:
    """
    Checks that the decorr_du is consistent with the numerical derivative.
    """
    result = True

    eps = 1e-4
    z_vals = np.arange(0.0, 2.01, 0.1)
    u_vals = np.arange(-2.0, 2.01, 0.1)

    for zt in z_vals:
        for ut in u_vals:
            numerical = (ecorr(zt, ut + eps) - ecorr(zt, ut)) / eps
            analytic = decorr_du(zt, ut)
            err = abs(numerical - analytic)
            if err > 0.01:
               result = False
    return result


def tabulate_bin_pvolumes(
        bin_edges: np.ndarray, fractional_area: float, zref: float, Pparam: float,
        cosmo: FlatLambdaCDM):
    """
    Computes P-weighted and standard comoving volume for redshift bins.
    """
    p_volumes = pweighted_volumes(bin_edges, cosmo, Pparam, zref)
    p_vol_bins = (p_volumes[1:] - p_volumes[:-1]) * fractional_area
    return p_vol_bins

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


def update_catalog(
        catalog: pd.DataFrame, zref: float, survey: SurveySpec, u: float, a: float) -> pd.DataFrame:
    """
    Updates the catalogs magnitudes and pv calculations for every galaxy.
    """
    cosmo = FlatLambdaCDM(H0=70, Om0 = 0.3)
    p_vol_min_total = pweighted_volume_simple(survey.zmin, cosmo, a, zref)
    p_vol_min = p_vol_min_total*survey.survey_fractional_area

    luminosity_distances = cosmo.luminosity_distance(catalog['z']).value

    catalog['pv'] = survey.survey_fractional_area * pweighted_volumes(catalog['z'], cosmo, a, zref) - p_vol_min
    catalog['dm'] = 25 + 5 * np.log10(cosmo.luminosity_distance(catalog['z']).value) + kpluse(catalog['z'], u, zref) #TODO: Check that this doesn't have to be parsecs
    mag = catalog['absmag'] + 5 * np.log10(luminosity_distances) + 25 + kpluse(catalog['z'], u,zref)
    catalog['absmag'] = catalog['absmag'] - mag + catalog['mag']
    catalog['z_max'] = z_at_mag(catalog['absmag'], u, survey)
    catalog.loc[catalog['z_max'] > survey.zmax, 'z_max'] = survey.zmax - 1e-5
    catalog['pvmax'] = survey.survey_fractional_area * pweighted_volumes(catalog['z_max'], cosmo, a, zref)

    return catalog


def estimate_lf(catalog: pd.DataFrame, zbin: pd.DataFrame, delta: np.ndarray[float], u: float, nlf: int):
    
    catalog['dabsmag_du'] = -dkpluse_du(catalog['z'], u)

    pvsumeff = 0.0
    pvsum = 0.0
    
    # For every catalog row, we need to work put the svsumeff and pvsum
    for i in range(len(catalog)):
        j = 0
        while catalog[i].z_max < zbin[j+1]:
            print("poo")


