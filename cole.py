"""
Script that will generate galaxies using the pure Cole method.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.interpolate import interp1d
import pylab as plt
import pandas as pd

from geom_calcs import calculate_area_of_rectangular_region

def ab_mag_to_luminosity(absolute_mags: np.ndarray[float]) -> np.ndarray[float]:
    """
    Converts absolute magnitudes to solar luminosties
    (https://resources.wolframcloud.com/FormulaRepository/resources/Luminosity-Formula-for-Absolute-Magnitude)
    """ 
    return np.log10(10**(0.4*(4.83 - absolute_mags)))

def convert_ap_to_ab(apparent_mags: np.ndarray[float], redshifts: np.ndarray[float], cosmology: FlatLambdaCDM) -> np.ndarray[float]:
    """
    Converts apparent magnitudes to absolute magnitudes
    """
    return apparent_mags -5*np.log10(cosmology.luminosity_distance(redshifts).value) - 25

def calculate_z_limit(cosmology, redshifts, absolute_mags, mag_limit) -> float | np.ndarray[float]:
    """
    Determines the mininum redshift value (zmin) where the galaxy would pass the faint magnitude
    limit. This depends on the classic formula m_a - m_b = 5log_10(d_a/d_b).

    If apparent_mag is larger than the magnitude limit then z will be z_max.
    If apparent_mag is less than the magnitude limit then z will be z_min.
    """
    distances = cosmology.luminosity_distance(redshifts)
    dist_to_z = interp1d(distances, redshifts, fill_value='extrapolate') # approximate inverse
    delta_m = mag_limit - absolute_mags
    distance_at_limit = (10**((np.array(delta_m) + 5)/5)) * u.pc
    z_lim = dist_to_z(distance_at_limit.to(u.Mpc).value)
    return z_lim


if __name__ == '__main__':
    ap_min, ap_max = 17., 19.65
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    area = calculate_area_of_rectangular_region(129*u.deg, 141*u.deg, -2*u.deg, 3*u.deg)
    fractional_area = (area.to(u.steradian)/(np.pi*4)).value

    gama = pd.read_csv('cut_9.dat', sep = '\s+')
    gama = gama[(gama['Rpetro'] < ap_max) & (gama['Rpetro'] > ap_min)]
    gama['abs'] = convert_ap_to_ab(gama['Rpetro'], gama['Z'], cosmo)
    gama['lumslog10'] = ab_mag_to_luminosity(gama['abs'])
    gama['zmins'] = calculate_z_limit(cosmo, gama['Z'], gama['abs'], ap_min)
    gama['zmaxs'] = calculate_z_limit(cosmo, gama['Z'], gama['abs'], ap_max)
    
    # working out max volumes 
    shell_volumes = cosmo.comoving_volume(gama['zmaxs']) - cosmo.comoving_volume(gama['zmins'])
    survey_volumes = shell_volumes * fractional_area
    max_volumes = survey_volumes.to(u.Mpc**3)
    gama['max_volumes'] = max_volumes.value
    #

    lum_bins = np.arange(gama['lumslog10'].min(), gama['lumslog10'].max() + 0.1, 0.1)
    lum_bins_mid = (lum_bins[:-1] + lum_bins[1:])/2
    estimator = []
    for i in range(len(lum_bins) - 1):
        local_bin = gama[(gama['lumslog10'] < lum_bins[i+1]) & (gama['lumslog10'] > lum_bins[i])]
        if len(local_bin) == 0:
            estimator.append(np.nan)
        else:
            v_maxs = cosmo.comoving_volume(local_bin['zmaxs']) - cosmo.comoving_volume(local_bin['zmins'])
            estimator.append(np.sum(1./(v_maxs.value*fractional_area)))
    estimator = np.array(estimator)
    plt.plot(lum_bins_mid, np.log10(estimator))
    plt.show()

    redshift_bins = np.arange(gama['Z'].min(), gama['Z'].max()+0.01, 0.01)
    mid_redshift_bins = (redshift_bins[:-1] + redshift_bins[1:])/2
    deltas = []
    for i in range(len(redshift_bins) -1 ):
        bin_volume = cosmo.comoving_volume(redshift_bins[i+1]) - cosmo.comoving_volume([redshift_bins[i]])
        bin_volume = bin_volume*fractional_area
        local_bin = gama[(gama['Z']>redshift_bins[i]) & (gama['Z'] < redshift_bins[i+1])]
        limit_ab = convert_ap_to_ab(ap_max, redshift_bins[i+1], cosmo)
        log10_lum_lim = ab_mag_to_luminosity(limit_ab)

        n_p = np.sum(estimator[lum_bins_mid > log10_lum_lim])

        delta = len(local_bin)/(n_p * bin_volume.value*fractional_area)
        deltas.append(delta)
    plt.plot(mid_redshift_bins, deltas)
    plt.show()
