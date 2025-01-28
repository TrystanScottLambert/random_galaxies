"""
Script that will generate galaxies using the pure Cole method.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import curve_fit
import pylab as plt
import pandas as pd
import tqdm

from geom_calcs import calculate_area_of_rectangular_region

def double_power_law(logL, logL_star, log_phi_star, alpha, beta):
    L = 10**np.array(logL, dtype=float)  # Ensure L is a float
    L_star = 10**float(logL_star)        # Ensure L_star is a float
    phi_star = 10**float(log_phi_star)   # Ensure phi_star is a float
    return np.log10(phi_star) + alpha * np.log10(L / L_star) - np.log10(1 + (L / L_star)**beta)

def double_power_law_not_log(logL, logL_star, log_phi_star, alpha, beta):
    return 10**(double_power_law(logL, logL_star, log_phi_star, alpha, beta))

def fit_double_power_law(logL, log_phi, initial_guess=None, maxfev=10000):
    """
    Fits a double power-law to the input data while masking invalid values.

    Parameters:
    -----------
    logL : array-like
        The log10 of luminosity values (independent variable).
    log_phi : array-like
        The log10 of phi values (dependent variable).
    initial_guess : list or None, optional
        Initial guesses for the parameters [logL_star, log_phi_star, alpha, beta].
    maxfev : int, optional
        Maximum number of function evaluations for the fitting process.

    Returns:
    --------
    popt : array
        Fitted parameters [logL_star, log_phi_star, alpha, beta].
    pcov : 2D array
        Covariance matrix of the fitted parameters.
    """
    # Mask invalid values (NaN, inf)
    valid_mask = np.isfinite(log_phi)
    logL_valid = logL[valid_mask]
    log_phi_valid = log_phi[valid_mask]

    if initial_guess is None:
        initial_guess = [10, -2, -1, 2]

    popt, pcov = curve_fit(double_power_law, logL_valid, log_phi_valid, p0=initial_guess, maxfev=maxfev)
    return popt, pcov


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
    return apparent_mags -5*np.log10(cosmology.luminosity_distance(redshifts).to(u.pc).value) + 5

def calculate_z_limit(cosmology, absolute_mags, mag_limit) -> float | np.ndarray[float]:
    """
    Determines the mininum redshift value (zmin) where the galaxy would pass the faint magnitude
    limit. This depends on the classic formula m_a - m_b = 5log_10(d_a/d_b).

    If apparent_mag is larger than the magnitude limit then z will be z_max.
    If apparent_mag is less than the magnitude limit then z will be z_min.
    """

    delta_m = mag_limit - absolute_mags
    distance_at_limit = (10**((np.array(delta_m) + 5)/5)) * u.pc

    max_redshift = z_at_value(cosmology.luminosity_distance, np.max(distance_at_limit))
    zs = np.linspace(0, max_redshift+0.1, 10000)
    distances = cosmology.luminosity_distance(zs)
    dist_to_z = interp1d(distances.to(u.Mpc).value, zs)
    z_lim = np.array(dist_to_z(distance_at_limit.to(u.Mpc).value))
    return z_lim


def estimate_lf(data: pd.DataFrame, volume_label: str, luminosity_label: str, lum_bins: np.ndarray) -> np.ndarray:
    """
    estimates the LF via a 1/V estimator given the volume label.
    """
    estimator = []
    for i in range(len(lum_bins) - 1):
        local_bin = data[(data[luminosity_label] < lum_bins[i+1]) & (data[luminosity_label] >= lum_bins[i])]
        if len(local_bin) == 0:
            estimator.append(np.nan)
        else:
            estimator.append(np.sum(1./(local_bin[volume_label])))
    return np.array(estimator)


if __name__ == '__main__':
    ap_min, ap_max = 17., 19.8
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    area = calculate_area_of_rectangular_region(129*u.deg, 141*u.deg, -2*u.deg, 3*u.deg)
    fractional_area = (area.to(u.steradian)/(np.pi*4)).value

    gama = pd.read_csv('cut_9.dat', sep = '\s+')
    #gama = gama[(gama['Rpetro'] < ap_max) & (gama['Rpetro'] > ap_min)]
    gama['abs'] = convert_ap_to_ab(gama['Rpetro'], gama['Z'], cosmo)
    gama['lumslog10'] = ab_mag_to_luminosity(gama['abs'])
    gama['zmins'] = np.zeros(len(gama))# calculate_z_limit(cosmo, gama['Z'], gama['abs'], ap_min)
    gama['zmaxs'] = calculate_z_limit(cosmo, gama['abs'], ap_max)

    # working out max volumes
    shell_volumes = cosmo.comoving_volume(gama['zmaxs']) - cosmo.comoving_volume(gama['zmins'])
    survey_volumes = shell_volumes * fractional_area
    max_volumes = survey_volumes.to(u.Mpc**3)
    gama['max_volumes'] = max_volumes.value
    ###

    #Luminosity bins
    lum_binwidth = 0.05
    lum_bins = np.arange(gama['lumslog10'].min(), gama['lumslog10'].max() + lum_binwidth, lum_binwidth)
    lum_bins_mid = (lum_bins[:-1] + lum_bins[1:])/2

    estimator = estimate_lf(gama, 'max_volumes', 'lumslog10', lum_bins)
    estimator *= 100
    popt, pcov = fit_double_power_law(lum_bins_mid, np.log10(estimator))
    log_l_star, log_phi_star, alpha,  beta = popt

    plt.step(lum_bins_mid, np.log10(estimator), where='mid')
    plt.plot(lum_bins_mid, double_power_law(lum_bins_mid, *popt))
    plt.show()

    # Redshift bins
    redshift_bin_width = 0.005
    redshift_bins = np.arange(gama['Z'].min(), gama['Z'].max()+redshift_bin_width, redshift_bin_width)
    mid_redshift_bins = (redshift_bins[:-1] + redshift_bins[1:])/2
    for j in range(2):
        deltas = []
        bin_volumes = []
        for i in range(len(redshift_bins) -1 ):
            bin_volume = cosmo.comoving_volume(redshift_bins[i+1]) - cosmo.comoving_volume([redshift_bins[i]])
            bin_volume = bin_volume * fractional_area

            local_bin = gama[(gama['Z']>=redshift_bins[i]) & (gama['Z'] < redshift_bins[i+1])]
            limit_ab = convert_ap_to_ab(ap_max, redshift_bins[i], cosmo)
            log10_lum_lim = ab_mag_to_luminosity(limit_ab)

            n_p, _ = quad(double_power_law_not_log, log10_lum_lim, np.inf, args=(log_l_star, log_phi_star, alpha, beta))

            delta = len(local_bin)/(n_p * bin_volume.value)
            deltas.append(delta)
            bin_volumes.append(bin_volume.value)
        deltas = np.array(deltas)
        deltas = np.ones(len(deltas))


        delta_func = interp1d(mid_redshift_bins, deltas.reshape(len(deltas)), fill_value=0, bounds_error=False)

        def integrand(redshift):
            volume_element = cosmo.differential_comoving_volume(redshift) * area.to(u.steradian)
            return delta_func(redshift) * volume_element.value

        volumes = np.array(bin_volumes)

        zs = np.linspace(0, 0.6, 1000)
        plt.plot(zs, delta_func(zs))
        plt.show()

        plt.plot(zs, quad(integrand, 0, zs), label='integral')
        #plt.plot(zs, cosmo.comoving_volume(zs), label='astropy')
        plt.legend()
        plt.show()

        v_dc_maxs = []
        for zmin, zmax in zip(gama['zmins'], gama['zmaxs']):
            cut = np.where((mid_redshift_bins < zmax))[0]
            v_dc_maxs.append(np.sum(deltas[cut] * volumes[cut]))
            print((np.sum(deltas[cut] * volumes[cut]), quad(integrand, 0, zmax)[0]))


        gama['v_dc_maxs'] = np.array(v_dc_maxs)
        print(gama['v_dc_maxs'])
        estimator = estimate_lf(gama, 'v_dc_maxs', 'lumslog10', lum_bins)
        popt, pcov = fit_double_power_law(lum_bins_mid[5:], np.log10(estimator[5:]))
        log_l_star, log_phi_star, alpha,  beta = popt
        plt.step(lum_bins_mid, np.log10(estimator), where='mid')
        plt.plot(lum_bins_mid, double_power_law(lum_bins_mid, *popt))
        plt.show()
        plt.plot(mid_redshift_bins, deltas, label=f'{j}')
    plt.legend()
    plt.show()

    number_of_clones = 400/deltas
