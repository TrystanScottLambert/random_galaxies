"""
Module which implements the main algorithm for cloning a sample of galaxies.
"""

from warnings import warn
from dataclasses import dataclass

import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.units import Quantity
from scipy.integrate import cumulative_trapezoid, quad
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from geom_calcs import SurveyGeometries, calculate_area_of_rectangular_region


def test_delta(z):
    return 1

def _term(redshift, a_value, exponent, z_p):
    """
    helper function for the k_correction. Within the summation in eq. 8 of Robotham+2011.
    """
    return a_value * ((redshift-z_p)**exponent)

def k_correction(redshifts):
    """
    Calculates the k_correction via eq. 8 of Robotham+2011

    This also includes the e correction so what is returned is the (k+e)(z)
    """
    Q_ZREF = 1.75 # when ZREF = 0 
    A = np.array([[0.2085], [1.0226], [0.5237], [3.5902], [2.3843]])
    Z_REF = 0
    Z_P = 0.2
    val = np.sum([_term(redshifts, a, i, Z_P) for i, a in enumerate(A)], axis=0)
    return val + Q_ZREF*(redshifts - Z_REF)


def generate_clones(redshift_array: np.ndarray[float], n_copies: np.ndarray[int], z_ranges: np.ndarray[tuple[float]]) -> np.ndarray:
    """
    Optimized version of cloning galaxies within the redshift_array.
    """
    total_new_galaxies = np.sum(n_copies)
    print('This is the number:',total_new_galaxies/len(redshift_array))

    # Preallocate array for the new galaxies
    new_galaxies = np.empty(total_new_galaxies, dtype=redshift_array.dtype)

    start_idx = 0
    for n, z_range in zip(n_copies, z_ranges):
        end_idx = start_idx + n
        new_galaxies[start_idx:end_idx] = np.random.uniform(z_range[0], z_range[1], n)
        start_idx = end_idx
    plt.hist(new_galaxies, bins=100)
    plt.show()
    return new_galaxies

@dataclass
class Survey:
    """
    Class which manages the observed survey including meta data such as geometry and cosmology
    """
    geometry: SurveyGeometries
    redshifts: np.ndarray[float]
    magnitudes: np.ndarray[float]
    faint_mag_limit: float
    bright_mag_limit: float
    binwidth: float
    n_clones: int = 100 # default value used in Farrow+2015
    randoms: np.ndarray[float] = None
    k_corrections: np.ndarray[float] = None

    def __post_init__(self) -> None:
        if np.max(self.magnitudes) > self.faint_mag_limit:
            raise ValueError(
                "There are magnitudes which are fainter than the faint magnitude limit."
            )
        if np.min(self.magnitudes) < self.bright_mag_limit:
            warn(
                "There are magnitudes which are brighter than the bright magnitude limit. The zmin for these galaxies will be set to 0."
                 )
        self.cosmo = self.geometry.cosmology  # I don't want to write this all the time

        self.absolute_mags = self.magnitudes -5*np.log10(self.cosmo.luminosity_distance(self.redshifts).value) - 25
        if self.k_corrections is not None:
            self.absolute_mags = self.absolute_mags - self.k_corrections

        self.bins = np.arange(np.min(self.redshifts), np.max(self.redshifts) + self.binwidth, self.binwidth)
        self.mid_bins = (self.bins[:-1] + self.bins[1:])/2
        self.n_g, _ = np.histogram(self.redshifts, bins = self.bins, density=True)

        #working out the zmin and zmax for each galaxy.
        self.z_maxs = self._calculate_z_limit(self.faint_mag_limit)
        self.z_mins = self._calculate_z_limit(self.bright_mag_limit)
        self.z_mins[self.z_mins < 0] = 0

        # get the max volumes for each galxy. This won't change
        self.max_volumes = self.geometry.calculate_survey_volume(self.z_mins, self.z_maxs)
    
        if self.randoms is None:
            # Intitial step (before we have delta(z) determined) set deltaz = 1
            # if delta is 1 then the number of copies is n_clones for all gals
            n_copies = self.n_clones * np.ones(len(self.redshifts))
            limits = np.array([self.z_mins, self.z_maxs]).T
            self.randoms = generate_clones(self.redshifts, n_copies.astype(int), limits)

    def _calculate_z_limit(self, mag_limit) -> float | np.ndarray[float]:
        """
        Determines the mininum redshift value (zmin) where the galaxy would pass the faint magnitude
        limit. This depends on the classic formula m_a - m_b = 5log_10(d_a/d_b).

        If apparent_mag is larger than the magnitude limit then z will be z_max.
        If apparent_mag is less than the magnitude limit then z will be z_min.
        """
        distances = self.cosmo.luminosity_distance(self.redshifts)
        dist_to_z = interp1d(distances, self.redshifts, fill_value='extrapolate') # approximate inverse

        delta_m = mag_limit - self.absolute_mags
        distance_at_limit = (10**((delta_m + 5)/5)) * u.pc

        z_lim = dist_to_z(distance_at_limit.to(u.Mpc).value)
        return z_lim

    def delta(self) -> np.ndarray[float]:
        """
        Generates the overdensity in redshift of every galaxy at each redshift.
        """
        n_r, _ = np.histogram(self.randoms, bins = self.bins, density=True)
        plt.plot(self.mid_bins, self.n_g/n_r)
        plt.axhline(0, color='r')
        plt.axhline(1, color='k')
        plt.show()
        delta_vals =  (self.n_g/n_r)# * self.n_clones
        delta_func = interp1d(self.mid_bins, delta_vals, fill_value='extrapolate')
        return delta_func

    @property
    def volume_dcs(self) -> np.ndarray[Quantity[u.Mpc**3]]:
        """
        Approximating the integrals of the overdensities. 
        """
        # Create a grid of integrand values that we will approximate.
        z_grid = np.linspace(0, 2, 500)
        dv_dz = self.cosmo.differential_comoving_volume(z_grid) #Mpc3/steradian
        dv_dz *= self.geometry.area
        delta_vals = test_delta(z_grid)
        integrand_vals = delta_vals * dv_dz.to(u.Mpc**3).value

        integrated_vals = [quad()]

        #approximate integral
        cummulative_integral = np.cumsum(integrand_vals) * np.gradient(z_grid) # cummulative trapezoid method
        cum_trap = cumulative_trapezoid(integrand_vals, z_grid, initial=0)
        print(len(cum_trap))
        print(len(cummulative_integral))
        z_grid_interp = interp1d(z_grid, cum_trap, fill_value='extrapolate')
        volumes = z_grid_interp(self.z_maxs) - z_grid_interp(self.z_mins)

        for volume, max_volume, zmin, zmax in zip(volumes, self.max_volumes, self.z_mins, self.z_maxs):
            print(volume, max_volume, zmin, zmax)
        return volumes.to(u.Mpc**3)

    @property
    def number_copies(self) -> np.ndarray[int]:
        """
        Generates the number of times each galaxy is copied.
        """
        vals = self.n_clones * (self.max_volumes/self.volume_dcs)
        vals[vals < 0] = 0
        plt.scatter(np.arange(len(vals)),vals)
        plt.show()

        return vals.astype(int) # has to be an int

    @property
    def clones(self) -> np.ndarray[float]:
        """
        Takes every galaxy at redshift i in redshift array and randomly distributes n_clones around
        that redshift. This is the simple implementation used by Cole+2011.

        returns the redshift array with the real galaxies + the cloned galaxies. 
        """
        z_limits = np.array([self.z_mins, self.z_maxs]).T
        return generate_clones(self.redshifts, self.number_copies, z_limits)

    @property
    def windowed_clones(self) -> np.ndarray[float]:
        """
        Takes every galaxy at redshift i and creates n number of clones around that redshift,
        applying the gaussian window described in Farrow+2015 (extension of Cole+2011).
        """

def generate_random_cat(survey: Survey, method: str = 'un-windowed') -> np.ndarray[float]:
    """
    Creates the random catalog by other using the windowed method or the non window method.
    """
    for i in range(2):
        if method == 'un-windowed':
            clones = survey.clones
            print(f'This should be 100: {len(clones)/len(survey.absolute_mags)}')
            survey.randoms = clones
        elif method == 'windowed':
            clones = survey.windowed_clones
            survey.randoms = clones
        else:
            raise ValueError('method must be either "un-windowed" or "windowed".')
        print(f'{i+1} iterations')
    return clones


if __name__ == "__main__":
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    test_zs = np.array([0.1, 0.11, 0.2, 0.12, 0.4])  # redshifts
    test_mags = np.array([19, 19.2, 18.4, 18.3, 19.7])  # magntidues


    area = calculate_area_of_rectangular_region(129*u.deg, 141*u.deg, -2*u.deg, 3*u.deg)
    survey = SurveyGeometries(cosmo, area)

    test_array = np.random.normal(0.2, 0.01, 1000)
    test_magnitudes = np.random.normal(18, 0.2, len(test_array))
    z_ranges = np.array([(z-0.005, z+0.005) for z in test_array])
    ns = np.ones(len(test_array)) * 400

    # Testing with actual GAMA data
    gama_z, gama_mag = np.loadtxt('cut_9.dat', usecols=(-2, -1), unpack=True, skiprows=1)
    gama_k_corrections = k_correction(gama_z)

    gama = Survey(survey, gama_z, gama_mag, 19.8, 17, 0.005, k_corrections=gama_k_corrections)
    clones = generate_random_cat(gama)
    bins = np.arange(0, 0.6, 0.001)
    plt.hist(clones, histtype='step', density=True, bins=bins, label='our randoms')
    published_clones_unwindoes = np.loadtxt('randoms_unpublished.csv', skiprows=1)
    plt.hist(published_clones_unwindoes, density=True, histtype='step', bins=bins, label='non-windowed')
    infile= '../../my_tools/group_finders/FoFR/gen_ran_out.randoms.csv'
    published_clones = np.loadtxt(infile, skiprows=1)
    plt.hist(published_clones, bins=bins, density=True, histtype='step', label='R')
    plt.hist(gama_z, bins=bins, density=True, histtype='step', label='GAMA-9 region')
    plt.xlabel('redshift')
    plt.legend()
    plt.show()
