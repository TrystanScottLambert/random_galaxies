"""
Module which implements the main algorithm for cloning a sample of galaxies.
"""

from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from astropy.units import Quantity
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
from geom_calcs import SurveyGeometries



def generate_clones(redshift_array: np.ndarray[float], n_copies: np.ndarray[int], z_ranges: np.ndarray[tuple[float]]) -> np.ndarray:
    """
    Optimized version of cloning galaxies within the redshift_array.
    """
    total_new_galaxies = np.sum(n_copies)

    # Preallocate array for the new galaxies
    new_galaxies = np.empty(total_new_galaxies, dtype=redshift_array.dtype)

    start_idx = 0
    for n, z_range in zip(n_copies, z_ranges):
        end_idx = start_idx + n
        new_galaxies[start_idx:end_idx] = np.random.uniform(z_range[0], z_range[1], n)
        start_idx = end_idx

    combined_array = np.concatenate((redshift_array, new_galaxies))
    return combined_array

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
    bins: np.ndarray[float]
    n_clones: int = 400 # default value used in Farrow+2015
    randoms: np.ndarray[float] = None

    def __post_init__(self) -> None:
        if np.max(self.magnitudes) > self.faint_mag_limit:
            raise ValueError(
                "There are magnitudes which are fainter than the faint magnitude limit."
            )
        if np.min(self.magnitudes) < self.bright_mag_limit:
            raise ValueError(
                "There are magnitudes which are brighter than the bright magnitude limit."
            )
        self.cosmo = self.geometry.cosmology  # I don't want to write this all the time
        self.mid_bins = np.array([(self.bins[i] + self.bins[i+1])/2 for i in range(len(self.bins) - 1)])
        self.n_g, _ = np.histogram(self.redshifts, bins = self.bins)

        #working out the zmin and zmax for each galaxy.
        self.z_maxs = self._calculate_z_limit(self.magnitudes, self.redshifts, self.faint_mag_limit)
        self.z_mins = self._calculate_z_limit(self.magnitudes, self.redshifts, self.bright_mag_limit)

        if self.randoms is None:
            # Intitial step (before we have delta(z) determined) set deltaz = 1
            # if delta is 1 then the number of copies is n_clones for all gals
            n_copies = self.n_clones * np.ones(len(self.redshifts))
            limits = np.array([self.z_mins, self.z_maxs]).T
            self.randoms = generate_clones(self.redshifts, n_copies.astype(int), limits)

    def _calculate_z_limit(
        self,
        apparent_mag: float | np.ndarray[float],
        redshift: float | np.ndarray[float],
        mag_limit: float,
    ) -> float | np.ndarray[float]:
        """
        Determines the mininum redshift value (zmin) where the galaxy would pass the faint magnitude
        limit. This depends on the classic formula m_a - m_b = 5log_10(d_a/d_b).

        If apparent_mag is larger than the magnitude limit then z will be z_max.
        If apparent_mag is less than the magnitude limit then z will be z_min.
        """
        distance = self.cosmo.luminosity_distance(redshift)
        delta_m = mag_limit - apparent_mag
        distance_at_limit = distance * 10 ** (delta_m / 5)
        z_min = z_at_value(self.cosmo.luminosity_distance, distance_at_limit)
        return z_min.value

    @property
    def max_volumes(self) -> np.ndarray[Quantity[u.Mpc**3]]:
        """
        Calculates the maximum possible volumes that are available to each galaxy.
        """
        max_volumes = self.geometry.calculate_survey_volume(self.z_mins, self.z_maxs)
        return max_volumes

    @property
    def delta(self) -> np.ndarray[float]:
        """
        Generates the overdensity in redshift of every galaxy at each redshift.
        """
        n_r, _ = np.histogram(self.randoms, bins = self.bins)
        delta_vals = self.n_clones * (self.n_g/n_r)
        delta_vals[delta_vals==np.nan] = 0
        delta_func = interp1d(self.mid_bins, delta_vals, fill_value='extrapolate')
        return delta_func

    @property
    def volume_dcs(self) -> np.ndarray[Quantity[u.Mpc**3]]:
        """
        Approximating the integrals of the overdensities. 
        """
        # Create a grid of integrand values that we will approximate.
        z_grid = np.linspace(np.min(self.z_mins), np.max(self.z_maxs), 5000)
        dv_dz = self.cosmo.differential_comoving_volume(z_grid).value
        delta_vals = self.delta(z_grid)
        integrand_vals = delta_vals * dv_dz

        #approximate integral
        cummulative_integral = np.cumsum(integrand_vals) * np.gradient(z_grid) # trapezoid method
        z_grid_interp = interp1d(z_grid, cummulative_integral, kind='linear', bounds_error=True)

        volumes = z_grid_interp(self.z_maxs) - z_grid_interp(self.z_mins)
        volumes *= self.geometry.area.value
        volumes = volumes * (u.Mpc**3 / u.steradian) * self.geometry.area.unit
        return volumes.to(u.Mpc**3)


    @property
    def number_copies(self) -> np.ndarray[int]:
        """
        Generates the number of times each galaxy is copied.
        """
        vals = self.n_clones * (self.max_volumes/self.volume_dcs)
        print(vals)
        return vals.astype(int) # has to be an int

    @property
    def clones(self) -> np.ndarray[float]:
        """
        Takes every galaxy at redshift i in redshift array and randomly distributes n_clones around
        that redshift. This is the simple implementation used by Cole+2011.

        returns the redshift array with the real galaxies + the cloned galaxies. 
        """
        zs = np.array([self.z_mins, self.z_maxs]).T
        return generate_clones(self.redshifts, self.number_copies, zs)

    @property
    def windowed_clones(self) -> np.ndarray[float]:
        """
        Takes every galaxy at redshift i and creates n number of clones around that redshift,
        applying the gaussian window described in Farrow+2015 (extension of Cole+2011).
        """


if __name__ == "__main__":
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    area = 110 * (u.deg**2)
    test_zs = np.array([0.1, 0.11, 0.2, 0.12, 0.4])  # redshifts
    test_mags = np.array([19, 19.2, 18.4, 18.3, 19.7])  # magntidues

    survey = SurveyGeometries(cosmo, area)

    test_array = np.random.normal(0.2, 0.01, 1000)
    test_magnitudes = np.random.normal(18, 0.2, len(test_array))
    z_ranges = np.array([(z-0.005, z+0.005) for z in test_array])
    ns = np.ones(len(test_array)) * 400

    test_survey = Survey(survey, test_array, test_magnitudes, 19.8, 17, np.arange(0, 0.4, 0.001))
    dc_volumes = test_survey.volume_dcs
    new_redshifts = test_survey.clones
