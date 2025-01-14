"""
Module which implements the main algorithm for cloning a sample of galaxies.
"""

from dataclasses import dataclass
import numpy as np
from astropy.units import Quantity
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
from geom_calcs import SurveyGeometries


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

    def __post_init__(self) -> None:
        if np.min(self.magnitudes) > self.faint_mag_limit:
            raise ValueError(
                "There are magnitudes which are fainter than the faint magnitude limit."
            )
        if np.max(self.magnitudes) < self.bright_mag_limit:
            raise ValueError(
                "There are magnitudes which are brighter than the bright magnitude limit."
            )
        self.cosmo = self.geometry.cosmology  # I don't want to write this all the time
        

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
    def z_limits(self) -> tuple[np.ndarray[float], np.ndarray[float]]:
        """
        working out the zmin and zmaxs for each galaxy
        """
        z_mins = self._calculate_z_limit(self.magnitudes, self.redshifts, self.faint_mag_limit)
        z_maxs = self._calculate_z_limit(self.magnitudes, self.redshifts, self.bright_mag_limit)
        return z_mins, z_maxs

    @property
    def max_volumes(self) -> np.ndarray[Quantity[u.Mpc**3]]:
        """
        Calculates the maximum possible volumes that are available to each galaxy.
        """
        z_mins, z_maxs = self.z_limits
        max_volumes = self.geometry.calculate_survey_volume(z_mins, z_maxs)
        return max_volumes



if __name__ == "__main__":
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    area = 110 * (u.deg**2)
    test_zs = [0.1, 0.11, 0.2, 0.12, 0.4]  # redshifts
    test_mags = [19, 19.2, 18.4, 18.3, 19.7]  # magntidues
    faint_mag_limit = 19.8
    bright_mag_limit = 17.0

    N_CLONES = 400
    survey = SurveyGeometries(cosmo, area)
