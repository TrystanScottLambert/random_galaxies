"""
Helper functions for calculating the volume around galaxies in a redshift survey.
"""

from dataclasses import dataclass
import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.units import Quantity

SKY_AREA = 4 * np.pi * u.steradian

@dataclass
class SurveyCosmology:
    """
    survey cosmology class which can be used to calculate volumes.
    """

    cosmology: FlatLambdaCDM
    area: Quantity

    def __post_init__(self):
        if not self.area.unit.is_equivalent(u.steradian):
            raise ValueError("Area must have area units.")

    def calculate_shell_volume(
        self, z_min: float | np.ndarray[float], z_max: float | np.ndarray[float]
    ) -> Quantity[u.Mpc**3]:
        """
        Calculates the simple comoving shell volume from z_min to z_max and returns the volume in
        Mpc^3
        """
        return self.cosmology.comoving_volume(z_max) - self.cosmology.comoving_volume(
            z_min
        )

    def calculate_survey_volume(
        self, z_min: float | np.ndarray[float], z_max: float | np.ndarray[float]
    ) -> Quantity[u.Mpc**3]:
        """
        Works out the volume that the actual survey footprint would encompass.
        """
        area_radians = self.area.to(u.steradian)  # make sure area is in steradians
        percentage = (
            area_radians / SKY_AREA
        )  # what percentage of the sky does the survey cover
        return percentage * self.calculate_shell_volume(z_min, z_max)


if __name__ == "__main__":
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    area = 110 * u.deg**2
    sc = SurveyCosmology(cosmo, area)
    shell = sc.calculate_shell_volume(0.3, 0.4)
    print(shell)
    survey = sc.calculate_survey_volume(1, 2)
