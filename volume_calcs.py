"""
Helper functions for calculating the volume around galaxies in a redshift survey.
"""

from dataclasses import dataclass
import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.units import Quantity

SKY_AREA = 4 * np.pi * u.steradian


def calculate_area_of_rectangular_region(
    alpha_1: Quantity, alpha_2: Quantity, dec_1: Quantity, dec_2: Quantity
) -> Quantity[u.deg**2]:
    """
    Calculates the on sky area of a rectangular region between alpha_1 and alpha_2 and dec_1
    and dec_2.
    """
    # check that the units are all degrees:
    for angle in [alpha_1, alpha_2, dec_1, dec_2]:
        if isinstance(angle, np.ndarray) and not isinstance(angle, Quantity):
            raise ValueError("Arrays must be passed as angular Quantities (*u.deg).")
        if not angle.unit.is_equivalent(u.deg):
            raise ValueError("Arguments must be in angular units (e.g. degrees).")

    area_in_steradians = (
        (alpha_1.to(u.rad).value - alpha_2.to(u.rad).value)
        * (np.sin(dec_1.to(u.rad).value) - np.sin(dec_2.to(u.rad).value))
        * (u.steradian)
    )
    return area_in_steradians.to(u.deg**2)


@dataclass
class SurveyCosmology:
    """
    Survey cosmology class which can be used to calculate volumes.
    """

    cosmology: FlatLambdaCDM
    area: Quantity

    def __post_init__(self) -> None:
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
    survey = sc.calculate_survey_volume(1, 2)
    x = calculate_area_of_rectangular_region(
        129 * u.deg, 141 * u.deg, -2 * u.deg, 3 * u.deg
    )
    alpha_1s = np.ones(2) * 129
    alpha_2s = np.ones(2) * 141
    dec_1s = np.ones(2) * -2
    dec_2s = np.ones(2) * 3
    y = calculate_area_of_rectangular_region(
        alpha_1s * u.deg, alpha_2s * u.deg, dec_1s * u.deg, dec_2s * u.deg
    )
