"""
Helper functions for dealing with the geometry such as volume and area.
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
class SurveyGeometries:
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
        return self.cosmology.comoving_volume(z_max) - self.cosmology.comoving_volume(z_min)

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
