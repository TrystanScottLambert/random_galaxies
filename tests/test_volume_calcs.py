"""
Tests for the volume_cals.py module.
"""

import unittest
import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from volume_calcs import SurveyCosmology


class TestSurveyCosmology(unittest.TestCase):
    """
    Testing the SurveyCosmology class and methods.
    """
    @classmethod
    def setUpClass(cls):
        """Default cosmology for all tests, with Om0 = 0.3"""
        cls.default_cosmology = FlatLambdaCDM(H0=100, Om0=0.3)

    def setUp(self):
        """Default survey with 100 deg^2 area"""
        self.default_survey = SurveyCosmology(
            cosmology=self.default_cosmology, area=100 * u.deg**2
        )

    def test_area_unit_validation(self):
        """Testing that non-area quantities return a value error."""
        with self.assertRaises(ValueError):
            SurveyCosmology(self.default_cosmology, area=10 * u.m)

    def test_calculate_shell_volume_single(self):
        """Testing that the shell volume works for single float values."""
        z_min = 0.1
        z_max = 0.2
        volume = self.default_survey.calculate_shell_volume(z_min, z_max)
        self.assertTrue(volume.unit.is_equivalent(u.Mpc**3))
        self.assertGreater(volume.value, 0)

    def test_calculate_shell_volume_array(self):
        """Testing that shell volume will work for an array of values."""
        z_min = np.array([0.1, 0.2, 0.3])
        z_max = np.array([0.2, 0.3, 0.4])
        volume = self.default_survey.calculate_shell_volume(z_min, z_max)
        self.assertTrue(volume.unit.is_equivalent(u.Mpc**3))
        self.assertTrue(np.all(volume.value > 0))

    def test_calculate_survey_volume_single(self):
        """Testing that survey volumes work for a single value."""
        z_min = 0.1
        z_max = 0.2
        volume = self.default_survey.calculate_survey_volume(z_min, z_max)
        self.assertTrue(volume.unit.is_equivalent(u.Mpc**3))
        self.assertGreater(volume.value, 0)

    def test_calculate_survey_volume_array(self):
        """Testing that survey volumes work for an array of values."""
        z_min = np.array([0.1, 0.2, 0.3])
        z_max = np.array([0.2, 0.3, 0.4])
        volume = self.default_survey.calculate_survey_volume(z_min, z_max)
        self.assertTrue(volume.unit.is_equivalent(u.Mpc**3))
        self.assertTrue(np.all(volume.value > 0))

    def test_calculate_shell_volume_known_values(self):
        """Comparing shell volumes to the celestial package in R. (Aaron's version.)"""
        survey = SurveyCosmology(self.default_cosmology, area=4 * np.pi * u.steradian)

        # Known volumes
        z_min, z_max = 0.3, 0.4
        expected_volume = 2.918226 * (u.Gpc**3)
        volume = survey.calculate_shell_volume(z_min, z_max)
        print(volume)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)

        z_min, z_max = 3, 4
        expected_volume = 160.8046 * (u.Gpc**3)
        volume = survey.calculate_shell_volume(z_min, z_max)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)

    def test_calculate_survey_volume_known_values(self):
        """
        Comparing survey volumes to the celestial package in R. (Aaron's version) by
        expecting half of the shell volume when half the sky area is passed (2pi steradian)
        """
        survey = SurveyCosmology(self.default_cosmology, area=2 * np.pi * u.steradian)

        # Known scaled volumes
        z_min, z_max = 0.3, 0.4
        expected_volume = (2.918226/ 2) * (u.Gpc**3)
        volume = survey.calculate_survey_volume(z_min, z_max)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)

        z_min, z_max = 3, 4
        expected_volume = (160.8046 / 2) * (u.Gpc**3)
        volume = survey.calculate_survey_volume(z_min, z_max)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)


if __name__ == "__main__":
    unittest.main()
