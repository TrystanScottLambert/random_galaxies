"""
Tests for the volume_cals.py module.
"""

import unittest
import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from geom_calcs import SurveyGeometries, calculate_area_of_rectangular_region


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
        self.default_survey = SurveyGeometries(
            cosmology=self.default_cosmology, area=100 * u.deg**2
        )

    def test_area_unit_validation(self):
        """Testing that non-area quantities return a value error."""
        with self.assertRaises(ValueError):
            SurveyGeometries(self.default_cosmology, area=10 * u.m)

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
        survey = SurveyGeometries(self.default_cosmology, area=4 * np.pi * u.steradian)

        # Known volumes
        z_min, z_max = 0.3, 0.4
        expected_volume = 2.918226 * (u.Gpc**3)
        volume = survey.calculate_shell_volume(z_min, z_max)
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
        survey = SurveyGeometries(self.default_cosmology, area=2 * np.pi * u.steradian)

        # Known scaled volumes
        z_min, z_max = 0.3, 0.4
        expected_volume = (2.918226/ 2) * (u.Gpc**3)
        volume = survey.calculate_survey_volume(z_min, z_max)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)

        z_min, z_max = 3, 4
        expected_volume = (160.8046 / 2) * (u.Gpc**3)
        volume = survey.calculate_survey_volume(z_min, z_max)
        self.assertAlmostEqual(volume.to(u.Gpc**3).value, expected_volume.value, places=3)


class TestCalculateAreaOfRectangularRegion(unittest.TestCase):
    """
    Testing the calculate_area_of_rectangular_region function.
    """
    def test_single_values(self):
        """Test the function with single values."""
        ra1 = 129 * u.deg
        ra2 = 141 * u.deg
        dec1 = -2 * u.deg
        dec2 = 3 * u.deg

        expected_area = 59.97867933 * u.deg**2
        result = calculate_area_of_rectangular_region(ra1, ra2, dec1, dec2)
        self.assertAlmostEqual(result.value, expected_area.value, places=5)
        self.assertEqual(result.unit, expected_area.unit)

    def test_array_values(self):
        """Test the function with arrays of values."""
        ra1 = np.array([129, 129]) * u.deg
        ra2 = np.array([141, 141]) * u.deg
        dec1 = np.array([-2, -2]) * u.deg
        dec2 = np.array([3, 3]) * u.deg

        expected_area = np.array([59.97867933, 59.97867933]) * u.deg**2
        result = calculate_area_of_rectangular_region(ra1, ra2, dec1, dec2)
        np.testing.assert_almost_equal(result.value, expected_area.value, decimal=5)
        self.assertEqual(result.unit, expected_area.unit)

    def test_invalid_units(self):
        """Test error handling for invalid units."""
        ra1 = 129 * u.m  # Invalid unit
        ra2 = 141 * u.deg
        dec1 = -2 * u.deg
        dec2 = 3 * u.deg

        with self.assertRaises(ValueError) as context:
            calculate_area_of_rectangular_region(ra1, ra2, dec1, dec2)
        self.assertIn("Arguments must be in angular units", str(context.exception))

    def test_missing_quantity(self):
        """Test error handling for missing Quantity (e.g., raw arrays)."""
        ra1 = np.array([129, 129])  # Not a Quantity
        ra2 = np.array([141, 141]) * u.deg
        dec1 = np.array([-2, -2]) * u.deg
        dec2 = np.array([3, 3]) * u.deg

        with self.assertRaises(ValueError) as context:
            calculate_area_of_rectangular_region(ra1, ra2, dec1, dec2)
        self.assertIn("Arrays must be passed as angular Quantities", str(context.exception))


if __name__ == "__main__":
    unittest.main()
