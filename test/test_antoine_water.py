import unittest
from ad_meal_prep_control.antoine_water import vapour_pressure_h2o

# The values for testing were retrieved from this url: https://en.wikipedia.org/wiki/Vapour_pressure_of_water


class TestVapourPressure(unittest.TestCase):
    """Testing if values are within 0.01% deviation from true vapour pressure."""

    def test_vapour_pressure(self):
        true_pressures = {
            0: 0.6113,
            20.0: 2.3388,
            40.0: 7.3814,
            60.0: 19.9320,
            80.0: 47.3730,
            100.0: 101.3200,
        }
        for temp, true_press in true_pressures.items():
            self.assertAlmostEqual(
                vapour_pressure_h2o(temp + 273.15),
                true_press * 1.0e-2,
                delta=true_press / 10000.0,
            )
