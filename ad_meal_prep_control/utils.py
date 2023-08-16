from dataclasses import dataclass
import numpy as np


@dataclass
class CHP:
    """
    Combined heat and power plant.

    Attributes:
        max_power:
            Maximum power of the plant in kW.
        thermal_efficiency:
            Thermal efficiency of the plant as percentage.
    """

    max_power: float = 100.0  # kW
    thermal_efficiency: float = 42.0  # %
    lhv = 50.01  # CH4: MJ/kg
    R_S = 518.4  # J/(kg*K)

    def __post_init__(self):
        # Convert lhv to kJ/kg
        self.lhv *= 1000.0

        # Convert thermal efficiency to decimal number
        self.thermal_efficiency /= 100.0

    def compute_vol_flow_rate(self, load: np.ndarray, press: float, temp: float):
        """
        Computes the actual volume flow rate of CH4 from the gas storage in m^3/s.

        Parameters:
            load:
                Array consisting of the plant utilization at different time steps in percent.
            press:
                Pressure of the gas storage in bars.
            temp:
                Temperature of the gas storage in K.

        Returns:
            Array of the volume flow rates in m^3/s.
        """
        mass_flow_rate = (
            self.max_power * load / (100.0 * self.lhv * self.thermal_efficiency)
        )
        vol_flow_rate = mass_flow_rate * self.R_S * temp / (press * 1.0e5)

        return vol_flow_rate
