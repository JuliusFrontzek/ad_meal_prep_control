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
            Thermal efficiency of the plant as a percentage.
    """

    max_power: float = 100.0
    thermal_efficiency: float = 0.36
    lhv = 50.01  # CH4: MJ/kg
    R_S = 518.4  # J/(kg*K)

    def __post_init__(self):
        # Convert lhv to kJ/kg
        self.lhv *= 1000.0

        # Convert thermal efficiency to decimal number
        assert (
            self.thermal_efficiency >= 0 and self.thermal_efficiency <= 1.0
        ), "Thermal efficiency must be between 0 and 1."

    def ch4_vol_flow_rate(self, load: np.ndarray, press: float, temp: float):
        """
        Computes the actual volume flow rate of CH4 from the gas storage in m^3/d.

        Parameters:
            load:
                Array consisting of the plant utilization at different time steps
                as a decimal value between 0 and 1.
            press:
                Pressure of the gas storage in bars.
            temp:
                Temperature of the gas storage in K.

        Returns:
            Array of the volume flow rates in m^3/d.
        """
        assert (
            np.min(load) >= 0 and np.max(load) <= 1.0
        ), "Load values must be between 0 and 1."

        mass_flow_rate = self.max_power * load / (self.lhv * self.thermal_efficiency)
        vol_flow_rate = mass_flow_rate * self.R_S * temp / (press * 1.0e5)

        return vol_flow_rate * 60 * 60 * 24


@dataclass
class Disturbances:
    """
    Class used to store the desired setup of one or multiple disturbances.

    Attributes:
        state_jumps:
            Dictionary describing the jump of individual states.
            Its keys describe the respective state indices (counting from 0).
            Its values are tuples consisting of the time index (counting from 0) at which
            the jump shall occure as well as a fraction of its nominal value (used for
            normalization) by which it shall jump.
        max_feeding_error:
            Tuple listing the maximum fraction by which each substrates' (in order) actual feeding
            may differ from the computed feeding value.
        feed_computation_stuck:
            Tuple describing the times at which the biogas plant is being fed with the substrates
            last computed before the disturbance occured.
            It consists of the starting time index (counting from 0) as well as the
            number of time steps at which the feeding is "stuck".
        clogged_feeding:
            Dictionary describing the times at which the biogas plant is not fed at all with certain
            substrates.
            Its key indicate the respective substrate index.
            Its values consist of tuples that contain the start time index and the number
            of time steps for which this particular substrate feeder is clogged.

    Example of a 'Disturbances' object initialization:
        disturbances = Disturbances(
            state_jumps={0: (5, 0.1), 1: (4, 0.2)},
            max_feeding_error=(0.1, 0.3, 0.2),
            feed_computation_stuck=(7, 2),
            clogged_feeding={1: (10, 2)},
        )
    """

    state_jumps: dict[int, tuple[int, float]] = None
    max_feeding_error: tuple[float] = None
    feed_computation_stuck: tuple[int, int] = None
    clogged_feeding: dict[int, tuple[int, int]] = None