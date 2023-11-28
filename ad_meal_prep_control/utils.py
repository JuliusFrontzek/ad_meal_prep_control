from dataclasses import dataclass
from typing import Union
import numpy as np
from casadi import SX
from enum import Enum, auto
import utils


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


@dataclass(kw_only=True)
class CostFunction:
    mterm: str
    lterm: str


@dataclass(kw_only=True)
class Bound:
    lower: bool
    variable: str
    value: float

    def __post_init__(self):
        assert (
            self.variable[0] in ["u", "x"] and self.variable[1] == "_"
        ), f"Variable {self.variable} currently not supported for bounds."
        self.variable_type = f"_{self.variable[0]}"

    @property
    def direction(self):
        if self.lower:
            return "lower"
        else:
            return "upper"


@dataclass(kw_only=True)
class NlConstraint:
    expression: str
    ub: float
    soft_constraint: bool
    penalty_term_cons: float
    maximum_violation: float = np.inf


class StateObserver(Enum):
    MHE = auto()
    STATEFEEDBACK = auto()


@dataclass(kw_only=True)
class ScenarioData:
    name: str
    external_gas_storage_model: bool
    mpc_n_horizon: int
    mpc_n_robust: int
    t_step: float  # Time in days
    n_days_steady_state: float
    n_days_mpc: float
    sub_names: list[str]
    disturbances: utils.Disturbances
    x0_true: np.ndarray
    Tx: np.ndarray
    Ty: np.ndarray
    u_max: dict[str, float]
    num_std_devs: float  # Lower and upper bound of uncertainties is determined by the number of standard deviations that we consider
    plot_vars: list[str]
    state_observer: StateObserver
    mhe_n_horizon: int = 5
    cost_func: CostFunction
    consider_substrate_costs: bool = True
    bounds: list[Bound] = None
    nl_cons: list[NlConstraint] = None
    rterm: str = None
    consider_uncertainty: bool = True
    simulate_steady_state: bool = True
    simulate_mpc: bool = True
    mpc_live_vis: bool = True
    pygame_vis: bool = False
    store_results: bool = True
    compile_nlp: bool = False
    vol_flow_rate: Union[np.ndarray, None] = None

    _state_names = [
        "S_ac",
        "S_ch4",
        "S_IC",
        "S_IN",
        "S_h2o",
        "X_ch_f",
        "X_ch_s",
        "X_pr",
        "X_li",
        "X_bac",
        "X_ac",
        "X_ash",
        "S_ion",
        "S_ac−",
        "S_hco3−",
        "S_nh3",
        "S_ch4_gas",
        "S_co2_gas",
        "V_CH4",
        "V_CO2",
    ]

    _meas_names = ["V´_g", "p_CH4", "p_CO2", "pH", "S_IN", "TS", "VS", "S_ac"]
