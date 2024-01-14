from dataclasses import dataclass, field
import numpy as np
from enum import Enum, auto
from ad_meal_prep_control import params_R3
from dataclasses_json import dataclass_json
from copy import deepcopy


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
        ch4_outflow_rate = mass_flow_rate * self.R_S * temp / (press * 1.0e5)

        return ch4_outflow_rate * 60 * 60 * 24


@dataclass(order=True)
class Time:
    hour: int
    minute: int = 0


@dataclass
class TimeSlot:
    start_time: Time
    end_time: Time


def typical_ch4_vol_flow_rate(max_power: float, n_steps: int, t_step: float):
    schedule = {
        0: (TimeSlot(Time(7), Time(15)), TimeSlot(Time(16), Time(22))),
        1: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(15), Time(22))),
        2: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(16), Time(22))),
        3: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(15), Time(22))),
        4: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(16), Time(23))),
        5: (TimeSlot(Time(9), Time(12)), TimeSlot(Time(17), Time(23))),
        6: (
            TimeSlot(Time(0), Time(1)),
            TimeSlot(Time(11), Time(12)),
            TimeSlot(Time(17), Time(23, 59)),
        ),
        # 0: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 1: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 2: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 3: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 4: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 5: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
        # 6: (TimeSlot(Time(6), Time(12)), TimeSlot(Time(15), Time(21))),
    }

    # Set up CHP
    chp = CHP(max_power=max_power)
    chp_load = np.zeros(n_steps)
    for i in range(n_steps):
        day_floating_point = (i * t_step) % 7
        total_minutes = round(day_floating_point * 24 * 60)

        day = total_minutes // (24 * 60)
        hour = (total_minutes % (24 * 60)) // 60
        minute = total_minutes % 60

        time_slots = schedule[day]
        for slot in time_slots:
            if (
                Time(hour, minute) >= slot.start_time
                and Time(hour, minute) < slot.end_time
            ):
                chp_load[i] = 1.0
                break

    return chp.ch4_vol_flow_rate(
        load=chp_load, press=params_R3.p_gas_storage, temp=params_R3.T_gas_storage
    )


@dataclass_json
@dataclass(kw_only=True)
class Disturbances:
    """
    Class used to store the desired setup of one or multiple disturbances.

    Attributes:
        state_jumps:
                                Dictionary describing the jump of individual states.
                                Its keys describe the respective state indices (counting from 0).
                                Its values are lists of tuples consisting of the time index (counting from 0) at which
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
        dictated_feeding:
                                Dictionary describing the substrates whose feed is dictated.
                                Keys: Strings that correspond to substrate names.
                                Values: Tuple consisting of start time, end time and normalized feed volume flow.

    Example of a 'Disturbances' object initialization:
        disturbances = Disturbances(
            state_jumps={0: [(5, 0.1)], 1: [(4, 0.2)]},
            max_feeding_error=(0.1, 0.3, 0.2),
            feed_computation_stuck=(7, 2),
            clogged_feeding={1: (10, 2)},
            dictated_feeding={"CATTLE_MANURE": (0.5, 1.2, 0.3)},
        )
    """

    state_jumps: dict[int, list[tuple[int, float]]] = None
    max_feeding_error: tuple[float] = None
    feed_computation_stuck: tuple[int, int] = None
    clogged_feeding: dict[int, tuple[int, int]] = None
    dictated_feeding: dict[str, tuple[float, float, float]] = None


@dataclass(kw_only=True)
class CostFunction:
    """
    Stores the l-term and m-term of the cost function.

    Attributes:
        lterm:
                Stage cost
        mterm:
                Terminal cost
    """

    lterm: str
    mterm: str


@dataclass(kw_only=True)
class Bound:
    """
    Stores the parameters required to set up a hard bound.
    Refer to https://www.do-mpc.com/en/latest/api/do_mpc.controller.MPC.html#bounds for more details abound do-mpc bounds.

    Attributes:
        lower:
                    If true, this is a lower bound. If false, this is an upper bound.
        variable:
                    The variable it refers to.
        value:
                    The value that constitues the lower/upper bound of the specified variable.
    """

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
    """
    Stores the parameters required to set up a nonlinear constraint.
    Refer to https://www.do-mpc.com/en/latest/api/do_mpc.controller.MPC.html#set-nl-cons for more details about do-mpc nl_cons.

    Attributes:
        expression:
                            The nonlinear constraint expression.
        ub:
                            The upper bound of the specified nonlinear expression.
        soft_constraint:
                            Whether or not the constraint shall be enforced as a soft or hard constraint.
                            If true, the soft constraint is implemented using slack variables within do-mpc.
        penalty_term_cons:
                            Penalty term constant
        maximum_violation:
                            Maximum violation that is allowed.
    """

    expression: str
    ub: float
    soft_constraint: bool
    penalty_term_cons: float
    maximum_violation: float = np.inf


class StateObserver(Enum):
    """
    Enum used to setup different scenarios.
    """

    MHE = auto()
    STATEFEEDBACK = auto()


@dataclass
class LimitedSubstrate:
    """
    Stores the parameters required to set up a limited substrate scenario, i.e. a scenario
    in which a certain amount of a substrate shall be fed within a specified amount of time.
    After that time, there is no more mass of that substrate available for feeding.

    Attributes:
        name:
                            Name of the limited substrate.
        amount_remaining:
                            The amount of the substrate remaining in kilograms.
        days_remaining:
                            The amount of days remaining to feed that substrate.
    """

    name: str
    amount_remaining: float
    days_remaining: float


@dataclass
class SetpointFunction:
    """
    This class holds information required to follow setpoints.

    Attributes:
        setpoints:
                        Setpoint values. 1D-Numpy array with one more element than time_points.
        time_points:
                        Time points [d] at which the setpoints change.
                        First time point specifies when the change happens from the 0th setpoint to the 1st.
    """

    setpoints: np.ndarray
    time_points: np.ndarray

    def __post_init__(self):
        assert np.all(
            np.diff(self.time_points) > 0
        ), f"Supplied time points {self.time_points} are not in strictly increasing order"
        assert (
            len(self.setpoints) == len(self.time_points) + 1
        ), "You must supply one more setpoint than the number of time points."

    def get_current_setpoint(self, time: float):
        current_time_idx = np.argmax(self.time_points > time)
        if np.all(self.time_points < time):
            current_time_idx = -1
        return self.setpoints[current_time_idx]


@dataclass(kw_only=True)
class ControllerParams:
    mpc_n_horizon: int
    mpc_n_robust: int
    num_std_devs: float
    cost_func: CostFunction
    gas_storage_bound_fraction: float = 0.05
    substrate_cost_formulation: str = "quadratic"
    bounds: list[Bound] = None
    nl_cons: list[NlConstraint] = None
    rterm: str = None
    ch4_set_point_function: SetpointFunction = None


@dataclass_json
@dataclass(kw_only=True)
class Scenario:
    name: str
    external_gas_storage_model: bool
    t_step: float  # Time in days
    n_days_steady_state: float
    n_days_mpc: float
    sub_names: list[str]
    disturbances: Disturbances
    x0_true: np.ndarray
    Tx: np.ndarray
    Ty: np.ndarray
    u_max: dict[str, float]
    plot_vars: list[str]
    state_observer: StateObserver
    mhe_n_horizon: int = 5
    controller_params: ControllerParams
    num_std_devs_sim: float
    simulate_steady_state: bool = True
    simulate_mpc: bool = True
    mpc_live_vis: bool = False
    pygame_vis: bool = False
    save_results: bool = True
    compile_nlp: bool = False
    P_el_chp: float = None
    limited_substrates: list[LimitedSubstrate] = None
    _state_names: list[str] = field(
        default_factory=lambda: [
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
    )

    _meas_names: list[str] = field(
        default_factory=lambda: [
            "V´_g",
            "p_CH4",
            "p_CO2",
            "pH",
            "S_IN",
            "TS",
            "VS",
            "S_ac",
        ]
    )

    Tu: np.ndarray = None


class ScenarioFactory:
    _x0_true_default = np.array(
        [
            0.0494,
            0.0116,
            4.97,
            0.963,
            957.0,
            1.48,
            1.48,
            0.948,
            0.412,
            1.92,
            0.521,
            1,
            0.0487,
            0.0493,
            4.54,
            0.0223,
            0.358,
            0.660,
            0.20 * params_R3.V_GAS_STORAGE_MAX,  # m^3
            0.20 * params_R3.V_GAS_STORAGE_MAX,  # m^3
        ]
    )

    _Tx_default = np.array(
        [
            0.137,
            0.0127,
            4.79,
            0.950,
            958.0,
            2.59,
            8.09,
            1.46,
            0.624,
            1.45,
            0.421,
            14.0,
            0.0487,
            0.137,
            4.42,
            0.0297,
            0.380,
            0.569,
            params_R3.V_GAS_STORAGE_MAX,
            params_R3.V_GAS_STORAGE_MAX,
        ]
    )

    _u_max_default = {"solid": 80.0, "liquid": 450.0}

    _Ty_default = np.array(
        [
            450.0,
            0.574,
            0.376,
            7.4,
            0.850,
            0.0422,
            0.668,
            0.0959,
        ]
    )

    _default_dict = {
        "t_step": 0.5 / 24,
        "n_days_steady_state": 30,
        "n_days_mpc": 30,
        "sub_names": [
            "CORN_SILAGE",
            "GRASS_SILAGE",
            "CATTLE_MANURE",
            "SUGAR_BEET_SILAGE",
        ],
        "disturbances": Disturbances(),
        "x0_true": _x0_true_default,
        "Tx": _Tx_default,
        "Ty": _Ty_default,
        "u_max": _u_max_default,
        "plot_vars": [],
        "state_observer": StateObserver.STATEFEEDBACK,
        "mhe_n_horizon": 5,
        "num_std_devs_sim": 2.0,
    }

    methanation_dict = deepcopy(_default_dict)
    methanation_dict["external_gas_storage_model"] = False

    cogeneration_dict = deepcopy(_default_dict)
    cogeneration_dict["external_gas_storage_model"] = True

    def create_scenario(
        self, scenario_type: str, controller_params: ControllerParams, **kwargs
    ) -> Scenario:
        # Create default dict based on scenario_type
        if scenario_type == "methanation":
            scenario_dict = self.methanation_dict
        elif scenario_type == "cogeneration":
            scenario_dict = self.cogeneration_dict
        else:
            raise NotImplementedError(f"Invalid scenario type {scenario_type}")

        # Edit kwargs with more specific options
        scenario_dict["controller_params"] = controller_params
        for key, value in kwargs.items():
            scenario_dict[key] = value
        return Scenario(**scenario_dict)
