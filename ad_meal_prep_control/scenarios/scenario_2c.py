from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
import numpy as np
from numpy.random import default_rng
import random

np.random.seed(seed=42)

# user input:
n_days_mpc = 28
n_std_dev = 1  # number std deviations
t_step = 0.5 / 24.0

fill_level_setpoint = 0.5
c_1 = 1e3
lterm = (f"{c_1} * (model.aux['v_gas_storage']/V_GAS_STORAGE_MAX - {fill_level_setpoint})**2"
         )

mterm = "model.tvp['dummy_tvp']"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

controller_params = ControllerParams(
    mpc_n_horizon=48,
    mpc_n_robust=1,
    num_std_devs=n_std_dev,
    cost_func=cost_func,
    substrate_cost_formulation="linear",
    gas_storage_bound_fraction=0.05,
)

# add gas storage measurement noise ("state jumps"):
rng = default_rng()
mpc_t_steps = int(n_days_mpc / t_step)

state_jumps_ch4 = []
state_jumps_co2 = []

for i in range(mpc_t_steps):
    if i % int(5 / t_step / 24) == 0:
        state_jumps_ch4.append((i, random.random() * 0.06 - 0.03))  # +/- 3%
        state_jumps_co2.append((i, random.random() * 0.06 - 0.03))
    else:
        state_jumps_ch4.append((i, random.random() * 0.02 - 0.01))  # +/- 1%
        state_jumps_co2.append((i, random.random() * 0.02 - 0.01))


kwargs = {
    "name": "Scenario_2c_dynamic",
    # "name": "Scenario_2c_uninhibited",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": 50.0,
    "t_step": t_step,
    "plot_vars": [
        "u_norm",
        "y_meas_1",
        "v_ch4_dot_tank_in",
        "y_meas_4",
    ],
    "disturbances": Disturbances(
        state_jumps={18: state_jumps_ch4, 19: state_jumps_co2},
        dictated_feeding={
            "CATTLE_MANURE_VERY_UNCERTAIN": [(5.0, 7.0, 1.292e-3),
                                             (9.0, 12.0, 2.585e-3),
                                             (15.0, 19.0, 3.877e-3)],
        },
        #max_feeding_error=0.05,
    ),
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": n_std_dev,
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
