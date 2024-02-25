from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np
from numpy.random import default_rng
import random

np.random.seed(seed=42)

lterm = "(0.5*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - 0.43)**2 + 50.*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - 0.43)**4)"

mterm = "model.tvp['dummy_tvp']"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

n_days_mpc = 30

controller_params = ControllerParams(
    mpc_n_horizon=40,
    mpc_n_robust=1,
    num_std_devs=2.0,
    cost_func=cost_func,
    substrate_cost_formulation="linear",
    gas_storage_bound_fraction=0.05,
)

t_step = 0.5 / 24.0

rng = default_rng()
mpc_t_steps = int(n_days_mpc / t_step)

state_jumps_ch4 = []
state_jumps_co2 = []

for i in range(mpc_t_steps):
    if i % int(5 / t_step / 24) == 0:
        state_jumps_ch4.append((i, random.random() * 0.06 - 0.03))
        state_jumps_co2.append((i, random.random() * 0.06 - 0.03))
    else:
        state_jumps_ch4.append((i, random.random() * 0.02 - 0.01))
        state_jumps_co2.append((i, random.random() * 0.02 - 0.01))


kwargs = {
    "name": "Scenario_2c_dynamic_test",
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
            "CATTLE_MANURE_VERY_UNCERTAIN": [(5.0, 10.0, 0.01)],
        },
        max_feeding_error=0.05,
    ),
    "n_days_mpc": n_days_mpc,
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
