from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np

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

kwargs = {
    "name": "Scenario_2a_dynamic",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": 50.0,
    "plot_vars": [
        "u_norm",
        "y_meas_1",
        "v_ch4_dot_tank_in",
        "y_meas_4",
    ],
    "t_step": 0.5 / 24.0,
    "disturbances": Disturbances(
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
