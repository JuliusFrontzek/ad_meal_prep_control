from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np

lterm = "(0.5*(model.x['x_19'] + model.x['x_20'] - 0.5)**2 + 25.*(model.x['x_19'] + model.x['x_20'] - 0.5)**4)"  # "10*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_out_mean'])/model.tvp['v_ch4_dot_tank_out_mean'])**2"
mterm = "model.tvp['dummy_tvp']"  # "(model.x['x_19'] + model.x['x_20'] - 0.5)**2"  # "100*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_out_mean'])/model.tvp['v_ch4_dot_tank_out_mean'])**2"


cost_func = CostFunction(lterm=lterm, mterm=mterm)

n_days_mpc = 30

rterms = [
    f"0.03*(model.u['u_norm'][{i}] - mpc.u_prev['u_norm'][{i}])**2" for i in range(4)
]
rterm = " + ".join(rterms)

controller_params = ControllerParams(
    mpc_n_horizon=24,
    mpc_n_robust=0,
    num_std_devs=2.0,
    cost_func=cost_func,
    substrate_cost_formulation="quadratic",
    # rterm=rterm,
)

kwargs = {
    "name": "Scenario_2a_test",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": P_el_chp,
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
