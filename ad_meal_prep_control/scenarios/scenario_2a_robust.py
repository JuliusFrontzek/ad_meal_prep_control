from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
    NlConstraint
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np

lterm = ("(0.5*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - 0.4)**2 +"
         "+ 50*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - 0.4)**4)"
         )

mterm = "model.tvp['dummy_tvp']"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

# user input:
n_days_mpc = 30         # length of simulation [d]
n_std_dev = 3           # effectively only used for plant
t_stp_ahead_pred = 8    # for controller plotting

controller_params = ControllerParams(
    mpc_n_horizon=40,
    mpc_n_robust=1,
    num_std_devs=n_std_dev,
    cost_func=cost_func,
    substrate_cost_formulation="linear",
    gas_storage_bound_fraction=0.05,
    #nl_cons=[NlConstraint(expression='7-model.aux["y_4"]', ub=0, soft_constraint=True, penalty_term_cons=10e0)]
)

kwargs = {
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
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": n_std_dev,
    "feedback": True,
    "mismatch": True,
    "t_stp_ahead_pred": t_stp_ahead_pred,
    "name": f"Scenario_2a_dynamic_robust_feedback_mismatch_{n_std_dev}std_{t_stp_ahead_pred}tsap",
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
