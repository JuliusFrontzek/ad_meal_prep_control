from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    SetpointFunction,
    Disturbances,
)
import numpy as np
import math

c_1 = 1e3
c_2 = 0
c_3 = 1e-1
lterm = f"{c_1}*((model.aux['v_dot_ch4_AD_norm_condition'] - model.tvp['v_ch4_dot_AD_setpoint'])/model.tvp['v_ch4_dot_AD_setpoint'])**2"
mterm = f"{c_2}*((model.aux['v_dot_ch4_AD_norm_condition'] - model.tvp['v_ch4_dot_AD_setpoint'])/model.tvp['v_ch4_dot_AD_setpoint'])**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

# user input:
n_days_mpc = 28         # length of simulation [d]
n_std_dev = 1           # used for plant and controller
t_stp_ahead_pred = 3    # for controller plotting
substrate_cost_formulation = "quadratic"  # linear or quadratic

ch4_set_point_function = SetpointFunction(
    setpoints=np.array([350.0, 550.0, 450.0, 350.0]),
    time_points=np.array([3, 6, 9]),
)

rterms = [
    f"{c_3}*(model.u['u_norm'][{i}] - mpc.u_prev['u_norm'][{i}])**2" for i in range(4)
]
rterm = " + ".join(rterms)

# No feedback
controller_params = ControllerParams(
    mpc_n_horizon=15,
    mpc_n_robust=0,
    num_std_devs=n_std_dev,
    cost_func=cost_func,
    substrate_cost_formulation="quadratic",
    ch4_set_point_function=ch4_set_point_function,
    rterm=rterm,
)

kwargs = {
    "name": "Scenario_1a_quadratic_no_feedback_mismatch_1std_ch",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": n_std_dev,
    "feedback": False,
    "mismatch": True,
    "t_stp_ahead_pred": t_stp_ahead_pred
}

scenario = ScenarioFactory().create_scenario(
    "methanation", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
