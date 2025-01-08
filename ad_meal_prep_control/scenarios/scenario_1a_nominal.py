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

lterm = "10*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2"
mterm = "100*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

n_days_mpc = 30
setpoints = np.array([350.0, 550.0, 450.0, 350.0])

ch4_set_point_function = SetpointFunction(
    setpoints=setpoints,
    time_points=np.array([3, 6, 9]),
)

rterms = [
    f"0.1*(model.u['u_norm'][{i}] - mpc.u_prev['u_norm'][{i}])**2" for i in range(4)
]
rterm = " + ".join(rterms)

controller_params = ControllerParams(
    mpc_n_horizon=15,
    mpc_n_robust=0,
    num_std_devs=5.0,
    cost_func=cost_func,
    substrate_cost_formulation="quadratic",
    ch4_set_point_function=ch4_set_point_function,
    rterm=rterm,
)


kwargs = {
    "name": "Scenario_1a_quadratic_nominal_feedback_mismatch_5std_3tsap",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": 5.0,
    "feedback": True,
    "mismatch": True,
    "t_stp_ahead_pred": 3
}

scenario = ScenarioFactory().create_scenario(
    "methanation", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
