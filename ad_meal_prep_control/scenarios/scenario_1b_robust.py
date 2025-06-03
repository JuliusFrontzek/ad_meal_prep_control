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

# user input:
n_days_mpc = 30         # length of simulation [d]
n_std_dev = 1.5         # used for plant and controller
t_stp_ahead_pred = 3    # for controller plotting
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
    mpc_n_robust=1,
    num_std_devs=n_std_dev,
    cost_func=cost_func,
    substrate_cost_formulation="quadratic",
    ch4_set_point_function=ch4_set_point_function,
    rterm=rterm,
)


kwargs = {
    "pygame_vis": False,
    "mpc_live_vis": False,
    "disturbances": Disturbances(
        dictated_feeding={
            "CATTLE_MANURE_VERY_UNCERTAIN": [
                (5.0, 10.0, 1.292e-3),
                (13.0, 17.0, 2.585e-3),
                (22.0, 26.0, 3.877e-3),
            ],
        },
        #max_feeding_error=0.05,
    ),
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": n_std_dev,
    "feedback": True,
    "mismatch": True,
    "t_stp_ahead_pred": t_stp_ahead_pred,
    "name": f"Scenario_1b_quadratic_robust_feedback_mismatch_{n_std_dev}std_{t_stp_ahead_pred}tsap",
}

scenario = ScenarioFactory().create_scenario(
    "methanation", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
