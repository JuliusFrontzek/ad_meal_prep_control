import simulation
from utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    SetpointFunction,
    Disturbances,
)
import numpy as np

lterm = "1000*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2 + 1*(model.aux['y_1_norm'] - 1.)**2"
mterm = "1000*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2 + 1*(model.aux['y_1_norm'] - 1.)**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

ch4_set_point_function = SetpointFunction(
    setpoints=np.array([450.0, 550.0, 650.0, 550.0]),
    time_points=np.array([0.5, 1.0, 1.5]),
)

controller_params = ControllerParams(
    mpc_n_horizon=15,
    mpc_n_robust=0,
    num_std_devs=1.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
    ch4_set_point_function=ch4_set_point_function,
)

kwargs = {
    "pygame_vis": True,
    "mpc_live_vis": True,
    "plot_vars": [
        "u_norm",
        "y_meas_1",
        "v_ch4_dot_tank_in",
        # "dictated_sub_feed_1",
        "y_meas_4",
    ],
    "disturbances": Disturbances(
        dictated_feeding={
            "CATTLE_MANURE": (0.2, 0.4, 1.0),
            # "GRASS_SILAGE": (0.1, 0.3, 0.3),
        }
    ),
}

scenario = ScenarioFactory().create_scenario(
    "methanation", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
