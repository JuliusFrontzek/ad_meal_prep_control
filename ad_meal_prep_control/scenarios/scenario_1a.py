import simulation
from utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    SetpointFunction,
    Disturbances,
)
import numpy as np

lterm = "30*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2"
mterm = "300*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_in_setpoint'])/model.tvp['v_ch4_dot_tank_in_setpoint'])**2 + 1*(model.aux['y_1_norm'] - 1.)**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

n_days_mpc = 30

ch4_set_point_function = SetpointFunction(
    setpoints=np.array(
        [[350.0, 450.0, 550.0] for _ in range(round(n_days_mpc / 3))]
    ).flatten(),
    time_points=np.array([i for i in range(1, n_days_mpc)]),
)

controller_params = ControllerParams(
    mpc_n_horizon=15,
    mpc_n_robust=0,
    num_std_devs=2.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
    ch4_set_point_function=ch4_set_point_function,
)


kwargs = {
    "name": "Scenario_1a",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "n_days_mpc": n_days_mpc,
}

scenario = ScenarioFactory().create_scenario(
    "methanation", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()