import simulation
from utils import ScenarioFactory, CostFunction, ControllerParams, SetpointFunction
import numpy as np

lterm = "100*((model.aux['y_1'] - model.tvp['v_ch4_dot_out_setpoint'])/model.tvp['v_ch4_dot_out_setpoint'])**2"
mterm = "1000*((model.aux['y_1'] - model.tvp['v_ch4_dot_out_setpoint'])/model.tvp['v_ch4_dot_out_setpoint'])**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

ch4_set_point_function = SetpointFunction(setpoints=np.array([450., 550., 650., 450.]), time_points=np.array([0.5, 1., 1.5]))

controller_params = ControllerParams(
    mpc_n_horizon=10,
    mpc_n_robust=0,
    num_std_devs=1.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
    ch4_set_point_function=ch4_set_point_function,
)

kwargs = { "pygame_vis": True, "mpc_live_vis": True, "plot_vars":["u_norm", "y_meas_1", "y_meas_4",],}

scenario = ScenarioFactory().create_scenario("methanation", controller_params=controller_params, **kwargs)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
