import simulation
from utils import ScenarioFactory, CostFunction, ControllerParams
from params_R3 import P_el_chp

lterm = "100*(model.aux['y_1_norm'] - 1.)**2"
mterm = "1000*(model.aux['y_1_norm'] - 1.)**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

controller_params = ControllerParams(
    mpc_n_horizon=10,
    mpc_n_robust=0,
    num_std_devs=1.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
)

kwargs = { "pygame_vis": True, "mpc_live_vis": True, "plot_vars":["u_norm", "y_meas_1", "y_meas_4",],"P_el_chp": P_el_chp}

scenario = ScenarioFactory().create_scenario("cogeneration", controller_params=controller_params, **kwargs)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
