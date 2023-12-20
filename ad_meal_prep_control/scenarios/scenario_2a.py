import simulation
from utils import ScenarioFactory, CostFunction, ControllerParams, Disturbances
from params_R3 import P_el_chp

# lterm = mterm = "model.tvp['dummy_tvp']"
lterm = "30*((model.aux['v_ch4_dot_tank_in'] - 500.)/500.)**2"
mterm = "300*((model.aux['v_ch4_dot_tank_in'] - 500.)/500.)**2 + 1*(model.aux['y_1_norm'] - 1.)**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

controller_params = ControllerParams(
    mpc_n_horizon=10,
    mpc_n_robust=0,
    num_std_devs=1.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
)

kwargs = {
    "name": "cogeneration_07_12",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": P_el_chp,
    "plot_vars": [
        "u_norm",
        "y_meas_1",
        "v_ch4_dot_tank_in",
        "dictated_sub_feed_1",
        "y_meas_4",
    ],
    "disturbances": Disturbances(
        dictated_feeding={
            "GRASS_SILAGE": (0.2, 0.4, 0.1),
            "CATTLE_MANURE": (3.0, 5.0, 0.05),
        }
    ),
    "n_days_mpc": 10,
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
