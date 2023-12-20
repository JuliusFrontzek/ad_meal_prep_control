from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np

lterm = "30*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_out_mean'])/model.tvp['v_ch4_dot_tank_out_mean'])**2"
mterm = "300*((model.aux['v_ch4_dot_tank_in'] - model.tvp['v_ch4_dot_tank_out_mean'])/model.tvp['v_ch4_dot_tank_out_mean'])**2 + 1*(model.aux['y_1_norm'] - 1.)**2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

n_days_mpc = 30

controller_params = ControllerParams(
    mpc_n_horizon=20,
    mpc_n_robust=0,
    num_std_devs=2.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
)

kwargs = {
    "name": "Scenario_2a",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": P_el_chp,
    "disturbances": Disturbances(
        state_jumps={18: [(5, 0.1)], 19: [(4, 0.2)]},
        dictated_feeding={
            "CATTLE_MANURE_VERY_UNCERTAIN": (5.0, 10.0, 0.1),
        },
    ),
    "n_days_mpc": n_days_mpc,
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
