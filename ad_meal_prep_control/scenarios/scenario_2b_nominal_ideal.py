from ad_meal_prep_control import simulation
from ad_meal_prep_control.utils import (
    ScenarioFactory,
    CostFunction,
    ControllerParams,
    Disturbances,
)
from ad_meal_prep_control.params_R3 import P_el_chp
import numpy as np

fill_level_setpoint = 0.4
lterm = (f"(5e2*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - {fill_level_setpoint})**2 +"
         f"+ 5e1*(model.x['x_19'] + model.x['x_20'] + model.aux['V_H2O']/V_GAS_STORAGE_MAX - {fill_level_setpoint})**4)"
         )

mterm = "model.tvp['dummy_tvp']"
cost_func = CostFunction(lterm=lterm, mterm=mterm)

# user input:
n_days_mpc = 2.5  # length of simulation [d]
n_std_dev = 0  # number std deviations

controller_params = ControllerParams(
    mpc_n_horizon=48,
    mpc_n_robust=0,
    num_std_devs=n_std_dev,
    cost_func=cost_func,
    substrate_cost_formulation="linear",
    gas_storage_bound_fraction=0.05,
)

kwargs = {
    "name": "Scenario_2b_nominal_ideal",
    "pygame_vis": False,
    "mpc_live_vis": False,
    "P_el_chp": 50.0,
    "t_step": 0.5 / 24.0,
    "plot_vars": [
            "u_norm",
            "y_meas_1",
            "v_ch4_dot_tank_in",
            "y_meas_4",
    ],
    "disturbances": Disturbances(
        dictated_feeding={
            "MAIZE_SILAGE": [
                            (1.0, 2.0, 0.02),
                            ],
        },
    ),
    "n_days_mpc": n_days_mpc,
    "num_std_devs_sim": 0,
}

scenario = ScenarioFactory().create_scenario(
    "cogeneration", controller_params=controller_params, **kwargs
)

sim = simulation.Simulation(scenario=scenario)
sim.setup()
sim.run()
