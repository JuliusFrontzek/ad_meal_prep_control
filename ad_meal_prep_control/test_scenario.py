import simulation
from utils import Disturbances
import numpy as np
import params_R3
from utils import StateObserver, Scenario, CostFunction, ControllerParams, typical_ch4_vol_flow_rate

disturbances = Disturbances()

x0_true = np.array(
    [
        0.0494667574155131,
        0.0116512808544296,
        4.97521803226548,
        0.963856429890969,
        957.102301745169,
        1.48089980608238,
        1.48089980608238,
        0.948575266222027,
        0.412040527872582,
        1.92558575222279,
        0.521526149938689,
        1,
        0.0487500000000000,
        0.0493342637010683,
        4.54552248672517,
        0.0223970128000483,
        0.358267793064052,
        0.660494133806800,
        0.44 * params_R3.V_GAS_STORAGE_MAX,  # m^3
        0.4 * params_R3.V_GAS_STORAGE_MAX,  # m^3
    ]
)

Tx = np.array(
    [
        0.137434457537417,
        0.0127484119685527,
        4.79932043710544,
        0.950816195802454,
        958.064331770733,
        2.59654516843011,
        8.09630748329204,
        1.46003537123105,
        0.624174795213073,
        1.45262583426474,
        0.421713306327953,
        14.0000000291803,
        0.0487500000000000,
        0.137097486806281,
        4.42830805698549,
        0.0297771563953578,
        0.380487873826158,
        0.569429468392225,
        params_R3.V_GAS_STORAGE_MAX,
        params_R3.V_GAS_STORAGE_MAX,
    ]
)

u_max = {
    "solid": 80.0,
    "liquid": 450.0,
}

Ty = np.array(
    [
        250.0,
        0.574083930894918,
        0.376314347120225,
        7.0,
        0.850445630702126,
        0.0422284958830547,
        0.668470313534998,
        0.0959467827166042,
    ]
)

n_days_steady_state = 0.5
n_days_mpc = 3
t_step = 0.5 / 24

mpc_n_horizon = 10
mpc_n_robust = 0
mhe_n_horizon = 5

n_steps_mpc = round(n_days_mpc / t_step)

lterm = "100*(model.aux['y_1_norm'] - 1.)**2"  # "100*(model.aux['y_1_norm'] - 1.) ** 2"
mterm = "1000*(model.aux['y_1_norm'] - 1.)**2"
# lterm = "fabs(model.aux['y_1_norm'] - 1.) + (model.aux['y_1_norm'] - 1.) ** 2"

cost_func = CostFunction(lterm=lterm, mterm=mterm)

controller_params = ControllerParams(
    mpc_n_horizon=mpc_n_horizon,
    mpc_n_robust=mpc_n_robust,
    num_std_devs=1.0,
    cost_func=cost_func,
    consider_substrate_costs=True,
)


test_scenario_data = Scenario(
    name="test_scenario",
    external_gas_storage_model=True,
    t_step=t_step,
    n_days_steady_state=n_days_steady_state,
    n_days_mpc=n_days_mpc,
    sub_names=[
        "CORN_SILAGE",
        "GRASS_SILAGE",
        "CATTLE_MANURE",
    ],  # "STANDARD_SUBSTRATE",
    disturbances=disturbances,
    x0_true=x0_true,
    Tx=Tx,
    Ty=Ty,
    u_max=u_max,
    plot_vars=[
        "u_norm",
        "y_meas_1",
        "y_meas_4",
    ],
    # + [f"x_{i+1}" for i in range(18)],
    state_observer=StateObserver.STATEFEEDBACK,
    mhe_n_horizon=mhe_n_horizon,
    controller_params=controller_params,
    num_std_devs_sim=1.0,
    simulate_steady_state=True,
    simulate_mpc=True,
    mpc_live_vis=True,
    pygame_vis=True,
    save_results=True,
    compile_nlp=False,
)

test_scenario = simulation.Simulation(scenario=test_scenario_data)

test_scenario.setup()
test_scenario.run()
