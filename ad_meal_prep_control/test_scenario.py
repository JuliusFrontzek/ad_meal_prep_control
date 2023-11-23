import scenario
from utils import Disturbances
import numpy as np
import params_R3
from utils import ScenarioType, StateObserver, CHP, ScenarioData, CostFunction

disturbances = Disturbances()

x0 = np.array(
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
        39.0 * 1.7 / params_R3.SCALEDOWN,  # m^3
        36.0 * 1.7 / params_R3.SCALEDOWN,  # m^3
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
    "solid": 80_000.0 / params_R3.SCALEDOWN,
    "liquid": 450_000.0 / params_R3.SCALEDOWN,
}

Ty = np.array(
    [
        250.0,  # 140.300279906936
        0.574083930894918,
        0.376314347120225,
        5.5,
        0.850445630702126,
        0.0422284958830547,
        0.668470313534998,
        0.0959467827166042,
    ]
)

n_days_steady_state = 0.5
n_days_mpc = 3
t_step = 0.5 / 24

mpc_n_horizon = 5
mpc_n_robust = 0
mhe_n_horizon = 10

total_steps = round((n_days_steady_state + n_days_mpc) / t_step)

n_steps_steady_state = round(n_days_steady_state / t_step)

# Set up CHP
chp = CHP(max_power=124.0 / params_R3.SCALEDOWN)
chp_load = np.zeros(total_steps + mpc_n_horizon)
for i in range(6):
    chp_load[n_steps_steady_state + i :: 48] = 1.0  # 6:00 - 12:00
    chp_load[n_steps_steady_state + 18 + i :: 48] = 1.0  # 15:00 - 21:00

vol_flow_rate = chp.ch4_vol_flow_rate(
    load=chp_load, press=params_R3.p_gas_storage, temp=params_R3.T_gas_storage
)

mterm = (
    lterm
) = "(model.aux['y_1_norm'] - 1.) ** 2"  # (self.model.aux["y_1_norm"] - 1.0) ** 2

cost_func = CostFunction(mterm=mterm, lterm=lterm)

test_scenario_data = ScenarioData(
    name="test_scenario",
    scenario_type=ScenarioType.METHANATION,
    mpc_n_horizon=mpc_n_horizon,
    mpc_n_robust=mpc_n_robust,
    t_step=t_step,
    n_days_steady_state=n_days_steady_state,
    n_days_mpc=n_days_mpc,
    sub_names=["CORN_SILAGE", "GRASS_SILAGE", "CATTLE_MANURE"],
    disturbances=disturbances,
    x0=x0,
    Tx=Tx,
    Ty=Ty,
    u_max=u_max,
    num_std_devs=0.5,
    plot_vars=[
        "u_norm",
        "y_meas_1",
        "y_meas_4",
    ]
    + [f"x_{i+1}" for i in range(18)],
    state_observer=StateObserver.MHE,
    mhe_n_horizon=mhe_n_horizon,
    cost_func=cost_func,
    consider_uncertainty=True,
    simulate_steady_state=False,
    simulate_mpc=True,
    mpc_live_vis=True,
    pygame_vis=True,
    store_results=True,
    compile_nlp=False,
    vol_flow_rate=vol_flow_rate,
)

test_scenario = scenario.Scenario(scenario_data=test_scenario_data)

test_scenario.setup()
test_scenario.run()
