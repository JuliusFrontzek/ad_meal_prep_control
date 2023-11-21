import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import sys
import os

sys.path.append(os.getcwd())
import ad_meal_prep_control.utils as utils
import params_R3
import copy
import pygame
import ad_meal_prep_control.visualization as vis
import state_estimator
import substrates
import time
from uncertainties import ufloat

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc
from do_mpc.tools import Timer

import matplotlib.pyplot as plt
from models import adm1_r3_frac_norm
from mpc import mpc_setup
from simulator import simulator_setup

""" User settings: """
show_animation = True
store_results = True

# pygame setup
pygame.init()
screen_size = (1280, 720)
screen = pygame.display.set_mode(screen_size)
clock = pygame.time.Clock()

bga = vis.BioGasPlantVis(150.0, screen)
data = vis.DataVis(screen)

# Model
model_type = "continuous"  # either 'discrete' or 'continuous'
model = do_mpc.model.Model(model_type)

compile_nlp = False

subs = [substrates.STANDARD_SUBSTRATE]
xi = [sub.xi for sub in subs]
uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in subs]

uncertain_xis[0] = (ufloat(xi[0][5], 0.0), ufloat(xi[0][7], 0.0), ufloat(xi[0][8], 0.0))

# Simulation
n_days_steady_state = 300
n_days_mpc = 7

t_step = 0.5 / 24  # Time in days
n_steps_steady_state = round(n_days_steady_state / t_step)
n_steps = round(n_days_mpc / t_step)


# MPC
n_horizon = 5

# Set up CHP
chp = utils.CHP(max_power=124.0 / params_R3.SCALEDOWN)
chp_load = np.zeros(n_steps + n_horizon)
for i in range(6):
    chp_load[i::48] = 1.0  # 6:00 - 12:00
    chp_load[18 + i :: 48] = 1.0  # 15:00 - 21:00

vol_flow_rate = chp.ch4_vol_flow_rate(
    load=chp_load, press=params_R3.p_gas_storage, temp=params_R3.T_gas_storage
)

# Set up disturbances
disturbances = utils.Disturbances(
    # state_jumps={0: (12, 0.5), 1: (13, 0.2)},
    # max_feeding_error=(0.01, 0.01, 0.01),
    # feed_computation_stuck=(2, 5),
    # clogged_feeding={1: (2, 10), 1: (2, 10)},
)

# Set the initial state of mpc and simulator
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

state_names = [
    "S_ac",
    "S_ch4",
    "S_IC",
    "S_IN",
    "S_h2o",
    "X_ch_f",
    "X_ch_s",
    "X_pr",
    "X_li",
    "X_bac",
    "X_ac",
    "X_ash",
    "S_ion",
    "S_ac−",
    "S_hco3−",
    "S_nh3",
    "S_ch4_gas",
    "S_co2_gas",
    "V_CH4",
    "V_CO2",
]

meas_names = ["V´_g", "p_CH4", "p_CO2", "pH", "S_IN", "TS", "VS", "S_ac"]

# Normalisation
# Define numeric values for normalization with steady state
# Tx = x0.copy()
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

Tu = np.array([u_max[sub.state] for sub in subs])

Ty = np.array(
    [
        140.300279906936,
        0.574083930894918,
        0.376314347120225,
        7.31094007321728,
        0.850445630702126,
        0.0422284958830547,
        0.668470313534998,
        0.0959467827166042,
    ]
)

xi_norm = list(np.array([val / Tx[:-2] for val in xi]).T)

x0_norm = np.copy(x0)
x0_norm /= Tx

# Set model
model = adm1_r3_frac_norm(xi_norm, Tu, Tx, Ty)

num_std_devs = 0.0  # Lower and upper bound of uncertainties is determined by the number of standard deviations that we consider

x_ch_nom = np.array([un_xi[0].nominal_value for un_xi in uncertain_xis])
x_pr_nom = np.array([un_xi[1].nominal_value for un_xi in uncertain_xis])
x_li_nom = np.array([un_xi[2].nominal_value for un_xi in uncertain_xis])

x_ch_std_dev = np.array([un_xi[0].std_dev for un_xi in uncertain_xis])
x_pr_std_dev = np.array([un_xi[1].std_dev for un_xi in uncertain_xis])
x_li_std_dev = np.array([un_xi[2].std_dev for un_xi in uncertain_xis])

x_ch_in = np.array(
    [
        x_ch_nom - num_std_devs * x_ch_std_dev,
        x_ch_nom,
        x_ch_nom + num_std_devs * x_ch_std_dev,
    ]
)
x_pr_in = np.array(
    [
        x_pr_nom - num_std_devs * x_pr_std_dev,
        x_pr_nom,
        x_pr_nom + num_std_devs * x_pr_std_dev,
    ]
)
x_li_in = np.array(
    [
        x_li_nom - num_std_devs * x_li_std_dev,
        x_li_nom,
        x_li_nom + num_std_devs * x_li_std_dev,
    ]
)

# Normalize uncertain xi's
x_ch_in /= Tx[5]
x_pr_in /= Tx[7]
x_li_in /= Tx[8]

mpc = mpc_setup(
    model=model,
    t_step=t_step,
    n_horizon=n_horizon,
    n_robust=0,
    x_ch_in=x_ch_in,
    x_pr_in=x_pr_in,
    x_li_in=x_li_in,
    compile_nlp=compile_nlp,
    vol_flow_rate=vol_flow_rate,
)

simulator = simulator_setup(
    model=model,
    t_step=t_step,
    x_ch_in=x_ch_in,
    x_pr_in=x_pr_in,
    x_li_in=x_li_in,
    vol_flow_rate=np.zeros(n_steps_steady_state),
)

estimator = state_estimator.StateEstimator(model)


# Feeding
constant_feeding = False

# Set x0 and u0
mpc.x0 = x0_norm
mpc.u0 = np.ones(len(subs)) * 0.5
simulator.x0 = x0_norm

plot_vars = [
    "u_norm",
    "x_13",
    "x_14",
    "x_8",
    "x_19",
    "x_20",
] + [f"y_{i+1}" for i in range(8)]
# plot_vars = [f"y_{i+1}" for i in range(8)]  # + ["x_19", "x_20"]
# plot_vars = ["u"] + [f"x_{i+1}" for i in range(18)]
# plot_vars = ["u", "x_13", "x_14", "x_8"]

if not constant_feeding:
    plot_vars = [var for var in plot_vars if not var.startswith("y")]

if constant_feeding:
    # Create arrays to save results
    results = {}
    for var in plot_vars:
        results[var] = np.empty(n_steps)
else:
    mpc.set_initial_guess()

    # Initialize graphic:
    graphics = do_mpc.graphics.Graphics(mpc.data)

    # Configure plot:
    fig, ax = plt.subplots(len(plot_vars), sharex=True)
    for idx, var in enumerate(plot_vars):
        graphics.add_line(var_type=f"_{var[0]}", var_name=var, axis=ax[idx])
        ax[idx].set_ylabel(var)

    # Update properties for all prediction lines:
    for line_i in graphics.pred_lines.full:
        line_i.set_linewidth(1)

    fig.align_ylabels()
    fig.tight_layout()
    plt.ion()

timer = Timer()


u_norm_computed = None

# Simulate until steady state
for k in range(n_steps_steady_state):
    screen.fill("white")

    u_norm_steady_state = np.array(
        [[8.0 / Tu[0] for _ in range(len(subs))]]
    ).T  # 1.0 / 1e4

    y_next = simulator.make_step(u_norm_steady_state)
    x0_norm = estimator.estimate_x(y_next)

    # simulator._x0.master[-2] = 50.0 / params_R3.SCALEDOWN
    # simulator._x0.master[-1] = 50.0 / params_R3.SCALEDOWN

np.savetxt("results.csv", simulator.data._x, delimiter=",")

x0_ss = simulator._x0.master

simulator = simulator_setup(
    model=model,
    t_step=t_step,
    x_ch_in=x_ch_in,
    x_pr_in=x_pr_in,
    x_li_in=x_li_in,
    vol_flow_rate=vol_flow_rate,
)

simulator.x0 = x0_ss

simulator._x0.master[-2] = 90.0 / params_R3.SCALEDOWN
simulator._x0.master[-1] = 40.0 / params_R3.SCALEDOWN


# MPC
for k in range(n_steps):
    # fill the screen with a color to wipe away anything from last frame
    screen.fill("white")

    timer.tic()
    u_norm_computed_old = copy.copy(u_norm_computed)
    if constant_feeding:
        u_norm_computed = np.array([[0.011111 for _ in range(len(subs))]]).T
    else:
        u_norm_computed = mpc.make_step(x0_norm)

    u_norm_actual = np.copy(u_norm_computed)

    # Manipulate the actual feed to the biogas plant 'u_norm_actual'
    # based on the set disturbances
    if not disturbances.feed_computation_stuck is None:
        stuck_start_idx = disturbances.feed_computation_stuck[0]
        stuck_end_idx = stuck_start_idx + disturbances.feed_computation_stuck[1]
        if k in range(stuck_start_idx, stuck_end_idx):
            u_norm_actual = copy.copy(u_norm_computed_old)
            u_norm_computed = copy.copy(u_norm_computed_old)

    if not disturbances.clogged_feeding is None:
        for sub_idx, val in disturbances.clogged_feeding.items():
            if k in range(val[0], val[0] + val[1]):
                u_norm_actual[sub_idx, 0] = 0.0

    if not disturbances.max_feeding_error is None:
        u_norm_actual *= (
            1.0
            + np.array(
                [
                    np.random.uniform(low=-1.0, high=1.0, size=u_norm_actual.shape[0])
                    * disturbances.max_feeding_error
                ]
            ).T
        )

    timer.toc()

    y_next = simulator.make_step(u_norm_actual)
    x0_norm = estimator.estimate_x(y_next)

    if not disturbances.state_jumps is None:
        for x_idx, val in disturbances.state_jumps.items():
            if k == val[0]:
                x0_norm[x_idx] += val[1]

    if show_animation:
        if constant_feeding:
            if "u" in plot_vars:
                results["u"][k] = u_norm_actual

            for var in plot_vars:
                if var.startswith("x"):
                    results[var][k] = y_next[int(var.lstrip("x_")) - 1]
                elif var.startswith("y"):
                    idx = int(var.lstrip("y_"))
                    results[var][k] = simulator.data._aux[-1, idx]
        else:
            graphics.plot_results(t_ind=k)
            graphics.plot_predictions(t_ind=k)
            graphics.reset_axes()
        plt.show()
        plt.pause(0.01)

        vis.visualize(
            bga,
            data,
            state_names,
            meas_names,
            Tx,
            Ty,
            x0_norm,
            simulator,
            u_norm_actual,
            params_R3.V_GAS_STORAGE_MAX,
        )

if constant_feeding:
    fig, ax = plt.subplots(len(plot_vars), sharex=True)
    for idx, var in enumerate(plot_vars):
        ax[idx].set_ylabel(var)
        ax[idx].plot(results[var])

    plt.show()

timer.info()
timer.hist()

pygame.quit()

# Store results:
if store_results:
    # do_mpc.data.save_results([mpc, simulator], f"testing_results")
    np.savetxt("results.csv", simulator.data._x, delimiter=",")
