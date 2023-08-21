import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import sys
import os
import ad_meal_prep_control.utils as utils
import params_R3
import copy

# import ad_meal_prep_control.state_estimator as state_estimator
import state_estimator
import substrates

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

# Model
# substrate_names = ["corn", "corn"]

model_type = "continuous"  # either 'discrete' or 'continuous'
model = do_mpc.model.Model(model_type)


# _, x_ch_nom, x_pr_nom, x_li_nom = substrates.xi_values(substrate_names)
# subs = substrates.get_subs(substrate_names)
corn = substrates.CORN
subs = [corn] + [corn.create_similar_substrate() for _ in range(2)]
xi = [sub.xi for sub in subs]
uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in subs]

# Simulation
n_steps = 336
t_step = 0.5 / 24  # Time in days

# MPC
n_horizon = 10

# Set up CHP
chp = utils.CHP(max_power=124.0)
chp_load = np.zeros(n_steps + n_horizon)
for i in range(12):
    chp_load[i::48] = 1.0  # 6:00 - 12:00
    chp_load[18 + i :: 48] = 1.0  # 15:00 - 21:00

vol_flow_rate = chp.compute_vol_flow_rate(
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
        0.0959467827166042,
        0.0125956088566735,
        4.64551507877156,
        0.850445630702126,
        957.771504116945,
        5.17942682112536,
        1.61926162860691e-09,
        1.49423713467601,
        0.627629518798448,
        1.96096023576704,
        0.535229147950488,
        13.9999999999977,
        0.0487500000000000,
        0.0956777093602567,
        4.22707761131655,
        0.0188914691525065,
        0.355199814936346,
        0.640296031589364,
        39.0,  # m^3
        36.0,  # m^3
    ]
)

# Normalisation
# Define numeric values for normalization with steady state
Tx = x0.copy()
feedVolFlowSS = 8.0 * 15000.0
Tu = np.array([feedVolFlowSS for _ in range(len(xi))])

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

xInNorm = list(np.array([val / Tx[:-2] for val in xi]).T)

x0[:-2] /= Tx[:-2]

# Set model
model = adm1_r3_frac_norm(xInNorm, Tu, Tx, Ty)

num_std_devs = 1  # Lower and upper bound of uncertainties is determined by the number of standard deviations that we consider

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
    x_ch_in=x_ch_in,
    x_pr_in=x_pr_in,
    x_li_in=x_li_in,
    vol_flow_rate=vol_flow_rate,
)

simulator = simulator_setup(
    model=model,
    t_step=t_step,
    x_ch_in=x_ch_in,
    x_pr_in=x_pr_in,
    x_li_in=x_li_in,
    vol_flow_rate=vol_flow_rate,
)
estimator = state_estimator.StateEstimator(model)


# Feeding
constant_feeding = False
# feed = np.zeros((n_steps, 1))

# feed[120:144] = 7.0 * 24 / Tu  # 42.0
# feed[264:312] = 3.0 * 24 / Tu  # 18.0

# feed[:, :] = feed[:, :]  # / 3.0

# Set x0
mpc.x0 = x0
simulator.x0 = x0

plot_vars = ["u_norm", "x_13", "x_14", "x_8", "x_19", "x_20"] + [
    f"y_{i+1}" for i in range(8)
]
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
for k in range(n_steps):
    timer.tic()
    u_norm_computed_old = copy.copy(u_norm_computed)
    if constant_feeding:
        # if feed.shape[1] == 1:
        #     u_norm = np.array([[feed[k]]])
        # else:
        u_norm_computed = np.array([feed[k]]).T
        # pass
    else:
        u_norm_computed = mpc.make_step(x0)

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

    print(u_norm_computed)
    print(u_norm_actual)
    timer.toc()

    y_next = simulator.make_step(u_norm_actual)
    x0 = estimator.make_step(y_next)

    if not disturbances.state_jumps is None:
        for x_idx, val in disturbances.state_jumps.items():
            if k == val[0]:
                x0[x_idx] += val[1]

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

if constant_feeding:
    fig, ax = plt.subplots(len(plot_vars), sharex=True)
    for idx, var in enumerate(plot_vars):
        ax[idx].set_ylabel(var)
        ax[idx].plot(results[var])

    plt.show()

timer.info()
timer.hist()

# input("Press any key to exit.")

# Store results:
if store_results:
    # do_mpc.data.save_results([mpc, simulator], f"testing_results")
    np.savetxt("results.csv", simulator.data._x, delimiter=",")
