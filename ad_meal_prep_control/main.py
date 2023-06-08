import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
import random
import ad_meal_prep_control.state_estimator as state_estimator
import substrates

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc
from do_mpc.tools import Timer

import matplotlib.pyplot as plt
import pickle
import time


from models import adm1_r3_frac, adm1_r4_frac, adm1_r3_frac_norm
from mpc import mpc_setup
from simulator import template_simulator

""" User settings: """
show_animation = True
store_results = True

# Model
model_type = "R3-frac"
substrate_names = ["corn"]
xi = substrates.xi_values(model_type, substrate_names)

if model_type == "R3-frac":
    # Set model
    model = adm1_r3_frac(xi)
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.0949255130768584,
            0.0127895628562427,
            4.65470911882496,
            0.850440543095153,
            957.771520363428,
            5.17942716544083,
            0,
            1.49423886792478,
            0.627629693919376,
            1.96096081782644,
            0.535260086844425,
            13.9999968004443,
            0.0487500000000000,
            0.0946545714143219,
            4.22865312776708,
            0.0185676898303715,
            0.362017286980919,
            0.652085696312695,
        ]
    )
elif model_type == "R4-frac":
    # Set model
    model = adm1_r4_frac(xi)
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.0221697506048759,
            0.554580163344390,
            0.557189956577887,
            0.959887916032620,
            2.03934924735466,
            14.7843348334834,
            4.14214818444786,
            1.28561345217435,
            1.01590412485130,
            0.0170000000000000,
            0.386609337719045,
            0.924038131164312,
            0.0,
            0.0,
        ]
    )
elif model_type == "R3-frac-norm":
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.0949255126142516,
            0.0127895628548919,
            4.65470911344473,
            0.850440540811826,
            957.771520370543,
            5.17942716513992,
            1.61926162860698e-09,
            1.49423886908271,
            0.627629693991017,
            1.96096082868930,
            0.535260089706592,
            13.9999999999951,
            0.0487500000000000,
            0.0946545709535683,
            4.22865312234033,
            0.0185676897742391,
            0.362017286954370,
            0.652085696384954,
        ]
    )

    # Normalisation
    # Define numeric values for normalization with steady state
    Tx = x0
    feedVolFlowSS = 8
    Tu = feedVolFlowSS
    Ty = np.array([137.6596, 0.5851, 0.3832, 7.3033, 0.8319, 0.0422, 0.6685, 0.0949])

    # Normalize simulation inputs
    uNorm = feedVolFlowSS / Tu
    x0SSNorm = x0 / Tx
    xInNorm = xi / Tx

    # Set model
    model = adm1_r3_frac_norm(xi, xInNorm, uNorm, Tu, Tx, Ty)
else:
    raise NotImplementedError(f"Model '{model_type}' not implemented.")

mpc = mpc_setup(model, model_type)
simulator = template_simulator(model)
estimator = state_estimator.StateEstimator(model)

# Simulation
n_steps = 336

# Feeding
constant_feeding = True
feed = np.zeros(n_steps)

feed[120:144] = 42.0 * 24
feed[264:312] = 18.0 * 24

# Set x0
mpc.x0 = x0
simulator.x0 = x0

plot_vars = ["u", "x_13", "x_14", "x_8"] + [f"y_{i+1}" for i in range(8)]
# plot_vars = ["u"] + [f"x_{i+1}" for i in range(18)]

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

for k in range(n_steps):
    timer.tic()
    if constant_feeding:
        u0 = np.array([[feed[k]]])
    else:
        u0 = mpc.make_step(x0)
    timer.toc()

    x_next = simulator.make_step(u0)
    x0 = estimator.make_step(x_next)

    if show_animation:
        if constant_feeding:
            if "u" in plot_vars:
                results["u"][k] = u0

            for var in plot_vars:
                if var.startswith("x"):
                    results[var][k] = x_next[int(var.lstrip("x_")) - 1]
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
