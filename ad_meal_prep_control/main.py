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
model_type = "R3-frac-norm"
substrate_names = ["corn"]
xi = substrates.xi_values(model_type, substrate_names)

if model_type == "R3-frac":
    # Set model
    model = adm1_r3_frac(xi)
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.0959467827268156,
            0.0125956088575309,
            4.64551507823071,
            0.850445630490791,
            957.771504117476,
            5.17942682165898,
            0,
            1.49423713477719,
            0.627629518812174,
            1.96096023676352,
            0.535229148231421,
            13.9999999815809,
            0.0487500000000000,
            0.0956777093701918,
            4.22707761079332,
            0.0188914691405196,
            0.355199814946965,
            0.640296031561076,
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
    x0Norm = np.array(
        [
            0.878474066127602,
            0.999999900040164,
            0.933730953851068,
            0.882336418918662,
            0.999999995095856,
            0.999999999988468,
            0.999984804739496,
            0.999999530600200,
            0.999998927072078,
            0.999996522412598,
            0.999996609020940,
            0.999985987995331,
            1.00000000000000,
            0.878406260523190,
            0.929943174578399,
            0.784593153260128,
            0.956768888763069,
            0.969419709921008,
        ]
    )

    # Normalisation
    # Define numeric values for normalization with steady state
    Tx = np.array(
        [
            0.109219823638302,
            0.0125956101155577,
            4.97521803226548,
            0.963856429890969,
            957.771508813996,
            5.17942682118517,
            1.61928623408307e-09,
            1.49423783607141,
            0.627630192200464,
            1.96096705520169,
            0.535230962909218,
            14.0001961708128,
            0.0487500000000000,
            0.108921934550271,
            4.54552248672517,
            0.0240780448734638,
            0.371249336286256,
            0.660494133806800,
        ]
    )
    feedVolFlowSS = 8.0
    Tu = feedVolFlowSS
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

    # Normalize simulation inputs
    uNorm = feedVolFlowSS / Tu
    # x0SSNorm = x0 / Tx
    xInNorm = xi / Tx

    x0 = x0Norm

    # Set model
    model = adm1_r3_frac_norm(xInNorm, uNorm, Tu, Tx, Ty)
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

feed[120:144] = 7.0 * 24 / Tu  # 42.0
feed[264:312] = 3.0 * 24 / Tu  # 18.0

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
