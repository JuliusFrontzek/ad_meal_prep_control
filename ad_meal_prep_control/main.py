import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
import random
import ad_meal_prep_control.state_estimator as state_estimator

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc
from do_mpc.tools import Timer

import matplotlib.pyplot as plt
import pickle
import time


from models import adm1_r3_frac, adm1_r4_frac
from mpc import mpc_setup
from simulator import template_simulator

""" User settings: """
show_animation = True
store_results = False

# Model
model_type = "R4"
if model_type == "R3":
    # Set model
    model = adm1_r3_frac()
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.0494667574155131,
            0.0116512808544296,
            4.97521803226548,
            0.963856429890969,
            957.102301745169 / 1000.0,
            0.740449903041190,
            2.22134970912357,
            0.948575266222027,
            0.412040527872582,
            1.92558575222279,
            0.521526149938689,
            0.0562500000000000 / 1000.0,
            0.00750000000000000,
            0.0493342637010683,
            4.54552248672517,
            0.0223970128000483,
            0.358267793064052,
            0.660494133806800,
        ]
    )
elif model_type == "R4":
    # Set model
    model = adm1_r4_frac()
    # Set the initial state of mpc and simulator
    x0 = np.array(
        [
            0.091,
            0.508,
            0.944,
            956.97 / 1000.0,
            0.5 * 3.26,
            0.5 * 3.26,
            0.956,
            0.413,
            2.569,
            1,
            0.315,
            0.78,
        ]
    )
else:
    raise NotImplementedError(f"Model '{model_type}' not implemented.")

mpc = mpc_setup(model, model_type)
simulator = template_simulator(model)
estimator = state_estimator.StateEstimator(model)

# Simulation
n_steps = 50

# Feeding
constant_feeding = False
feed = 42.0 * np.ones(n_steps)


mpc.x0 = x0
simulator.x0 = x0

if constant_feeding:
    fig, ax = plt.subplots(6, sharex=True)
    ax[0].set_ylabel("")
else:
    mpc.set_initial_guess()

    # Initialize graphic:
    graphics = do_mpc.graphics.Graphics(mpc.data)

    # Configure plot:
    plot_vars = ("x_10", "x_11", "u")
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
    y_next = simulator.make_step(u0)
    x0 = estimator.make_step(y_next)

    if show_animation:
        if constant_feeding:
            pass
        else:
            graphics.plot_results(t_ind=k)
            graphics.plot_predictions(t_ind=k)
        graphics.reset_axes()
        plt.show()
        plt.pause(0.01)

timer.info()
timer.hist()

input("Press any key to exit.")

# Store results:
if store_results:
    do_mpc.data.save_results([mpc, simulator], "admr1-r3-frac_results")
