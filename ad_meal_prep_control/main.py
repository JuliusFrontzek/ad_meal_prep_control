import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
import random

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc
from do_mpc.tools import Timer

import matplotlib.pyplot as plt
import pickle
import time


from models import adm1_r3_frac
from mpc import mpc_setup
from simulator import template_simulator

""" User settings: """
show_animation = True
store_results = False


model = adm1_r3_frac()
mpc = mpc_setup(model)
simulator = template_simulator(model)
estimator = do_mpc.estimator.StateFeedback(model)

# Set the initial state of mpc and simulator:
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

mpc.x0 = x0
simulator.x0 = x0

mpc.set_initial_guess()

# Initialize graphic:
graphics = do_mpc.graphics.Graphics(mpc.data)


fig, ax = plt.subplots(3, sharex=True)
# Configure plot:
graphics.add_line(var_type="_x", var_name="x_17", axis=ax[0])
graphics.add_line(var_type="_x", var_name="x_18", axis=ax[1])
graphics.add_line(var_type="_u", var_name="u", axis=ax[2])
ax[0].set_ylabel("x_17")
ax[1].set_ylabel("x_18")
ax[2].set_ylabel("u")
# Update properties for all prediction lines:
for line_i in graphics.pred_lines.full:
    line_i.set_linewidth(1)

# label_lines = graphics.result_lines["_x", "C_a"] + graphics.result_lines["_x", "C_b"]
# ax[0].legend(label_lines, ["C_a", "C_b"])
# label_lines = graphics.result_lines["_x", "T_R"] + graphics.result_lines["_x", "T_K"]
# ax[1].legend(label_lines, ["T_R", "T_K"])

fig.align_ylabels()
fig.tight_layout()
plt.ion()

timer = Timer()

for k in range(50):
    timer.tic()
    u0 = mpc.make_step(x0)
    timer.toc()
    y_next = simulator.make_step(u0)
    x0 = estimator.make_step(y_next)

    if show_animation:
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
