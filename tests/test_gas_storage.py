from ad_meal_prep_control.models.gas_storage import gas_storage_model, GasConditions

import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
import random
import ad_meal_prep_control.state_estimator as state_estimator
from ad_meal_prep_control.simulator import template_simulator
import unittest

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc
from do_mpc.tools import Timer

import matplotlib.pyplot as plt
import pickle
import time


""" User settings: """
show_animation = True
save_results = True

p0 = 1.0133  # [bar]
p_gas_storage = p0 * 1.001  # slightly larger than atmospheric pressure [bar]
p_norm = p0  # [bar]
T_norm = 273.15  # [K]
p_co2_phase_change = 0.3  # [bar]
p_ch4_phase_change = 0.3  # [bar]
p_h2o = 0.123440  # [bar] (Vapour pressure at 50°C)
T = 311  # operating temperature [K]

gas_conditions = GasConditions(
    p0=p0,
    p_gas_storage=p_ch4_phase_change,
    p_norm=p_norm,
    T_norm=T_norm,
    p_co2_phase_change=p_co2_phase_change,
    p_ch4_phase_change=p_ch4_phase_change,
    p_h2o=p_h2o,
    T=T,
)

model = gas_storage_model(gas_conditions)


simulator = template_simulator(model)
estimator = state_estimator.StateEstimator(model)

# Simulation
n_steps = 50

# Set x0
simulator.x0 = np.array([1.0, 0.0])

plot_vars = ["u", "x_1", "x_2"]

# Initialize graphic:
graphics = do_mpc.graphics.Graphics(simulator.data)

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

u = np.array([[10.0], [10.0]])

# for k in range(n_steps):
#     timer.tic()

#     timer.toc()

#     x_next = simulator.make_step(u)
#     x0 = estimator.make_step(x_next)

#     if show_animation:
#         graphics.plot_results(t_ind=k)
#         graphics.reset_axes()
#         plt.show()
#         plt.pause(0.01)

# timer.info()
# timer.hist()


class TestGasStorage(unittest.TestCase):
    def test_ch4_inflow(self):
        """
        Testing if a net inflow of CH4 into the tank, the CH4 volume inside the tank actually increases.
        """
        gas_conditions = GasConditions(
            p0=p0,
            p_gas_storage=p0,
            p_norm=p_norm,
            T_norm=T_norm,
            p_co2_phase_change=0.0,
            p_ch4_phase_change=p_ch4_phase_change,
            p_h2o=p_h2o,
            T=T_norm,
        )

        model = gas_storage_model(gas_conditions)
        simulator = template_simulator(model)

        x0 = np.array([1.0, 0.0])
        simulator.x0 = x0
        u = np.array([[10.0], [0.0]])

        for _ in range(n_steps):
            simulator.make_step(u)

        self.assertGreater(simulator.data._x[-1][0], x0[0])

    def test_ch4_outflow(self):
        pass
