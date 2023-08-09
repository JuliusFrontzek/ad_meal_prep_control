import numpy as np
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
import subprocess

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def mpc_setup(model):
    num_states = 20

    mpc = do_mpc.controller.MPC(model)

    setup_mpc = {
        "n_horizon": 30,
        "n_robust": 0,
        "open_loop": 0,
        "t_step": 0.5 / 24,
        "state_discretization": "collocation",
        "collocation_type": "radau",
        "collocation_deg": 2,
        "collocation_ni": 1,
        # Use MA27 linear solver in ipopt for faster calculations:
        # "nlpsol_opts": {
        #     "ipopt.linear_solver": "MA57",
        #     "ipopt.hsllib": "/home/julius/Documents/coinhsl-2023.05.26/builddir/libcoinhsl.so",
        # },
    }

    mpc.set_param(**setup_mpc)

    mpc.set_param(store_full_solution=True)

    # Scaling of units for better conditioning of optimization problem
    # mpc.scaling["_u", "u"] = 100

    mterm = (
        model.x[f"x_{14}"] - 0.95
    ) ** 2  # + (model.x[f"x_{num_states}"] - 0.8) ** 2
    lterm = (
        model.x[f"x_{14}"] - 0.95
    ) ** 2  # + (        model.x[f"x_{num_states}"] - 0.8    ) ** 2
    # mterm = (model.aux["y_2"] - 0.4) ** 2
    # lterm = (model.aux["y_2"] - 0.4) ** 2
    mpc.set_objective(lterm=lterm, mterm=mterm)

    # mpc.set_rterm(u=0.1)

    # # Hard constraints
    mpc.bounds["lower", "_u", "u_norm"] = 0.0
    mpc.bounds["upper", "_u", "u_norm"] = 10.0

    for i in range(num_states):
        mpc.bounds["lower", "_x", f"x_{i+1}"] = 0.0

    # # Soft constraints with slack variables
    # mpc.set_nl_cons(
    #     "some_constraint",
    #     model.x["x_2"],
    #     ub=50.0,
    #     soft_constraint=True,
    #     penalty_term_cons=1e2,
    # )

    mpc.setup()

    return mpc
