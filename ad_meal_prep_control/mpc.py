import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from params_R3 import *

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def mpc_setup(
    model: do_mpc.model.Model,
    n_horizon: int,
    t_step: float,
    x_ch_in: np.ndarray,
    x_pr_in: np.ndarray,
    x_li_in: np.ndarray,
    vol_flow_rate: np.ndarray,
    compile_nlp: bool,
):
    num_states = 20

    mpc = do_mpc.controller.MPC(model)

    setup_mpc = {
        "n_horizon": n_horizon,
        "n_robust": 0,
        "open_loop": 0,
        "t_step": t_step,
        "state_discretization": "collocation",
        "collocation_type": "radau",
        "collocation_deg": 2,
        "collocation_ni": 1,
        # Use MA27 linear solver in ipopt for faster calculations:
        # "nlpsol_opts": {
        #     "ipopt.linear_solver": "MA27",
        #     "ipopt.hsllib": "/home/julius/Documents/coinhsl-2023.05.26/builddir/libcoinhsl.so",
        # },
    }

    mpc.set_param(**setup_mpc)

    mpc.set_param(store_full_solution=True)

    tvp_template = mpc.get_tvp_template()

    def tvp_fun(t_now):
        t_now_idx = int(t_now / t_step)
        for k in range(n_horizon + 1):
            tvp_template["_tvp", k, "v_ch4_dot_out"] = vol_flow_rate[t_now_idx + k]

        return tvp_template

    mpc.set_tvp_fun(tvp_fun)

    # Scaling of units for better conditioning of optimization problem
    # mpc.scaling["_u", "u"] = 100

    # mterm = (model.x[f"x_{19}"] - 80.0) ** 2
    # lterm = (model.x[f"x_{19}"] - 80.0) ** 2
    mterm = (model.aux["y_4"] - 7.5) ** 2
    lterm = (model.aux["y_4"] - 7.5) ** 2
    mpc.set_objective(lterm=lterm, mterm=mterm)

    mpc.set_rterm(u_norm=0.1)

    # Hard constraints
    mpc.bounds["lower", "_u", "u_norm"] = 0.0
    mpc.bounds["upper", "_u", "u_norm"] = 1.0
    # mpc.bounds["lower", "_u", "u"] = 0.0
    # mpc.bounds["upper", "_u", "u"] = 10.0

    # NL constraints
    mpc.set_nl_cons(
        "max_vol_gas_storage",
        model._aux_expression["v_gas_storage"],
        ub=V_GAS_STORAGE_MAX,
        soft_constraint=True,
        penalty_term_cons=1e7,
    )

    # Bounds for states
    # for i in range(num_states):
    #     mpc.bounds["lower", "_x", f"x_{i+1}"] = 0.0

    mpc.bounds["lower", "_x", "x_19"] = 0.0
    mpc.bounds["lower", "_x", "x_20"] = 0.0

    # # Soft constraints with slack variables
    # mpc.set_nl_cons(
    #     "some_constraint",
    #     model.x["x_2"],
    #     ub=50.0,
    #     soft_constraint=True,
    #     penalty_term_cons=1e2,
    # )

    mpc.set_uncertainty_values(x_ch_in=x_ch_in, x_pr_in=x_pr_in, x_li_in=x_li_in)

    mpc.setup()
    if compile_nlp:
        mpc.compile_nlp(
            overwrite=True, compiler_command="gcc -fPIC -shared -O3 nlp.c -o nlp.so"
        )

    return mpc
