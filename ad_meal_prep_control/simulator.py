import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
import random

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def simulator_setup(
    model: do_mpc.model.Model,
    t_step: float,
    xi_ch_norm: np.ndarray,
    xi_pr_norm: np.ndarray,
    xi_li_norm: np.ndarray,
    ch4_outflow_rate: np.ndarray,
):
    num_states = model._x.size
    simulator = do_mpc.simulator.Simulator(model)

    params_simulator = {
        "integration_tool": "cvodes",
        "abstol": 1e-10,
        "reltol": 1e-10,
        "t_step": t_step,
    }

    simulator.set_param(**params_simulator)

    if num_states == 20:  # i.e. if we consider the gas storage
        tvp_num = simulator.get_tvp_template()

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            tvp_num["v_ch4_dot_out", 0] = ch4_outflow_rate[t_now_idx]
            return tvp_num

        simulator.set_tvp_fun(tvp_fun)

    p_num = simulator.get_p_template()
    p_num["xi_ch_norm"] = xi_ch_norm
    p_num["xi_pr_norm"] = xi_pr_norm
    p_num["xi_li_norm"] = xi_li_norm

    def p_fun(t_now):
        return p_num

    simulator.set_p_fun(p_fun)

    simulator.setup()

    return simulator
