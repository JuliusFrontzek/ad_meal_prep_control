import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
import random
from typing import Union

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def simulator_setup(
    model: do_mpc.model.Model,
    t_step: float,
    x_ch_in: np.ndarray,
    x_pr_in: np.ndarray,
    x_li_in: np.ndarray,
    vol_flow_rate: Union[np.ndarray, None] = None,
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
            t_now_idx = int(t_now / t_step)
            if vol_flow_rate is not None:
                tvp_num["v_ch4_dot_out", 0] = vol_flow_rate[t_now_idx]
            else:
                tvp_num["v_ch4_dot_out", 0] = 0.0
            return tvp_num

        simulator.set_tvp_fun(tvp_fun)

    # uncertain parameter realization in simulator -> drawn from uniform distribution
    # in between min and max values of uncertainties
    p_num = simulator.get_p_template()
    p_num["x_ch_in"] = np.array(
        [
            random.uniform(np.min(x_ch_in[:, i]), np.max(x_ch_in[:, i]))
            for i in range(x_ch_in.shape[1])
        ]
    )
    p_num["x_pr_in"] = np.array(
        [
            random.uniform(np.min(x_pr_in[:, i]), np.max(x_pr_in[:, i]))
            for i in range(x_pr_in.shape[1])
        ]
    )
    p_num["x_li_in"] = np.array(
        [
            random.uniform(np.min(x_li_in[:, i]), np.max(x_li_in[:, i]))
            for i in range(x_li_in.shape[1])
        ]
    )

    def p_fun(t_now):
        return p_num

    simulator.set_p_fun(p_fun)

    simulator.setup()

    return simulator
