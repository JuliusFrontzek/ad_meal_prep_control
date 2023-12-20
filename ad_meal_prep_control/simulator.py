import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from ad_meal_prep_control.utils import Disturbances

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
    disturbances: Disturbances,
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

    tvp_template = simulator.get_tvp_template()

    def dictated_sub_tvp_setup(t_now: float):
        if disturbances.dictated_feeding is not None:
            for feed_idx, dictated_sub in enumerate(
                disturbances.dictated_feeding.values()
            ):
                if t_now >= dictated_sub[0] and t_now < dictated_sub[1]:
                    tvp_template["dictated_sub_feed", feed_idx] = dictated_sub[2]
                else:
                    tvp_template["dictated_sub_feed", feed_idx] = 0.0

    if num_states == 20:  # i.e. if we consider the gas storage

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            tvp_template["v_ch4_dot_tank_out", 0] = ch4_outflow_rate[t_now_idx]
            dictated_sub_tvp_setup(t_now)
            return tvp_template

    else:

        def tvp_fun(t_now):
            dictated_sub_tvp_setup(t_now)
            return tvp_template

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
