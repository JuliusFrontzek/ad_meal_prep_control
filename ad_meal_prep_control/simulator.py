import numpy as np
from casadi import *
from casadi.tools import *
import pdb
import sys
import os

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def template_simulator(model):
    simulator = do_mpc.simulator.Simulator(model)

    params_simulator = {
        "integration_tool": "cvodes",
        "abstol": 1e-10,
        "reltol": 1e-10,
        "t_step": 0.5 / 24,
    }

    simulator.set_param(**params_simulator)

    # tvp_num = simulator.get_tvp_template()

    # def tvp_fun(t_now):
    #     return tvp_num

    # simulator.set_tvp_fun(tvp_fun)

    # p_num = simulator.get_p_template()
    # p_num["alpha"] = 1
    # p_num["beta"] = 1

    # def p_fun(t_now):
    #     return p_num

    # simulator.set_p_fun(p_fun)

    simulator.setup()

    return simulator
