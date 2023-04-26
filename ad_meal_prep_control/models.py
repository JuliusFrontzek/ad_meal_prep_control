import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def adm1_r3_frac():
    """
    ADM1-R3-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")

    # Time-invariant parameters
    # Read from some file
    # Temporary: random numbers
    c = [model.set_variable(var_type="_p", var_name=f"c_{i+1}") for i in range(33)]

    # Stoichiometric constants
    a = [
        [
            0.6555,
            0.0818,
            0.2245,
            -0.0169,
            -0.0574,
            -1.0,
            np.nan,
            np.nan,
            np.nan,
            0.1125,
            np.nan,
        ],
        [
            0.6555,
            0.0818,
            0.2245,
            -0.0169,
            -0.0574,
            np.nan,
            -1.0,
            np.nan,
            np.nan,
            0.1125,
            np.nan,
        ],
        [
            0.9947,
            0.0696,
            0.1029,
            0.1746,
            -0.4767,
            np.nan,
            np.nan,
            -1.0,
            np.nan,
            0.1349,
            np.nan,
        ],
        [
            1.7651,
            0.1913,
            -0.6472,
            -0.0244,
            -0.4470,
            np.nan,
            np.nan,
            np.nan,
            -1.0,
            0.1621,
            np.nan,
        ],
        [
            -26.5447,
            6.7367,
            18.4808,
            -0.1506,
            0.4778,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1.0,
        ],
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.18,
            np.nan,
            0.77,
            0.05,
            -1.0,
            np.nan,
        ],
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            0.18,
            np.nan,
            0.77,
            0.05,
            np.nan,
            -1.0,
        ],
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1.0,
            np.nan,
            np.nan,
            np.nan,
        ],
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            -1.0,
            np.nan,
            np.nan,
        ],
        [
            np.nan,
            -1.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            c[31],
            np.nan,
        ],
        [
            np.nan,
            np.nan,
            -1.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            c[31],
            np.nan,
        ],
    ]

    # Input
    u = model.set_variable(var_type="_u", var_name="u")

    # Time-variant parameters
    theta = [
        model.set_variable(var_type="_tvp", var_name=f"theta_{i+1}", shape=(1, 1))
        for i in range(9)
    ]
    xi = [
        model.set_variable(var_type="_tvp", var_name=f"xi_{i+1}", shape=(1, 1))
        for i in range(13)
    ]

    # States
    x = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(18)
    ]

    # Aux
    phi = x[12] + (x[3] - x[15]) / 17.0 - x[14] / 44.0 - x[13] / 60.0

    s_h_plus = -phi / 2.0 + 0.5 * np.sqrt(phi**2 + c[3])

    i_ac = (
        c[2]
        / (c[2] + s_h_plus ** c[1])
        * x[3]
        / (x[3] + c[7])
        * theta[7]
        / (theta[7] + x[15])
    )

    # Differential equations
    model.set_rhs(
        "x_1",
        c[0] * (xi[0] - x[0]) * u
        + a[0][0] * theta[0] * x[5]
        + a[0][1] * theta[1] * x[6]
        + a[0][2] * theta[2] * x[7]
        + a[0][3] * theta[3] * x[8]
        - a[0][4] * theta[5] * x[0] * x[9] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_2",
        c[0] * (xi[1] - x[1]) * u
        + a[1][0] * theta[0] * x[5]
        + a[1][1] * theta[1] * x[6]
        + a[1][2] * theta[2] * x[7]
        + a[1][3] * theta[3] * x[8]
        + a[1][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac
        - c[4] * x[1]
        + c[5] * x[16],
    )
    model.set_rhs(
        "x_3",
        c[0] * (xi[2] - x[2]) * u
        + a[2][0] * theta[0] * x[5]
        + a[2][1] * theta[1] * x[6]
        + a[2][2] * theta[2] * x[7]
        - a[2][3] * theta[3] * x[8]
        + a[2][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac
        - c[4] * x[2]
        + c[4] * x[14]
        + c[6] * x[17],
    )

    model.set_rhs("x_12", c[0] * (xi[11] - x[11]) * u)
    model.set_rhs("x_13", c[0] * (xi[12] - x[12]) * u)
    model.set_rhs("x_14", c[28] * (x[0] - x[13]) - c[8] * x[13] * s_h_plus)
    model.set_rhs("x_15", c[29] * (x[2] - x[14]) - c[9] * x[14] * s_h_plus)
    model.set_rhs("x_16", c[30] * (x[3] - x[15]) - c[10] * x[15] * s_h_plus)

    # Build the model
    model.setup()

    return model
