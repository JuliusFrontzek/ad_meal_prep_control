import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from ad_meal_prep_control.params_R3 import *

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def adm1_r3_frac(xi: np.ndarray):
    """
    ADM1-R3-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")
    num_inputs = xi.shape[1]

    # Time-invariant parameters
    c = np.array(
        [
            1 / Vl,
            nac,
            10 ** (-1.5 * (pHULac + pHLLac) / (pHULac - pHLLac)),
            4 * KW,
            kla,
            kla * Kch4 * R * T,
            kla * Kco2 * R * T,
            KS_IN,
            k_AB_ac,
            k_AB_co2,
            k_AB_IN,
            kla * Vl / Vg,
            kp / p0 * (R * T / Mch4) ** 2,
            2 * kp / p0 * (R * T) ** 2 / (Mch4 * Mco2),
            kp / p0 * (R * T / Mco2) ** 2,
            kp / p0 * R * T / Mch4 * (2 * ph2o - p0),
            kp / p0 * R * T / Mco2 * (2 * ph2o - p0),
            kp / p0 * (ph2o - p0) * ph2o,
            R * T / Mch4,
            R * T / Mco2,
            rho,
            -kp / (Vg * p0) * (R * T / Mch4) ** 2,
            -2 * kp / (Vg * p0) * (R * T) ** 2 / (Mch4 * Mco2),
            -kp / (Vg * p0) * (R * T / Mco2) ** 2,
            -kp / (Vg * p0) * R * T / Mch4 * (2 * ph2o - p0),
            -kp / (Vg * p0) * R * T / Mco2 * (2 * ph2o - p0),
            -kla * Vl / Vg * Kch4 * R * T - kp / (Vg * p0) * (ph2o - p0) * ph2o,
            -kla * Vl / Vg * Kco2 * R * T - kp / (Vg * p0) * (ph2o - p0) * ph2o,
            k_AB_ac * K_a_ac,
            k_AB_co2 * K_a_co2,
            k_AB_IN * K_a_IN,
            Vl / Vg,
            -kp / (Vg * p0) * (ph2o - p0) * ph2o,
        ]
    )

    # Stoichiometric constants
    a = np.array(
        [
            [
                0.6555,
                0.081837,
                0.2245,
                0.016932,
                0.057375,
                -1,
                np.nan,
                np.nan,
                np.nan,
                0.11246,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                0.6555,
                0.081837,
                0.2245,
                0.016932,
                0.057375,
                np.nan,
                -1,
                np.nan,
                np.nan,
                0.11246,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                0.9947,
                0.069636,
                0.10291,
                0.17456,
                0.47666,
                np.nan,
                np.nan,
                -1,
                np.nan,
                0.13486,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                1.7651,
                0.19133,
                0.64716,
                0.024406,
                0.44695,
                np.nan,
                np.nan,
                np.nan,
                -1,
                0.1621,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                26.5447,
                6.7367,
                18.4808,
                0.15056,
                0.4778,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                1,
                np.nan,
                np.nan,
                np.nan,
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
                0.18,
                np.nan,
                0.77,
                0.05,
                -1,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
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
                0.18,
                np.nan,
                0.77,
                0.05,
                np.nan,
                -1,
                np.nan,
                np.nan,
                np.nan,
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
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                -1,
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
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                -1,
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
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                -1,
                np.nan,
                np.nan,
            ],
            [
                np.nan,
                -1,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
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
                -1,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                c[31],
            ],
        ]
    ).transpose()

    # Input
    u = model.set_variable(var_type="_u", var_name="u", shape=(num_inputs, 1))

    # Time-variant parameters
    theta = np.array([kchF, kchS, kpr, kli, kdec, mu_m_ac, K_S_ac, K_I_nh3, fracChFast])

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

    # Set measurement equations as expressions
    model.set_expression(
        "y_1",
        c[12] * x[16] ** 2
        + c[13] * x[16] * x[17]
        + c[14] * x[17] ** 2
        + c[15] * x[16]
        + c[16] * x[17]
        + c[17],
    )
    model.set_expression("y_2", c[18] * x[16])
    model.set_expression("y_3", c[19] * x[17])
    model.set_expression("y_4", -np.log10(s_h_plus))
    model.set_expression("y_5", x[3])
    model.set_expression("y_6", 1.0 - 1.0 / c[20] * x[4])
    model.set_expression("y_7", 1.0 - 1.0 / (c[20] - x[4]) * x[11])
    model.set_expression("y_8", x[0])

    if num_inputs > 1:
        sum_u = cumsum(u)[1]
    else:
        sum_u = u

    # Differential equations
    model.set_rhs(
        "x_1",
        c[0] * (xi[0, :] @ u - sum_u * x[0])
        + a[0][0] * theta[0] * x[5]
        + a[0][1] * theta[1] * x[6]
        + a[0][2] * theta[2] * x[7]
        + a[0][3] * theta[3] * x[8]
        - a[0][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_2",
        c[0] * (xi[1, :] @ u - sum_u * x[1])
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
        c[0] * (xi[2, :] @ u - sum_u * x[2])
        + a[2][0] * theta[0] * x[5]
        + a[2][1] * theta[1] * x[6]
        + a[2][2] * theta[2] * x[7]
        - a[2][3] * theta[3] * x[8]
        + a[2][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac
        - c[4] * x[2]
        + c[4] * x[14]
        + c[6] * x[17],
    )
    model.set_rhs(
        "x_4",
        c[0] * (xi[3, :] @ u - sum_u * x[3])
        - a[3][0] * theta[0] * x[5]
        - a[3][1] * theta[1] * x[6]
        + a[3][2] * theta[2] * x[7]
        - a[3][3] * theta[3] * x[8]
        - a[3][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_5",
        c[0] * (xi[4, :] @ u - sum_u * x[4])
        - a[4][0] * theta[0] * x[5]
        - a[4][1] * theta[1] * x[6]
        - a[4][2] * theta[2] * x[7]
        - a[4][3] * theta[3] * x[8]
        + a[4][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_6",
        c[0] * (theta[8] * xi[5, :] @ u - sum_u * x[5])
        - theta[0] * x[5]
        + a[5][5] * theta[4] * x[9]
        + a[5][6] * theta[4] * x[10],
    )
    model.set_rhs(
        "x_7", c[0] * ((1.0 - theta[8]) * xi[5, :] @ u - sum_u * x[6]) - theta[1] * x[6]
    )
    model.set_rhs(
        "x_8",
        c[0] * (xi[7, :] @ u - sum_u * x[7])
        - theta[2] * x[7]
        + a[7][5] * theta[4] * x[9]
        + a[7][6] * theta[4] * x[10],
    )
    model.set_rhs(
        "x_9",
        c[0] * (xi[8, :] @ u - sum_u * x[8])
        - theta[3] * x[8]
        + a[8][5] * theta[4] * x[9]
        + a[8][6] * theta[4] * x[10],
    )
    model.set_rhs(
        "x_10",
        c[0] * (xi[9, :] @ u - sum_u * x[9])
        + a[9][0] * theta[0] * x[5]
        + a[9][1] * theta[1] * x[6]
        + a[9][2] * theta[2] * x[7]
        + a[9][3] * theta[3] * x[8]
        - theta[4] * x[9],
    )
    model.set_rhs(
        "x_11",
        c[0] * (xi[10, :] @ u - sum_u * x[10])
        + theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac
        - theta[4] * x[10],
    )
    model.set_rhs("x_12", c[0] * (xi[11, :] @ u - sum_u * x[11]))
    model.set_rhs("x_13", c[0] * (xi[12, :] @ u - sum_u * x[12]))
    model.set_rhs("x_14", c[28] * (x[0] - x[13]) - c[8] * x[13] * s_h_plus)
    model.set_rhs("x_15", c[29] * (x[2] - x[14]) - c[9] * x[14] * s_h_plus)
    model.set_rhs("x_16", c[30] * (x[3] - x[15]) - c[10] * x[15] * s_h_plus)
    model.set_rhs(
        "x_17",
        c[21] * x[16] ** 3
        + c[22] * x[16] ** 2 * x[17]
        + c[23] * x[16] * x[17] ** 2
        + c[24] * x[16] ** 2
        + c[25] * x[16] * x[17]
        + c[11] * x[1]
        + c[26] * x[16],
    )
    model.set_rhs(
        "x_18",
        c[23] * x[17] ** 3
        + c[22] * x[16] * x[17] ** 2
        + c[21] * x[16] ** 2 * x[17]
        + c[25] * x[17] ** 2
        + c[24] * x[16] * x[17]
        + c[11] * x[2]
        - c[11] * x[14]
        + c[27] * x[17],
    )

    # Build the model
    model.setup()

    return model
