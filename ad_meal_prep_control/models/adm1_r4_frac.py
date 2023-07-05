import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from params_R4 import *
import substrates

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def adm1_r4_frac(xi: np.ndarray):
    """
    ADM1-R4-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")
    num_inputs = xi.shape[1]

    # Time-invariant parameters
    c = np.array(
        [
            1 / Vl,
            kla,
            kla * Kch4 * R * T,
            kla * Kco2 * R * T,
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
            Vl / Vg,
            -kp / (Vg * p0) * (ph2o - p0) * ph2o,
        ]
    )

    # Stoichiometric constants
    a = np.array(
        [
            [
                0.2482,
                0.6809,
                0.0207,
                0.0456,
                1.0,
                np.nan,
                np.nan,
                np.nan,
                0.1372,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                0.2482,
                0.6809,
                0.0207,
                0.0456,
                np.nan,
                1.0,
                np.nan,
                np.nan,
                0.1372,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                0.3221,
                0.7954,
                0.1689,
                0.4588,
                np.nan,
                np.nan,
                1.0,
                np.nan,
                0.1723,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                0.6393,
                0.5817,
                0.0344,
                0.4152,
                np.nan,
                np.nan,
                np.nan,
                1.0,
                0.2286,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                0.18,
                np.nan,
                0.77,
                0.05,
                1.0,
                np.nan,
                np.nan,
                np.nan,
            ],
            [
                1.0,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                c[21],
                np.nan,
            ],
            [
                np.nan,
                1.0,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                c[21],
            ],
        ]
    ).transpose()

    a[3, :] = a[3, :] / 1000.0
    a[9, :] = a[9, :] / 1000.0

    # Input
    u = model.set_variable(var_type="_u", var_name="u", shape=(num_inputs, 1))

    # Time-variant parameters
    theta = np.array([kchF, kchS, kpr, kli, kdec, fracChFast])

    # inlet concentrations Mais A siehe Sörens Diss. [g/l] bzw. [mol/l]
    # S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ash, S_ch4_g, S_co2_g
    # xi = np.array([0, 0, 0.592, 960.512, 23.398, 0, 4.75, 1.381, 0, 17, 0, 0])  # [g/l]
    # xi = substrates.xi_values("R4-frac", "corn")

    # States
    x = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(14)
    ]

    # model.set_meas(
    #     "y_1",
    #     c[5] * x[10] ** 2
    #     + c[6] * x[10] * x[11]
    #     + c[7] * x[11] ** 2
    #     + c[8] * x[10]
    #     + c[9] * x[11]
    #     + c[10],
    # )

    # Set measurement equations as expressions
    model.set_expression(
        "y_1",
        (
            c[5] * x[10] ** 2
            + c[6] * x[10] * x[11]
            + c[7] * x[11] ** 2
            + c[8] * x[10]
            + c[9] * x[11]
            + c[10]
        )
        / 24.0,
    )
    model.set_expression("y_2", c[11] * x[10])
    model.set_expression("y_3", c[12] * x[11])
    model.set_expression("y_4", x[2])
    model.set_expression("y_5", 1.0 - 1.0 / c[13] * x[3])
    model.set_expression("y_6", 1.0 - 1.0 / (c[13] - x[3]) * x[9])

    if num_inputs > 1:
        sum_u = cumsum(u)[1]
    else:
        sum_u = u

    # Differential equations
    model.set_rhs(
        "x_1",
        c[0] * (xi[0, :] @ u - sum_u * x[0])
        + a[0][0] * theta[0] * x[4]
        + a[0][1] * theta[1] * x[5]
        + a[0][2] * theta[2] * x[6]
        + a[0][3] * theta[3] * x[7]
        - c[1] * x[0]
        + c[2] * x[10],
    )
    model.set_rhs(
        "x_2",
        c[0] * (xi[1, :] @ u - sum_u * x[1])
        + a[1][0] * theta[0] * x[4]
        + a[1][1] * theta[1] * x[5]
        + a[1][2] * theta[2] * x[6]
        + a[1][3] * theta[3] * x[7]
        - c[1] * x[1]
        + c[3] * x[11],
    )
    model.set_rhs(
        "x_3",
        c[0] * (xi[2, :] @ u - sum_u * x[2])
        - a[2][0] * theta[0] * x[4]
        - a[2][1] * theta[1] * x[5]
        + a[2][2] * theta[2] * x[6]
        - a[2][3] * theta[3] * x[7],
    )
    model.set_rhs(
        "x_4",
        c[0] * (xi[3, :] @ u - sum_u * x[3])
        - a[3][0] * theta[0] * x[4]
        - a[3][1] * theta[1] * x[5]
        - a[3][2] * theta[2] * x[6]
        - a[3][3] * theta[3] * x[7],
    )
    model.set_rhs(
        "x_5",
        c[0] * (theta[5] * xi[4, :] @ u - sum_u * x[4])
        - theta[0] * x[4]
        + a[4][4] * theta[4] * x[8],
    )
    model.set_rhs(
        "x_6",
        c[0] * ((1.0 - theta[5]) * xi[4, :] @ u - sum_u * x[5]) - theta[1] * x[5],
    )
    model.set_rhs(
        "x_7",
        c[0] * (xi[6, :] @ u - sum_u * x[6])
        - theta[2] * x[6]
        + a[6][4] * theta[4] * x[8],
    )
    model.set_rhs(
        "x_8",
        c[0] * (xi[7, :] @ u - sum_u * x[7])
        - theta[3] * x[7]
        + a[7][4] * theta[4] * x[8],
    )
    model.set_rhs(
        "x_9",
        c[0] * (xi[8, :] @ u - sum_u * x[8])
        + a[8][0] * theta[0] * x[4]
        + a[8][1] * theta[1] * x[5]
        + a[8][2] * theta[2] * x[6]
        + a[8][3] * theta[3] * x[7]
        - theta[4] * x[8],
    )
    model.set_rhs("x_10", c[0] * (xi[9, :] @ u - sum_u * x[9]))
    model.set_rhs(
        "x_11",
        c[14] * x[10] ** 3
        + c[15] * x[10] ** 2 * x[11]
        + c[16] * x[10] * x[11] ** 2
        + c[17] * x[10] ** 2
        + c[18] * x[10] * x[11]
        + c[19] * x[10]
        + c[4] * x[0],
    )
    model.set_rhs(
        "x_12",
        c[16] * x[11] ** 3
        + c[15] * x[10] * x[11] ** 2
        + c[14] * x[10] ** 2 * x[11]
        + c[18] * x[11] ** 2
        + c[17] * x[10] * x[11]
        + c[20] * x[11]
        + c[4] * x[1],
    )
    model.set_rhs(
        "x_13",
        model.aux["y_1"] * model.aux["y_2"] / (model.aux["y_2"] + model.aux["y_3"]),
    )
    model.set_rhs(
        "x_14",
        model.aux["y_1"] * model.aux["y_3"] / (model.aux["y_2"] + model.aux["y_3"]),
    )

    # Build the model
    model.setup()

    return model
