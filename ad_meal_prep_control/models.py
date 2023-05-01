import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from params import *

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
            kp / p0 * (R * T / Mco2),
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
                0.0818,
                0.2245,
                0.0169,
                0.0574 / 1000.0,
                1.0,
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
                0.0169,
                0.0574 / 1000.0,
                np.nan,
                1.0,
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
                0.4767 / 1000.0,
                np.nan,
                np.nan,
                1.0,
                np.nan,
                0.1349,
                np.nan,
            ],
            [
                1.7651,
                0.1913,
                0.6472,
                0.0244,
                0.4470 / 1000.0,
                np.nan,
                np.nan,
                np.nan,
                1.0,
                0.1621,
                np.nan,
            ],
            [
                26.5447,
                6.7367,
                18.4808,
                0.1506,
                0.4778 / 1000.0,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                1.0,
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
                1.0,
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
                1.0,
            ],
            [
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                1.0,
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
                1.0,
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
                1.0,
                np.nan,
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
                c[31],
                np.nan,
            ],
            [
                np.nan,
                np.nan,
                1.0,
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
    ).transpose()

    # Input
    u = model.set_variable(var_type="_u", var_name="u")

    # Time-variant parameters
    theta = np.array([kchF, kchS, kpr, kli, kdec, mu_m_ac, K_S_ac, K_I_nh3, fracChFast])

    # inlet concentrations Mais A siehe Sörens Diss. [g/l] bzw. [mol/l]
    # S_ac, S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ac, X_ash, S_ion, S_ac-, S_hco3-, S_nh3, S_ch4_g, S_co2_g
    xi = np.array(
        [
            2.14500583382249,
            0,
            0,
            0.592294500000000,
            960.511750000000 / 1000.0,
            21.2528071170000,
            0,
            4.74931351500000,
            1.38096405000000,
            0,
            0,
            0.0562500000000000 / 1000.0,
            0.00750000000000000,
            0,
            0,
            0,
            0,
            0,
        ]
    )  # [g/l]

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
    model.set_expression("y_4", s_h_plus)
    model.set_expression("y_5", x[3])
    model.set_expression("y_6", 1.0 - 1.0 / c[20] * x[4])
    model.set_expression("y_7", 1.0 - 1.0 / (c[20] - x[4]) * x[11])
    model.set_expression("y_8", x[0])

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
    model.set_rhs(
        "x_4",
        c[0] * (xi[3] - x[3]) * u
        - a[3][0] * theta[0] * x[5]
        - a[3][1] * theta[1] * x[6]
        + a[3][2] * theta[2] * x[7]
        - a[3][3] * theta[3] * x[8]
        - a[3][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_5",
        c[0] * (xi[4] - x[4]) * u
        - a[4][0] * theta[0] * x[5]
        - a[4][1] * theta[1] * x[6]
        - a[4][2] * theta[2] * x[7]
        - a[4][3] * theta[3] * x[8]
        + a[4][4] * theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac,
    )
    model.set_rhs(
        "x_6",
        c[0] * (theta[8] * xi[5] - x[5]) * u
        - theta[0] * x[5]
        + a[5][5] * theta[4] * x[9]
        + a[5][6] * theta[4] * x[10],
    )
    model.set_rhs("x_7", c[0] * ((1.0 - theta[8]) * xi[5] - x[6]) * u - theta[1] * x[6])
    model.set_rhs(
        "x_8",
        c[0] * (xi[7] - x[7]) * u
        - theta[2] * x[7]
        + a[7][5] * theta[4] * x[9]
        + a[7][6] * theta[4] * x[10],
    )
    model.set_rhs(
        "x_9",
        c[0] * (xi[8] - x[8]) * u
        - theta[3] * x[8]
        + a[8][5] * theta[4] * x[9]
        + a[8][6] * theta[4] * x[10],
    )
    model.set_rhs(
        "x_10",
        c[0] * (xi[9] - x[9]) * u
        + a[9][0] * theta[0] * x[5]
        + a[9][1] * theta[1] * x[6]
        + a[9][2] * theta[2] * x[7]
        + a[9][3] * theta[3] * x[8]
        - theta[4] * x[9],
    )
    model.set_rhs(
        "x_11",
        c[0] * (xi[10] - x[10]) * u
        + theta[5] * x[0] * x[10] / (theta[6] + x[0]) * i_ac
        - theta[4] * x[10],
    )
    model.set_rhs("x_12", c[0] * (xi[11] - x[11]) * u)
    model.set_rhs("x_13", c[0] * (xi[12] - x[12]) * u)
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
        - c[11] * x[13]
        + c[27] * x[17],
    )
    # Build the model
    model.setup()

    return model


def adm1_r4_frac():
    """
    ADM1-R4-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")

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
            kp / p0 * (R * T / Mco2),
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

    # Input
    u = model.set_variable(var_type="_u", var_name="u")

    # Time-variant parameters
    theta = np.array([kchF, kchS, kpr, kli, kdec, mu_m_ac, K_S_ac, K_I_nh3, fracChFast])

    # inlet concentrations Mais A siehe Sörens Diss. [g/l] bzw. [mol/l]
    # S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ash, S_ch4_g, S_co2_g
    xi = np.array([0, 0, 0.592, 960.512, 23.398, 0, 4.75, 1.381, 0, 17, 0, 0])  # [g/l]

    # States
    x = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(12)
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
        c[5] * x[10] ** 2
        + c[6] * x[10] * x[11]
        + c[7] * x[11] ** 2
        + c[8] * x[10]
        + c[9] * x[11]
        + c[10],
    )
    model.set_expression("y_2", c[11] * x[10])
    model.set_expression("y_3", c[12] * x[11])
    model.set_expression("y_4", x[2])
    model.set_expression("y_5", 1.0 - 1.0 / c[13] * x[3])
    model.set_expression("y_6", 1.0 - 1.0 / (c[13] - x[3]) * x[9])

    # Differential equations
    model.set_rhs(
        "x_1",
        c[0] * (xi[0] - x[0]) * u
        + a[0][0] * theta[0] * x[4]
        + a[0][1] * theta[1] * x[5]
        + a[0][2] * theta[2] * x[6]
        + a[0][3] * theta[3] * x[7]
        - c[1] * x[0]
        + c[2] * x[10],
    )
    model.set_rhs(
        "x_2",
        c[0] * (xi[1] - x[1]) * u
        + a[1][0] * theta[0] * x[4]
        + a[1][1] * theta[1] * x[5]
        + a[1][2] * theta[2] * x[6]
        + a[1][3] * theta[3] * x[7]
        - c[1] * x[1]
        + c[3] * x[11],
    )
    model.set_rhs(
        "x_3",
        c[0] * (xi[2] - x[2]) * u
        - a[2][0] * theta[0] * x[4]
        - a[2][1] * theta[1] * x[5]
        + a[2][2] * theta[2] * x[6]
        - a[2][3] * theta[3] * x[7],
    )
    model.set_rhs(
        "x_4",
        c[0] * (xi[3] - x[3]) * u
        - a[3][0] * theta[0] * x[4]
        - a[3][1] * theta[1] * x[5]
        - a[3][2] * theta[2] * x[6]
        - a[3][3] * theta[3] * x[7],
    )
    model.set_rhs(
        "x_5",
        c[0] * (theta[5] * xi[4] - x[4]) * u
        - theta[0] * x[4]
        + a[4][4] * theta[4] * x[8],
    )
    model.set_rhs("x_6", c[0] * ((1.0 - theta[5]) * xi[4] - x[5]) * u - theta[1] * x[5])
    model.set_rhs(
        "x_7", c[0] * (xi[6] - x[6]) * u - theta[2] * x[6] + a[6][4] * theta[4] * x[8]
    )
    model.set_rhs(
        "x_8", c[0] * (xi[7] - x[7]) * u - theta[3] * x[7] + a[7][4] * theta[4] * x[8]
    )
    model.set_rhs(
        "x_9",
        c[0] * (xi[8] - x[8]) * u
        + a[8][0] * theta[0] * x[4]
        + a[8][1] * theta[1] * x[5]
        + a[8][2] * theta[2] * x[6]
        + a[8][3] * theta[3] * x[7]
        - theta[4] * x[8],
    )
    model.set_rhs("x_10", c[0] * (xi[9] - x[9]) * u)
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
    # Build the model
    model.setup()

    return model
