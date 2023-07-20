import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from params_R3 import *
import substrates

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def adm1_r3_frac_norm(
    xi_norm: np.ndarray,
    Tu: np.ndarray,
    Tx: np.ndarray,
    Ty: np.ndarray,
):
    """
    Normalized ADM1-R3-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")
    num_inputs = xi_norm.shape[1]

    # Time-invariant parameters
    # c = np.array(
    #     [
    #         1 / Vl,
    #         nac,
    #         10 ** (-1.5 * (pHULac + pHLLac) / (pHULac - pHLLac)),
    #         4 * KW,
    #         kla,
    #         kla * Kch4 * R * T,
    #         kla * Kco2 * R * T,
    #         KS_IN,
    #         k_AB_ac,
    #         k_AB_co2,
    #         k_AB_IN,
    #         kla * Vl / Vg,
    #         kp / p0 * (R * T / Mch4) ** 2,
    #         2 * kp / p0 * (R * T) ** 2 / (Mch4 * Mco2),
    #         kp / p0 * (R * T / Mco2) ** 2,
    #         kp / p0 * R * T / Mch4 * (2 * ph2o - p0),
    #         kp / p0 * R * T / Mco2 * (2 * ph2o - p0),
    #         kp / p0 * (ph2o - p0) * ph2o,
    #         R * T / Mch4,
    #         R * T / Mco2,
    #         rho,
    #         -kp / (Vg * p0) * (R * T / Mch4) ** 2,
    #         -2 * kp / (Vg * p0) * (R * T) ** 2 / (Mch4 * Mco2),
    #         -kp / (Vg * p0) * (R * T / Mco2) ** 2,
    #         -kp / (Vg * p0) * R * T / Mch4 * (2 * ph2o - p0),
    #         -kp / (Vg * p0) * R * T / Mco2 * (2 * ph2o - p0),
    #         -kla * Vl / Vg * Kch4 * R * T - kp / (Vg * p0) * (ph2o - p0) * ph2o,
    #         -kla * Vl / Vg * Kco2 * R * T - kp / (Vg * p0) * (ph2o - p0) * ph2o,
    #         k_AB_ac * K_a_ac,
    #         k_AB_co2 * K_a_co2,
    #         k_AB_IN * K_a_IN,
    #         Vl / Vg,
    #         -kp / (Vg * p0) * (ph2o - p0) * ph2o,
    #     ]
    # )

    c = np.array(
        [
            0.0100000000000000,
            3,
            3.16227766016838e-20,
            8.31508422381744e-14,
            200,
            5.68912300000000,
            129.298250000000,
            0.00170000000000000,
            10000000000.0000,
            10000000000.0000,
            10000000000.0000,
            2000,
            128895.359323054,
            93742.0795076758,
            17044.0144559410,
            -70332.1614249235,
            -25575.3314272449,
            -3072.00828974637,
            1.61622812500000,
            0.587719318181818,
            1000,
            -12889.5359323054,
            -9374.20795076758,
            -1704.40144559410,
            7033.21614249235,
            2557.53314272449,
            250.309598974637,
            -985.781671025363,
            173780.082874937,
            4937.07339753436,
            11.1028665270807,
            10,
            307.200828974637,
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
    uNorm = model.set_variable(var_type="_u", var_name="uNorm", shape=(num_inputs, 1))

    # Time-variant parameters
    theta = np.array([kchF, kchS, kpr, kli, kdec, mu_m_ac, K_S_ac, K_I_nh3, fracChFast])

    # States
    x_norm = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(18)
    ]

    # Aux
    PhiNorm = (
        Tx[12] * x_norm[12]
        + (Tx[3] * x_norm[3] - Tx[15] * x_norm[15]) / 17
        - Tx[14] * x_norm[14] / 44
        - Tx[13] * x_norm[13] / 60
    )

    # equivalent proton concentration
    SHPlusNorm = -0.5 * PhiNorm + 0.5 * sqrt(PhiNorm**2 + c[3])

    # overall inhibition factor
    IacNorm = (
        c[2]
        / (c[2] + SHPlusNorm ** (c[1]))
        * x_norm[3]
        / (x_norm[3] + c[7] / Tx[3])
        * theta[7]
        / (theta[7] + Tx[15] * x_norm[15])
    )

    # Aux
    PhiNorm = (
        Tx[12] * x_norm[12]
        + (Tx[3] * x_norm[3] - Tx[15] * x_norm[15]) / 17
        - Tx[14] * x_norm[14] / 44
        - Tx[13] * x_norm[13] / 60
    )

    # equivalent proton concentration:
    SHPlusNorm = -PhiNorm / 2.0 + 0.5 * np.sqrt(PhiNorm**2 + c[3])

    # overall inhibition factor:
    IacNorm = (
        c[2]
        / (c[2] + SHPlusNorm ** (c[1]))
        * x_norm[3]
        / (x_norm[3] + c[7] / Tx[3])
        * theta[7]
        / (theta[7] + Tx[15] * x_norm[15])
    )

    # Set measurement equations as expressions
    model.set_expression(
        "y_1",
        c[12] * Tx[16] ** 2 / Ty[0] * x_norm[16] ** 2
        + c[13] * Tx[16] * Tx[17] / Ty[0] * x_norm[16] * x_norm[17]
        + c[14] * Tx[17] ** 2 / Ty[0] * x_norm[17] ** 2
        + c[15] * Tx[16] / Ty[0] * x_norm[16]
        + c[16] * Tx[17] / Ty[0] * x_norm[17]
        + c[17] / Ty[0],
    )
    model.set_expression("y_2", c[18] * Tx[16] / Ty[1] * x_norm[16])
    model.set_expression("y_3", c[19] * Tx[17] / Ty[2] * x_norm[17])
    model.set_expression("y_4", -1 / Ty[3] * log10(SHPlusNorm))
    model.set_expression("y_5", Tx[3] / Ty[4] * x_norm[3])
    model.set_expression("y_6", 1 / Ty[5] * (1 - Tx[4] * x_norm[4] / c[20]))
    model.set_expression(
        "y_7", 1 / Ty[6] * (1 - Tx[11] * x_norm[11] / (c[20] - Tx[4] * x_norm[4]))
    )
    model.set_expression("y_8", Tx[0] / Ty[7] * x_norm[0])

    # model.set_expression(
    #     "y_1",
    #     c[12] * Tx[16] ** 2 / Ty[0] * x_norm[16] ** 2
    #     + c[13] * Tx[16] * Tx[17] / Ty[0] * x_norm[16] * x_norm[17]
    #     + c[14] * Tx[17] ** 2 / Ty[0] * x_norm[17] ** 2
    #     + c[15] * Tx[16] / Ty[0] * x_norm[16]
    #     + c[16] * Tx[17] / Ty[0] * x_norm[17]
    #     + c[17] / Ty[0],
    # )  # volFlow
    # model.set_expression("y_2", c[18] * Tx[16] / Ty[1] * x_norm[16])  # pch4
    # model.set_expression("y_3", c[19] * Tx[17] / Ty[2] * x_norm[17])  # pco2
    # model.set_expression("y_4", -1 / Ty[3] * log10(SHPlusNorm))  # pH
    # model.set_expression(
    #     "y_5", Tx[3] / Ty[4] * x_norm[3] - Tx[15] / Ty[4] * x_norm[15]
    # )  # Snh4
    # model.set_expression("y_6", 1 / Ty[5] * (1 - Tx[4] * x_norm[4] / c[20]))  # TS
    # model.set_expression(
    #     "y_7", 1 / Ty[6] * (1 - Tx[11] * x_norm[11] / (c[20] - Tx[4] * x_norm[4]))
    # )  # VS
    # model.set_expression("y_8", Tx[0] / Ty[7] * x_norm[0])

    # if num_inputs > 1:
    #     sum_u = cumsum(u)[1]
    # else:
    #     sum_u = u

    # Differential equations

    model.set_rhs(
        "x_1",
        c[0] * (xi_norm[0] - x_norm[0]) * Tu * uNorm
        + a[0, 0] * theta[0] * Tx[5] / Tx[0] * x_norm[5]
        + a[0, 1] * theta[1] * Tx[6] / Tx[0] * x_norm[6]
        + a[0, 2] * theta[2] * Tx[7] / Tx[0] * x_norm[7]
        + a[0, 3] * theta[3] * Tx[8] / Tx[0] * x_norm[8]
        - a[0, 4]
        * theta[5]
        * Tx[10]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm,
    )
    model.set_rhs(
        "x_2",
        c[0] * (xi_norm[1] - x_norm[1]) * Tu * uNorm
        + a[1, 0] * theta[0] * Tx[5] / Tx[1] * x_norm[5]
        + a[1, 1] * theta[1] * Tx[6] / Tx[1] * x_norm[6]
        + a[1, 2] * theta[2] * Tx[7] / Tx[1] * x_norm[7]
        + a[1, 3] * theta[3] * Tx[8] / Tx[1] * x_norm[8]
        + a[1, 4]
        * theta[5]
        * Tx[0]
        * Tx[10]
        / Tx[1]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm
        - c[4] * x_norm[1]
        + c[5] * Tx[16] / Tx[1] * x_norm[16],
    )
    model.set_rhs(
        "x_3",
        c[0] * (xi_norm[2] - x_norm[2]) * Tu * uNorm
        + a[2, 0] * theta[0] * Tx[5] / Tx[2] * x_norm[5]
        + a[2, 1] * theta[1] * Tx[6] / Tx[2] * x_norm[6]
        + a[2, 2] * theta[2] * Tx[7] / Tx[2] * x_norm[7]
        - a[2, 3] * theta[3] * Tx[8] / Tx[2] * x_norm[8]
        + a[2, 4]
        * theta[5]
        * Tx[0]
        * Tx[10]
        / Tx[2]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm
        - c[4] * x_norm[2]
        + c[4] * Tx[14] / Tx[2] * x_norm[14]
        + c[6] * Tx[17] / Tx[2] * x_norm[17],
    )
    model.set_rhs(
        "x_4",
        c[0] * (xi_norm[3] - x_norm[3]) * Tu * uNorm
        - a[3, 0] * theta[0] * Tx[5] / Tx[3] * x_norm[5]
        - a[3, 1] * theta[1] * Tx[6] / Tx[3] * x_norm[6]
        + a[3, 2] * theta[2] * Tx[7] / Tx[3] * x_norm[7]
        - a[3, 3] * theta[3] * Tx[8] / Tx[3] * x_norm[8]
        - a[3, 4]
        * theta[5]
        * Tx[0]
        * Tx[10]
        / Tx[3]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm,
    )
    model.set_rhs(
        "x_5",
        c[0] * (xi_norm[4] - x_norm[4]) * Tu * uNorm
        - a[4, 0] * theta[0] * Tx[5] / Tx[4] * x_norm[5]
        - a[4, 1] * theta[1] * Tx[6] / Tx[4] * x_norm[6]
        - a[4, 2] * theta[2] * Tx[7] / Tx[4] * x_norm[7]
        - a[4, 3] * theta[3] * Tx[8] / Tx[4] * x_norm[8]
        + a[4, 4]
        * theta[5]
        * Tx[0]
        * Tx[10]
        / Tx[4]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm,
    )
    model.set_rhs(
        "x_6",
        c[0] * (theta[8] * xi_norm[5] - x_norm[5]) * Tu * uNorm
        - theta[0] * x_norm[5]
        + a[5, 5] * theta[4] * Tx[9] / Tx[5] * x_norm[9]
        + a[5, 6] * theta[4] * Tx[10] / Tx[5] * x_norm[10],
    )
    model.set_rhs(
        "x_7",
        c[0] * ((1 - theta[8]) * Tx[5] / Tx[6] * xi_norm[5] - x_norm[6]) * Tu * uNorm
        - theta[1] * x_norm[6],
    )
    model.set_rhs(
        "x_8",
        c[0] * (xi_norm[7] - x_norm[7]) * Tu * uNorm
        - theta[2] * x_norm[7]
        + a[7, 5] * theta[4] * Tx[9] / Tx[7] * x_norm[9]
        + a[7, 6] * theta[4] * Tx[10] / Tx[7] * x_norm[10],
    )
    model.set_rhs(
        "x_9",
        c[0] * (xi_norm[8] - x_norm[8]) * Tu * uNorm
        - theta[3] * x_norm[8]
        + a[8, 5] * theta[4] * Tx[9] / Tx[8] * x_norm[9]
        + a[8, 6] * theta[4] * Tx[10] / Tx[8] * x_norm[10],
    )
    model.set_rhs(
        "x_10",
        c[0] * (xi_norm[9] - x_norm[9]) * Tu * uNorm
        + a[9, 0] * theta[0] * Tx[5] / Tx[9] * x_norm[5]
        + a[9, 1] * theta[1] * Tx[6] / Tx[9] * x_norm[6]
        + a[9, 2] * theta[2] * Tx[7] / Tx[9] * x_norm[7]
        + a[9, 3] * theta[3] * Tx[8] / Tx[9] * x_norm[8]
        - theta[4] * x_norm[9],
    )
    model.set_rhs(
        "x_11",
        c[0] * (xi_norm[10] - x_norm[10]) * Tu * uNorm
        + theta[5]
        * Tx[0]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm
        - theta[4] * x_norm[10],
    )
    model.set_rhs("x_12", c[0] * (xi_norm[11] - x_norm[11]) * Tu * uNorm)
    model.set_rhs("x_13", c[0] * (xi_norm[12] - x_norm[12]) * Tu * uNorm)
    model.set_rhs(
        "x_14",
        c[28] * (Tx[0] / Tx[13] * x_norm[0] - x_norm[13])
        - c[8] * x_norm[13] * SHPlusNorm,
    )
    model.set_rhs(
        "x_15",
        c[29] * (Tx[2] / Tx[14] * x_norm[2] - x_norm[14])
        - c[9] * x_norm[14] * SHPlusNorm,
    )
    model.set_rhs(
        "x_16",
        c[30] * (Tx[3] / Tx[15] * x_norm[3] - x_norm[15])
        - c[10] * x_norm[15] * SHPlusNorm,
    )
    model.set_rhs(
        "x_17",
        c[21] * Tx[16] ** 2 * x_norm[16] ** 3
        + c[22] * Tx[16] * Tx[17] * x_norm[16] ** 2 * x_norm[17]
        + c[23] * Tx[17] ** 2 * x_norm[16] * x_norm[17] ** 2
        + c[24] * Tx[16] * x_norm[16] ** 2
        + c[25] * Tx[17] * x_norm[16] * x_norm[17]
        + c[11] * Tx[1] / Tx[16] * x_norm[1]
        + c[26] * x_norm[16],
    )
    model.set_rhs(
        "x_18",
        c[23] * Tx[17] ** 2 * x_norm[17] ** 3
        + c[22] * Tx[16] * Tx[17] * x_norm[16] * x_norm[17] ** 2
        + c[21] * Tx[16] ** 2 * x_norm[16] ** 2 * x_norm[17]
        + c[25] * Tx[17] * x_norm[17] ** 2
        + c[24] * Tx[16] * x_norm[16] * x_norm[17]
        + c[11] * Tx[2] / Tx[17] * x_norm[2]
        - c[11] * Tx[14] / Tx[17] * x_norm[14]
        + c[27] * x_norm[17],
    )

    # model.set_rhs(
    #     "x_1",
    #     c[0] * (xi_norm[0] - x_norm[0]) * Tu * uNorm
    #     + a[0, 0] * theta[0] * Tx[5] / Tx[0] * x_norm[5]
    #     + a[0, 1] * theta[1] * Tx[6] / Tx[0] * x_norm[6]
    #     + a[0, 2] * theta[2] * Tx[7] / Tx[0] * x_norm[7]
    #     + a[0, 3] * theta[3] * Tx[8] / Tx[0] * x_norm[8]
    #     - a[0, 4]
    #     * theta[5]
    #     * Tx[10]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm,
    # )
    # model.set_rhs(
    #     "x_2",
    #     c[0] * (xi_norm[1] - x_norm[1]) * Tu * uNorm
    #     + a[1, 0] * theta[0] * Tx[5] / Tx[1] * x_norm[5]
    #     + a[1, 1] * theta[1] * Tx[6] / Tx[1] * x_norm[6]
    #     + a[1, 2] * theta[2] * Tx[7] / Tx[1] * x_norm[7]
    #     + a[1, 3] * theta[3] * Tx[8] / Tx[1] * x_norm[8]
    #     + a[1, 4]
    #     * theta[5]
    #     * Tx[0]
    #     * Tx[10]
    #     / Tx[1]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm
    #     - c[4] * x_norm[1]
    #     + c[5] * Tx[16] / Tx[1] * x_norm[16],
    # )
    # model.set_rhs(
    #     "x_3",
    #     c[0] * (xi_norm[2] - x_norm[2]) * Tu * uNorm
    #     + a[2, 0] * theta[0] * Tx[5] / Tx[2] * x_norm[5]
    #     + a[2, 1] * theta[1] * Tx[6] / Tx[2] * x_norm[6]
    #     + a[2, 2] * theta[2] * Tx[7] / Tx[2] * x_norm[7]
    #     - a[2, 3] * theta[3] * Tx[8] / Tx[2] * x_norm[8]
    #     + a[2, 4]
    #     * theta[5]
    #     * Tx[0]
    #     * Tx[10]
    #     / Tx[2]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm
    #     - c[4] * x_norm[2]
    #     + c[4] * Tx[14] / Tx[2] * x_norm[14]
    #     + c[6] * Tx[17] / Tx[2] * x_norm[17],
    # )
    # model.set_rhs(
    #     "x_4",
    #     c[0] * (xi_norm[3] - x_norm[3]) * Tu * uNorm
    #     - a[3, 0] * theta[0] * Tx[5] / Tx[3] * x_norm[5]
    #     - a[3, 1] * theta[1] * Tx[6] / Tx[3] * x_norm[6]
    #     + a[3, 2] * theta[2] * Tx[7] / Tx[3] * x_norm[7]
    #     - a[3, 3] * theta[3] * Tx[8] / Tx[3] * x_norm[8]
    #     - a[3, 4]
    #     * theta[5]
    #     * Tx[0]
    #     * Tx[10]
    #     / Tx[3]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm,
    # )
    # model.set_rhs(
    #     "x_5",
    #     c[0] * (xi_norm[4] - x_norm[4]) * Tu * uNorm
    #     - a[4, 0] * theta[0] * Tx[5] / Tx[4] * x_norm[5]
    #     - a[4, 1] * theta[1] * Tx[6] / Tx[4] * x_norm[6]
    #     - a[4, 2] * theta[2] * Tx[7] / Tx[4] * x_norm[7]
    #     - a[4, 3] * theta[3] * Tx[8] / Tx[4] * x_norm[8]
    #     + a[4, 4]
    #     * theta[5]
    #     * Tx[0]
    #     * Tx[10]
    #     / Tx[4]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm,
    # )
    # model.set_rhs(
    #     "x_6",
    #     c[0] * (theta[8] * xi_norm[5] - x_norm[5]) * Tu * uNorm
    #     - theta[0] * x_norm[5]
    #     + a[5, 5] * theta[4] * Tx[9] / Tx[5] * x_norm[9]
    #     + a[5, 6] * theta[4] * Tx[10] / Tx[5] * x_norm[10],
    # )
    # model.set_rhs(
    #     "x_7",
    #     c[0] * ((0 - theta[8]) * Tx[5] / Tx[6] * xi_norm[5] - x_norm[6]) * Tu * uNorm
    #     - theta[1] * x_norm[6],
    # )
    # model.set_rhs(
    #     "x_8",
    #     c[0] * (xi_norm[7] - x_norm[7]) * Tu * uNorm
    #     - theta[2] * x_norm[7]
    #     + a[7, 5] * theta[4] * Tx[9] / Tx[7] * x_norm[9]
    #     + a[7, 6] * theta[4] * Tx[10] / Tx[7] * x_norm[10],
    # )
    # model.set_rhs(
    #     "x_9",
    #     c[0] * (xi_norm[8] - x_norm[8]) * Tu * uNorm
    #     - theta[3] * x_norm[8]
    #     + a[8, 5] * theta[4] * Tx[9] / Tx[8] * x_norm[9]
    #     + a[8, 6] * theta[4] * Tx[10] / Tx[8] * x_norm[10],
    # )
    # model.set_rhs(
    #     "x_10",
    #     c[0] * (xi_norm[9] - x_norm[9]) * Tu * uNorm
    #     + a[9, 0] * theta[0] * Tx[5] / Tx[9] * x_norm[5]
    #     + a[9, 1] * theta[1] * Tx[6] / Tx[9] * x_norm[6]
    #     + a[9, 2] * theta[2] * Tx[7] / Tx[9] * x_norm[7]
    #     + a[9, 3] * theta[3] * Tx[8] / Tx[9] * x_norm[8]
    #     - theta[4] * x_norm[9],
    # )
    # model.set_rhs(
    #     "x_11",
    #     c[0] * (xi_norm[10] - x_norm[10]) * Tu * uNorm
    #     + theta[5]
    #     * Tx[0]
    #     * x_norm[0]
    #     * x_norm[10]
    #     / (theta[6] + Tx[0] * x_norm[0])
    #     * IacNorm
    #     - theta[4] * x_norm[10],
    # )
    # model.set_rhs("x_12", c[0] * (xi_norm[11] - x_norm[11]) * Tu * uNorm)
    # model.set_rhs("x_13", c[0] * (xi_norm[12] - x_norm[12]) * Tu * uNorm)
    # model.set_rhs(
    #     "x_14",
    #     c[28] * (Tx[0] / Tx[13] * x_norm[0] - x_norm[13])
    #     - c[8] * x_norm[13] * SHPlusNorm,
    # )
    # model.set_rhs(
    #     "x_15",
    #     c[29] * (Tx[2] / Tx[14] * x_norm[2] - x_norm[14])
    #     - c[9] * x_norm[14] * SHPlusNorm,
    # )
    # model.set_rhs(
    #     "x_16",
    #     c[30] * (Tx[3] / Tx[15] * x_norm[3] - x_norm[15])
    #     - c[10] * x_norm[15] * SHPlusNorm,
    # )
    # model.set_rhs(
    #     "x_17",
    #     c[21] * Tx[16] ** 2 * x_norm[16] ** 3
    #     + c[22] * Tx[16] * Tx[17] * x_norm[16] ** 2 * x_norm[17]
    #     + c[23] * Tx[17] ** 2 * x_norm[16] * x_norm[17] ** 2
    #     + c[24] * Tx[16] * x_norm[16] ** 2
    #     + c[25] * Tx[17] * x_norm[16] * x_norm[17]
    #     + c[11] * Tx[1] / Tx[16] * x_norm[1]
    #     + c[26] * x_norm[16],
    # )
    # model.set_rhs(
    #     "x_18",
    #     c[23] * Tx[17] ** 2 * x_norm[17] ** 3
    #     + c[22] * Tx[16] * Tx[17] * x_norm[16] * x_norm[17] ** 2
    #     + c[21] * Tx[16] ** 2 * x_norm[16] ** 2 * x_norm[17]
    #     + c[25] * Tx[17] * x_norm[17] ** 2
    #     + c[24] * Tx[16] * x_norm[16] * x_norm[17]
    #     + c[11] * Tx[2] / Tx[17] * x_norm[2]
    #     - c[11] * Tx[14] / Tx[17] * x_norm[14]
    #     + c[27] * x_norm[17],
    # )

    # Build the model
    model.setup()

    return model
