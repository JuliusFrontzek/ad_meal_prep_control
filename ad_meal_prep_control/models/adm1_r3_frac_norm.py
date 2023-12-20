import numpy as np
from casadi import *
from casadi.tools import *
import sys
import os
from ad_meal_prep_control.params_R3 import *
from ad_meal_prep_control.antoine_water import vapour_pressure_h2o

rel_do_mpc_path = os.path.join("..", "..")
sys.path.append(rel_do_mpc_path)
import do_mpc


def adm1_r3_frac_norm(
    xi_norm: list,
    Tu: np.ndarray,
    Tx: np.ndarray,
    Ty: np.ndarray,
    external_gas_storage_model: bool,
    num_dictated_subs: int,
    limited_subs_indices: list[int],
):
    """
    Normalized ADM1-R3-frac model.
    """
    model_type = "continuous"
    model = do_mpc.model.Model(model_type, "SX")

    num_inputs = xi_norm[0].shape[0] - num_dictated_subs

    xi_norm[5] = model.set_variable(
        var_type="_p", var_name="xi_ch_norm", shape=(len(xi_norm[0]), 1)
    )
    xi_norm[7] = model.set_variable(
        var_type="_p", var_name="xi_pr_norm", shape=(len(xi_norm[0]), 1)
    )
    xi_norm[8] = model.set_variable(
        var_type="_p", var_name="xi_li_norm", shape=(len(xi_norm[0]), 1)
    )

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
    u_norm = model.set_variable(var_type="_u", var_name="u_norm", shape=(num_inputs, 1))

    if num_dictated_subs:
        dictated_sub_feed = model.set_variable(
            var_type="_tvp", var_name="dictated_sub_feed", shape=(num_dictated_subs, 1)
        )
        u_norm = vertcat(u_norm, dictated_sub_feed)

    # Time-variant parameters (time-invariant for now -> making them invariant and estimating them: TODO for Simon)
    theta = np.array([kchF, kchS, kpr, kli, kdec, mu_m_ac, K_S_ac, K_I_nh3, fracChFast])

    if external_gas_storage_model:
        num_states = 20
    else:
        num_states = 18

    # States
    x_norm = [
        model.set_variable(var_type="_x", var_name=f"x_{i+1}", shape=(1, 1))
        for i in range(num_states)
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
    y_norm = []
    y_norm.append(
        model.set_expression(
            "y_1_norm",
            c[12] * Tx[16] ** 2 / Ty[0] * x_norm[16] ** 2
            + c[13] * Tx[16] * Tx[17] / Ty[0] * x_norm[16] * x_norm[17]
            + c[14] * Tx[17] ** 2 / Ty[0] * x_norm[17] ** 2
            + c[15] * Tx[16] / Ty[0] * x_norm[16]
            + c[16] * Tx[17] / Ty[0] * x_norm[17]
            + c[17] / Ty[0],
        )
    )
    y_norm.append(model.set_expression("y_2_norm", c[18] * Tx[16] / Ty[1] * x_norm[16]))
    y_norm.append(model.set_expression("y_3_norm", c[19] * Tx[17] / Ty[2] * x_norm[17]))
    y_norm.append(model.set_expression("y_4_norm", -1 / Ty[3] * log10(SHPlusNorm)))
    y_norm.append(model.set_expression("y_5_norm", Tx[3] / Ty[4] * x_norm[3]))
    y_norm.append(
        model.set_expression("y_6_norm", 1 / Ty[5] * (1 - Tx[4] * x_norm[4] / c[20]))
    )
    y_norm.append(
        model.set_expression(
            "y_7_norm",
            1 / Ty[6] * (1 - Tx[11] * x_norm[11] / (c[20] - Tx[4] * x_norm[4])),
        )
    )
    y_norm.append(model.set_expression("y_8_norm", Tx[0] / Ty[7] * x_norm[0]))

    y = []
    for i in range(8):
        y.append(model.set_expression(f"y_{i+1}", Ty[i] * y_norm[i]))

    model.set_variable(var_type="_tvp", var_name="v_ch4_dot_tank_in_setpoint")
    model.set_variable(var_type="_tvp", var_name="dummy_tvp")

    p_h2o = vapour_pressure_h2o(T)
    p_gas_total_fermenter = model.set_expression(
        "p_gas_total_fermenter", y_norm[1] * Ty[1] + y_norm[2] * Ty[2] + p_h2o
    )
    v_total_dot_tank_in = model.set_expression(
        "v_total_dot_tank_in",
        y_norm[0] * Ty[0] * p_gas_total_fermenter / p_gas_storage * T_gas_storage / T,
    )
    v_ch4_dot_tank_in = model.set_expression(
        f"v_ch4_dot_tank_in",
        v_total_dot_tank_in * y_norm[1] * Ty[1] / p_gas_total_fermenter,
    )

    if external_gas_storage_model:
        v_ch4_dot_tank_out = model.set_variable(
            var_type="_tvp", var_name="v_ch4_dot_tank_out"
        )
        v_h2o = model.set_expression(
            "V_H2O",
            1.0
            / (1.0 - p_h2o / p_gas_storage)
            * p_h2o
            / p_gas_storage
            * (Tx[18] * x_norm[18] + Tx[19] * x_norm[19]),
        )

        v_gas_storage = model.set_expression(
            "v_gas_storage", Tx[18] * x_norm[18] + Tx[19] * x_norm[19] + v_h2o
        )

        y_h2o = model.set_expression(
            "y_h2o",
            v_h2o / (v_gas_storage),
        )
        y_co2 = model.set_expression(
            "y_co2",
            Tx[19] * x_norm[19] / (v_gas_storage),
        )

    if num_inputs > 1:
        sum_u = cumsum((SX(Tu).T * u_norm.T).T)[len(xi_norm[0]) - 1]
    else:
        sum_u = Tu * u_norm

    # Measurements
    u_meas = model.set_meas("u_meas", Tu * u_norm, meas_noise=False)
    y_meas = []
    for idx, y in enumerate(y):
        y_meas.append(model.set_meas(f"y_meas_{idx+1}", y, meas_noise=True))

    # Differential equations
    model.set_rhs(
        "x_1",
        c[0] * (SX(xi_norm[0]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[0])
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
        c[0] * (SX(xi_norm[1]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[1])
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
        c[0] * (SX(xi_norm[2]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[2])
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
        c[0] * (SX(xi_norm[3]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[3])
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
        c[0] * (SX(xi_norm[4]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[4])
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
        c[0]
        * (theta[8] * SX(xi_norm[5]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[5])
        - theta[0] * x_norm[5]
        + a[5, 5] * theta[4] * Tx[9] / Tx[5] * x_norm[9]
        + a[5, 6] * theta[4] * Tx[10] / Tx[5] * x_norm[10],
    )
    model.set_rhs(
        "x_7",
        c[0]
        * (
            (1 - theta[8]) * Tx[5] / Tx[6] * SX(xi_norm[5]).T @ (SX(Tu).T * u_norm.T).T
            - sum_u * x_norm[6]
        )
        - theta[1] * x_norm[6],
    )
    model.set_rhs(
        "x_8",
        c[0] * (SX(xi_norm[7]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[7])
        - theta[2] * x_norm[7]
        + a[7, 5] * theta[4] * Tx[9] / Tx[7] * x_norm[9]
        + a[7, 6] * theta[4] * Tx[10] / Tx[7] * x_norm[10],
    )
    model.set_rhs(
        "x_9",
        c[0] * (SX(xi_norm[8]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[8])
        - theta[3] * x_norm[8]
        + a[8, 5] * theta[4] * Tx[9] / Tx[8] * x_norm[9]
        + a[8, 6] * theta[4] * Tx[10] / Tx[8] * x_norm[10],
    )
    model.set_rhs(
        "x_10",
        c[0] * (SX(xi_norm[9]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[9])
        + a[9, 0] * theta[0] * Tx[5] / Tx[9] * x_norm[5]
        + a[9, 1] * theta[1] * Tx[6] / Tx[9] * x_norm[6]
        + a[9, 2] * theta[2] * Tx[7] / Tx[9] * x_norm[7]
        + a[9, 3] * theta[3] * Tx[8] / Tx[9] * x_norm[8]
        - theta[4] * x_norm[9],
    )
    model.set_rhs(
        "x_11",
        c[0] * (SX(xi_norm[10]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[10])
        + theta[5]
        * Tx[0]
        * x_norm[0]
        * x_norm[10]
        / (theta[6] + Tx[0] * x_norm[0])
        * IacNorm
        - theta[4] * x_norm[10],
    )
    model.set_rhs(
        "x_12",
        c[0] * (SX(xi_norm[11]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[11]),
    )
    model.set_rhs(
        "x_13",
        c[0] * (SX(xi_norm[12]).T @ (SX(Tu).T * u_norm.T).T - sum_u * x_norm[12]),
    )
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

    if external_gas_storage_model:
        model.set_rhs(
            "x_19",
            (v_ch4_dot_tank_in - v_ch4_dot_tank_out) / Tx[18],
        )  # V_CH4
        model.set_rhs(
            "x_20",
            (
                v_total_dot_tank_in * y_norm[2] * Ty[2] / p_gas_total_fermenter
                - y_co2
                / (1.0 - y_co2)
                / (1.0 - y_co2 * y_h2o / ((1.0 - y_co2) * (1.0 - y_h2o)))
                * (v_ch4_dot_tank_out * (1.0 + y_h2o / (1.0 - y_h2o)))
            )
            / Tx[18],
        )  # V_CO2

    # Build the model
    model.setup()

    return model
