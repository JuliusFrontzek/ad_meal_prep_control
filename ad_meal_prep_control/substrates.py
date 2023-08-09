import numpy as np
import casadi
from do_mpc.model import Model
import substrate_uncertainties


class ADM1_R3_FRAC_XI:
    """
    xi values in the following order:
        S_ac, S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ac, X_ash, S_ion, S_ac-, S_hco3-, S_nh3, S_ch4_g, S_co2_g

    unit: [g/l]
    """

    corn = np.array(
        [
            2.14500583382249,
            0.0,
            0.0,
            0.592294500000000,
            960.511750000000,
            21.2528071170000,
            0,
            4.74931351500000,
            1.38096405000000,
            0.0,
            0.0,
            14,
            0.0487500000000000,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )


class ADM1_R3_FRAC_NORM_XI:
    """
    xi values in the following order:
        S_ac, S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ac, X_ash, S_ion, S_ac-, S_hco3-, S_nh3, S_ch4_g, S_co2_g

    unit: [1]
    """

    corn = np.array(
        [
            2.145005833822492,
            0,
            0,
            0.592294500000000,
            960.5117500000000,
            21.252807116999996,
            0,
            4.749313514999999,
            1.380964050000000,
            0,
            0,
            14,
            0.048750000000000,
            0,
            0,
            0,
            0,
            0,
        ]
    )


class ADM1_R4_FRAC_XI:
    """
    xi values in the following order:
        S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ash, S_ch4_g, S_co2_g

    unit: [g/l]
    """

    corn = np.array([0, 0, 0.592, 960.512, 23.398, 0, 4.75, 1.381, 0, 17, 0, 0])
    manure = np.array([0, 0, 1.710, 960.512, 16.1, 0, 10.8, 1.4, 0, 17, 0, 0])
    shit = np.array([0, 0, 2.510, 960.512, 10.1, 0, 30.8, 1.4, 0, 17, 0, 0])

    corn[3] = corn[3] / 1000.0
    corn[9] = corn[9] / 1000.0

    manure[3] = manure[3] / 1000.0
    manure[9] = manure[9] / 1000.0

    shit[3] = shit[3] / 1000.0
    shit[9] = shit[9] / 1000.0


def xi_values(
    model: Model,
    substrate_names: list[str],
    model_type: str = "R3-frac-norm",
    include_uncertainties: bool = True,
) -> np.ndarray:
    assert isinstance(substrate_names, list), f"'substrate_names' must be of type list."
    if model_type == "R3-frac":
        xi_type = ADM1_R3_FRAC_XI
    elif model_type == "R4-frac":
        xi_type = ADM1_R4_FRAC_XI
    elif model_type == "R3-frac-norm":
        xi_type = ADM1_R3_FRAC_NORM_XI
    else:
        raise NotImplementedError(f"Model '{model_type}' not implemented.")

    xi = []

    for substrate in substrate_names:
        try:
            xi.append(getattr(xi_type, substrate).T)
        except AttributeError:
            raise AttributeError(
                f"The specified model '{model_type}' does not have xi values specified for the substrate '{substrate}'."
            )

    xi = casadi.horzcat(*xi)

    return xi
