import numpy as np
import casadi


class ADM1_R3_FRAC_XI:
    """
    xi values in the following order:
        S_ac, S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ac, X_ash, S_ion, S_ac-, S_hco3-, S_nh3, S_ch4_g, S_co2_g

    unit: g/l
    """

    corn = np.array(
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
    )


class ADM1_R4_FRAC_XI:
    """
    xi values in the following order:
        S_ch4, S_IC, S_IN, S_h2o, X_chF, X_chS, X_pr, X_li, X_bac, X_ash, S_ch4_g, S_co2_g

    unit: g/l
    """

    corn = np.array([0, 0, 0.592, 960.512, 23.398, 0, 4.75, 1.381, 0, 17, 0, 0])
    manure = np.array([0, 0, 1.710, 960.512, 16.1, 0, 10.8, 1.4, 0, 17, 0, 0])


def xi_values(model_type: str, substrate_names: list[str]) -> np.ndarray:
    assert isinstance(substrate_names, list), f"'substrate_names' must be of type list."
    if model_type == "R3-frac":
        model = ADM1_R3_FRAC_XI
    elif model_type == "R4-frac":
        model = ADM1_R4_FRAC_XI
    else:
        raise NotImplementedError(f"Model '{model_type}' not implemented.")

    xi = []

    for substrate in substrate_names:
        try:
            xi.append(getattr(model, substrate).T)
        except AttributeError:
            raise AttributeError(
                f"The specified model '{model_type}' does not have xi values specified for the substrate '{substrate}'."
            )

    xi = casadi.horzcat(*xi)

    return xi