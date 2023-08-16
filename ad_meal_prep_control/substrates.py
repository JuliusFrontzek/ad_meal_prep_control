from __future__ import annotations
import numpy as np
import casadi
from do_mpc.model import Model
from dataclasses import dataclass
from typing import ClassVar
import substrate_uncertainties
from uncertainties import ufloat
from copy import deepcopy


@dataclass
class Substrate:
    """
    Attributes:
        nominal_values:
                                Nominal values of XP, XL, XA, BMP, TS in this order
                                Note: The nominal values for XP, XL, XA and TS are to
                                be supplied as percentages.
        variation_coefficients:
                                Variation coefficients of XP, XL, XA, BMP, TS in this order.
                                Note: The variation coefficients are all to be supplied as percentages.
        xi:
                                Xi values describing the composition of the substrate.
                                The uncertain ones (Carbohydrates, proteins, lipids) will be
                                ignored later when used with do-mpc.

    Order: XP, XL, XA, BMP, TS
    """

    nominal_values: np.ndarray
    variation_coefficients: np.ndarray
    xi: list

    def __post_init__(self):
        # Conversion of percentages to decimal numbers where appropriate
        self.nominal_values /= 100.0
        self.nominal_values[3] *= 100.0
        self.variation_coefficients /= 100.0

        self._set_std_devs()

    def _set_std_devs(self) -> np.ndarray:
        self.std_devs = self.variation_coefficients * self.nominal_values

    def get_uncertain_xi_ch_pr_li(self):
        """
        Returns the uncertain values for xi_ch, xi_pr and xi_li.
        """

        # Create ufloat values for basic components of computation
        xp = ufloat(self.nominal_values[0], self.std_devs[0])
        xl = ufloat(self.nominal_values[1], self.std_devs[1])
        xa = ufloat(self.nominal_values[2], self.std_devs[2])
        bmp = ufloat(self.nominal_values[3], self.std_devs[3])
        ts = ufloat(self.nominal_values[4], self.std_devs[4])
        xc = substrate_uncertainties.xc_(xa, xp, xl)

        fq_ges = substrate_uncertainties.f_q_qes_(bmp)
        fq_ch = substrate_uncertainties.fq_ch_(xc, fq_ges, xa, xp, xl)

        xi_ch = substrate_uncertainties.x_ch_in_(fq_ch, xc, ts, 1000.0)
        xi_pr = substrate_uncertainties.x_pr_in_(1.0, xp, ts, 1000.0)
        xi_li = substrate_uncertainties.x_li_in_(1.0, xl, ts, 1000.0)

        return xi_ch, xi_pr, xi_li

    def create_similar_substrate(self) -> Substrate:
        # Create new substrate
        sub = deepcopy(self)

        # Modify its values slightly
        sub.nominal_values *= np.random.normal(
            loc=1.0, scale=0.1, size=sub.nominal_values.shape[0]
        )
        sub.variation_coefficients *= np.random.normal(
            loc=1.0, scale=0.1, size=sub.variation_coefficients.shape[0]
        )
        xi = np.array(sub.xi)
        xi *= np.random.normal(loc=1.0, scale=0.1, size=xi.shape[0])
        sub.xi = xi.tolist()
        sub._set_std_devs()

        return sub


CORN = Substrate(
    nominal_values=np.array([7.25, 3.27, 4.4, 357.0, 3.518]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 9.0, 2.14901476466728]
    ),
    xi=[
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
    ],
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


def xi_values(
    substrate_names: list[str], model_type: str = "R3-frac-norm"
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

    xi = list(np.array(xi).T)

    x_ch_nom = xi[5]
    x_pr_nom = xi[7]
    x_li_nom = xi[8]

    return xi, x_ch_nom, x_pr_nom, x_li_nom


def get_subs(substrate_names: list[str]):
    subs = []
    for substrate in substrate_names:
        try:
            subs.append(globals()[substrate.upper()])
        except:
            raise

    return subs
