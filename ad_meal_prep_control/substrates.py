from __future__ import annotations
import numpy as np
from dataclasses import dataclass
import substrate_uncertainties
from uncertainties import ufloat
from copy import deepcopy
from typing import Union


@dataclass
class Substrate:
    """
    Attributes:
        nominal_values:
                                Case 1:
                                    Nominal BMP and TS only -> nominal values for XP, XL, XA will be
                                    computed from the supplied xi values as well as the nominal BMP and TS.
                                Case 2:
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
        state:
                                The physical state of the substrate.
                                Must be either solid or liquid.

    Order: XP, XL, XA, BMP, TS
    """

    nominal_values: np.ndarray
    variation_coefficients: np.ndarray
    xi: list
    state: str

    def __post_init__(self):
        assert self.state in [
            "solid",
            "liquid",
        ], f"The physical state of an input must be either solid or liquid, not '{self.state}'."

        # Conversion of percentages to decimal numbers where appropriate
        if self.nominal_values.shape[0] == 5:
            self.nominal_values /= 100.0  # unit change [%] -> [-]
            self.nominal_values[3] *= 100.0  # undo unit change for BMP

        # Computation of remaining nominal values
        else:
            rho_fm = 1000.0
            bmp_nominal = self.nominal_values[0]
            ts_nominal = self.nominal_values[1] / 100.0  # unit change [%] -> [-]
            self.nominal_values = np.empty(5)

            # BMP
            self.nominal_values[3] = bmp_nominal

            # TS
            self.nominal_values[4] = ts_nominal

            # XP
            self.nominal_values[0] = self.xi[7] / (ts_nominal * rho_fm)

            # XL
            self.nominal_values[1] = self.xi[8] / (ts_nominal * rho_fm)

            # XA
            self.nominal_values[2] = self.xi[11] / (ts_nominal * rho_fm)

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
    state="solid",
)

CORN_SILAGE = Substrate(
    nominal_values=np.array([357.0, 31.6818754946664]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 9, 2.14901476466728]
    ),
    xi=[
        0.479732628205128,
        0,
        0,
        0.304492136205810,
        683.181245053336,
        221.044354457402,
        0,
        24.6510280564042,
        11.5108145186325,
        0.287465749624490,
        0.0151297762960258,
        14.2232290261481,
        -0.00992965566117833,
        0.478861243573483,
        0,
        0.0127401924325350,
        0,
        0,
    ],
    state="solid",
)

GRASS_SILAGE = Substrate(
    nominal_values=np.array([315.0, 39.3633764654239]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 6, 2.14901476466728]
    ),
    xi=[
        0.592421666666667,
        0,
        0,
        0.532977220532501,
        606.366235345761,
        205.314797104047,
        0,
        42.2696292679289,
        15.4197008294951,
        0.333138561121863,
        0.0175336084800981,
        42.9615950522776,
        -0.0214952155643065,
        0.591345594068231,
        0,
        0.0223001895429972,
        0,
        0,
    ],
    state="solid",
)

CROP_STRAW = Substrate(
    nominal_values=np.array([240.0, 34.6393249615227]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 5, 2.14901476466728]
    ),
    xi=[
        0.116910000000000,
        0,
        0,
        0.185686366820576,
        653.606750384773,
        174.363204087472,
        0,
        10.9013441888264,
        6.27023373775149,
        0.318426575098357,
        0.0167592934262293,
        11.2073810906402,
        -0.00897714096612082,
        0.116697645093754,
        0,
        0.00776926483182942,
        0,
        0,
    ],
    state="solid",
)

CATTLE_MANURE = Substrate(
    nominal_values=np.array([230.0, 8.44120640954394]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 7, 2.14901476466728]
    ),
    xi=[
        2.61861507142857,
        0,
        0,
        1.25605273489858,
        915.587935904561,
        18.1866451576761,
        0,
        13.7032984788022,
        3.21560794289618,
        0.0609005003485669,
        0.00320528949202984,
        20.3062742548427,
        -0.0303205187788519,
        2.61385862836992,
        0,
        0.0525542424425818,
        0,
        0,
    ],
    state="liquid",
)

SWINE_MANURE = Substrate(
    nominal_values=np.array([230.0, 5.98037364701355]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 4, 2.14901476466728]
    ),
    xi=[
        0.445798333333333,
        0,
        0,
        2.04031388795041,
        940.196263529865,
        9.00018594920682,
        0,
        11.4319940663316,
        3.06378373285972,
        0.0407603892852647,
        0.00214528364659288,
        16.8980635382779,
        -0.112601361807505,
        0.444988586833636,
        0,
        0.0853683509832653,
        0,
        0,
    ],
    state="liquid",
)

CHICKEN_DRY_MANURE = Substrate(
    nominal_values=np.array([272.0, 52.5728239460278]),
    variation_coefficients=np.array(
        [5.22776012611444, 12.8464621768132, 17.4321599553856, 6, 2.14901476466728]
    ),
    xi=[
        0.716442857142857,
        0,
        0,
        2.22410843045102,
        474.271760539722,
        78.3900890457954,
        0,
        134.976604663233,
        23.6609869497833,
        0.347698693319360,
        0.0182999312273347,
        159.729614913584,
        -0.118910256716285,
        0.715141512897205,
        0,
        0.0930584603853837,
        0,
        0,
    ],
    state="solid",
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
