from __future__ import annotations
import numpy as np
from dataclasses import dataclass
import substrate_uncertainties
from uncertainties import ufloat
from copy import deepcopy


@dataclass
class Substrate:
    """
    Attributes:
        name:
                                Name of the substrate.
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
        cost:
                                Cost of the substrate measured as €/t fresh matter.
        amount_remaining:
                                The amount of the substrate remaining in kilograms.
        days_remaining:
                                The amount of days remaining to feed that substrate.
        u_target:
                                The amount of substrate that should be fed on average per time step in kilograms.

    Order: XP, XL, XA, BMP, TS
    """

    name: str
    nominal_values: np.ndarray
    variation_coefficients: np.ndarray
    xi: list
    state: str
    cost: float

    def __post_init__(self):
        assert self.state in [
            "solid",
            "liquid",
        ], f"The physical state of an input must be either solid or liquid, not '{self.state}'."

        assert self.cost >= 0.0, f"Substrate cost must be larger than or equal to 0"
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

        self.limited = False

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

    def set_limit(self, amount_remaining, days_remaining, t_step):
        """
        Limit the amount of substrate that can be fed and add a time limit to it.

        Parameters:
            amount_remaining:   The amount of the substrate remaining in kilograms.
            days_remaining:     The amount of days remaining to feed that substrate.
            t_step:             The time step size in days.
        """
        self.amount_remaining = amount_remaining
        self.days_remaining = days_remaining

        n_steps_total = self.days_remaining / t_step

        self.u_target = self.amount_remaining / n_steps_total
        self.limited = True


STANDARD_SUBSTRATE = Substrate(
    name="STANDARD_SUBSTRATE",
    nominal_values=np.array(
        [
            7.25,
            3.27,
            4.4,
            357.0,
            3.518,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.22776012611444,
            12.8464621768132,
            17.4321599553856,
            9.0,
            2.14901476466728,
        ]
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
    cost=10.0,
)

CORN_SILAGE = Substrate(
    name="CORN_SILAGE",
    nominal_values=np.array(
        [
            357.0,
            33.7286380429205,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.22776012611444,
            12.8464621768132,
            17.4321599553856,
            9.0,
            2.14901476466728,
        ]
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
    cost=40.0,
)

GRASS_SILAGE = Substrate(
    name="GRASS_SILAGE",
    nominal_values=np.array(
        [
            315.0,
            31.7409815648268,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.22776012611444,
            12.8464621768132,
            17.4321599553856,
            6.0,
            2.14901476466728,
        ]
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
    cost=40.0,
)


CATTLE_MANURE = Substrate(
    name="CATTLE_MANURE",
    nominal_values=np.array(
        [
            230.0,
            8.0840310707442,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.22776012611444,
            12.8464621768132,
            17.4321599553856,
            7.0,
            2.14901476466728,
        ]
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
    cost=1.0,
)

SUGAR_BEET_SILAGE = Substrate(
    name="SUGAR_BEET_SILAGE",
    nominal_values=np.array(
        [
            657.0,
            31.8334469458509,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.22776012611444,
            12.8464621768132,
            17.4321599553856,
            3.0,
            2.14901476466728,
        ]
    ),
    xi=[
        0.386290000000000,
        0,
        0,
        0.388448549267835,
        681.665530541491,
        230.960309794950,
        0,
        13.4892573835735,
        0.968888766166786,
        0.280578693186050,
        0.0147672996413710,
        22.9884766310881,
        -0.0154667579678940,
        0.385588344224328,
        0,
        0.0162529953301199,
        0,
        0,
    ],
    state="solid",
    cost=50.0,
)
