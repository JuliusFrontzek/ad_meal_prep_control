from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from ad_meal_prep_control import substrate_uncertainties
from uncertainties import ufloat
from copy import deepcopy
import pandas as pd


@dataclass(kw_only=True)
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

    xi_names = [
        "Sac",
        "Sch4",
        "SIC",
        "SIN",
        "Sh2o",
        "Xch,f",
        "Xch,s",
        "Xpr",
        "Xli",
        "Xbac",
        "Xac",
        "Xash",
        "Sion",
        "Sac-",
        "Shco3-",
        "Snh3",
        "Sch4,gas",
        "Sco2,gas",
    ]

    def __post_init__(self):
        assert self.state in [
            "solid",
            "liquid",
        ], f"The physical state of an input must be either solid or liquid, not '{self.state}'."

        assert self.cost >= 0.0, f"Substrate cost must be larger than or equal to 0"
        # Conversion of percentages to decimal numbers where appropriate
        if self.nominal_values.shape[0] == 5: # (case 2)
            self.nominal_values /= 100.0  # unit change [%] -> [-]
            self.nominal_values[3] *= 100.0  # undo unit change for BMP

        # Computation of remaining nominal values (case 1)
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

    def to_latex(self):
        xi = deepcopy(self.xi)
        uncertain_xis = self.get_uncertain_xi_ch_pr_li()
        xi[5] = uncertain_xis[0].nominal_value
        xi[7] = uncertain_xis[1].nominal_value
        xi[8] = uncertain_xis[2].nominal_value

        xi_names = [r"$\xi_{" + str(i + 1) + r"}$" for i in range(len(self.xi_names))]

        df_dict = {r"$\xi$": xi_names, r"$[g/l]$": xi}
        df = pd.DataFrame(df_dict)
        df = df.loc[(df != 0).all(axis=1)]

        df = df.round(decimals=3)
        df[r"$[g/l]$"] = df[r"$[g/l]$"].astype(str)
        df.index.name = None
        df.T.to_latex(
            buf=f"./tables/{self.name}.tex",
            caption=f"{self.name.lower().replace('_', ' ')}" + r" $\xi$ values",
            label=f"tab:{self.name}",
            header=False,
            index=True,
            position="H",
        )


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
            5.52,
            10.04,
            7.4,
            14.51,
            1.94,
        ]
    ),
    xi=[
        10.32,
        0,
        0,
        0.764,
        662.7136196,
        239.7542835,
        0,
        26.3338422,
        7.991710817,
        0.306324523,
        0.016122343,
        14.83951389,
        - 0.028103425,
        1.019751485,
        0,
        6.50282E-06,
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
            372.0,
            31.7409815648268,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.52,
            10.04,
            7.4,
            14.51,
            1.94,
        ]
    ),
    xi=[
        10.44,
        0,
        0,
        1.574,
        682.5901844,
        199.9129515,
        0,
        42.28337446,
        7.633251677,
        0.26796237,
        0.014103283,
        35.3441635,
        - 0.001592555,
        5.460218849,
        0,
        0.000133961,
        0,
        0,

    ],
    state="solid",
    cost=35.0,
)


CATTLE_MANURE = Substrate(
    name="CATTLE_MANURE",
    nominal_values=np.array(
        [
            246.0,
            8.08403107074424,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.52,
            10.04,
            7.4,
            14.51,
            1.94,
        ]
    ),
    xi=[
        4.534935181,
        0,
        0,
        1.303128034,
        919.1596893,
        20.81796874,
        0,
        13.31346487,
        2.006047557,
        0.058613232,
        0.003084907,
        19.14217214,
        0.023350928,
        4.534182693,
        0,
        0.415286236,
        0,
        0,

    ],
    state="liquid",
    cost=20.0,
)

CATTLE_MANURE_VERY_UNCERTAIN = deepcopy(CATTLE_MANURE)
CATTLE_MANURE_VERY_UNCERTAIN.name = "CATTLE_MANURE_VERY_UNCERTAIN"
CATTLE_MANURE_VERY_UNCERTAIN.std_devs *= 2.5

SUGAR_BEET_SILAGE = Substrate(
    name="SUGAR_BEET_SILAGE",
    nominal_values=np.array(
        [
            389.0,
            39.2745959278011,
        ]
    ),
    variation_coefficients=np.array(
        [
            5.52,
            10.04,
            7.4,
            14.51,
            1.94,
        ]
    ),
    xi=[
        8.17,
        0,
        0,
        0.07,
        607.2540407,
        327.3251763,
        0,
        9.572991838,
        0.607880085,
        0.346182296,
        0.018220121,
        28.34354217,
        0.01227284,
        0.990980133,
        0,
        7.50076E-07,
        0,
        0,

    ],
    state="solid",
    cost=50.0,
)


def multiple_substrates_to_latex(subs: list[Substrate]):

    column_names = [r"$\xi_{" + f"{i}" + r"}$" for i in range(1, 19)]

    df = pd.DataFrame(columns=column_names)

    for sub in subs:
        xi = deepcopy(sub.xi)
        uncertain_xis = sub.get_uncertain_xi_ch_pr_li()
        xi[5] = uncertain_xis[0].nominal_value
        xi[7] = uncertain_xis[1].nominal_value
        xi[8] = uncertain_xis[2].nominal_value

        xi_names = deepcopy(sub.xi_names)

        xi_names = [r"$\xi_{" + xi_name[1:] + r"}$" for xi_name in xi_names]

        df_dict = {}
        for col_name, xi_ in zip(column_names, xi):
            df_dict[col_name] = xi_

        sub_name = sub.name.lower()
        if sub_name == "cattle_manure":
            sub_name = "manure"
        elif sub_name == "cattle_manure_very_uncertain":
            sub_name = "dictated_cattle_manure"
        sub_name = sub_name

        sub_name_list = sub_name.split("_")

        sub_name = (
            r"\begin{tabular}[c]{@{}l@{}}"
            + r"\\".join(sub_name_list)
            + r"\end{tabular}"
        )

        # sub_name = r"\begin{tabular}[c]{@{}l@{}}corresponding\\ state\end{tabular}"

        df_temp = pd.DataFrame(df_dict, index=[sub_name])
        df = pd.concat([df, df_temp])

    df = df.loc[:, (df != 0).all(axis=0)]
    df = df.round(decimals=3)

    df_dict = {}
    for col_name, xi_name in zip(column_names, subs[0].xi_names):
        xi_name = "$" + xi_name[0] + "_\mathrm{" + xi_name[1:] + "}$"
        df_dict[col_name] = xi_name
    df_temp = pd.DataFrame(
        df_dict,
        index=[r"\begin{tabular}[c]{@{}l@{}}corresponding\\ state\end{tabular}"],
    )

    df = pd.concat([df_temp, df])

    df = df.dropna(axis=1)

    df = df.T
    df = df.astype(str)
    df.index.name = None
    df.to_latex(
        buf=f"./tables/xi_values.tex",
        caption=r"Values of inlet concentrations '$\xi$' [\SI{}{\gram \per \liter}]",
        label=f"tab:xi_values",
        header=True,
        index=True,
        position="H",
    )


if __name__ == "__main__":
    multiple_substrates_to_latex(
        [
            CORN_SILAGE,
            GRASS_SILAGE,
            SUGAR_BEET_SILAGE,
            CATTLE_MANURE,
            CATTLE_MANURE_VERY_UNCERTAIN,
        ]
    )
    # CATTLE_MANURE.to_latex()
    # CORN_SILAGE.to_latex()
    # GRASS_SILAGE.to_latex()
    # SUGAR_BEET_SILAGE.to_latex()
    # CATTLE_MANURE_VERY_UNCERTAIN.to_latex()
