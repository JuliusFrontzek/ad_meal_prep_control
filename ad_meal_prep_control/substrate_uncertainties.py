from uncertainties import ufloat
import numpy as np
import pandas as pd


def xc_(xa, xp, xl):
    return 1000.0 - xa - xp - xl


def delta_xc_(delta_xa, delta_xp, delta_xl):
    return np.sqrt(delta_xa**2 + delta_xp**2 + delta_xl**2)


def fq_ch_(xc, fq_ges, xa, xp, xl):
    return xc ** (-1) * (fq_ges * (1000.0 - xa) - xp - xl)


def delta_fq_ch_(
    xc, fq_ges, xa, xp, xl, delta_xc, delta_fq_ges, delta_xa, delta_xp, delta_xl
):
    return np.sqrt(
        (-(xc ** (-2.0)) * (fq_ges * (1000 - xa) - xp - xl) * delta_xc) ** 2
        + (xc ** (-1.0) * (1000 - xa) * delta_fq_ges) ** 2
        + (xc ** (-1.0) * fq_ges * delta_xa) ** 2
        + (xc ** (-1.0) * delta_xp) ** 2
        + (xc ** (-1.0) * delta_xl) ** 2
    )


def x_ch_in_(fq_ch, xc, ts, rho_fm):
    return fq_ch * xc * ts * rho_fm


def x_pr_in_(fq_pr, xp, ts, rho_fm):
    return fq_pr * xp * ts * rho_fm


def x_li_in_(fq_li, xl, ts, rho_fm):
    return fq_li * xl * ts * rho_fm


def f_q_qes_(bmp):
    """
    Bestimmt Gesamtfermentierbarkeit.

    Nur für landwirtschaftliche Substrate gültig.
    """
    BMP_STOIC = 420  # 420l_CH4/kg_FoTS nach Weißbach für landwirtschaftliche Substrate
    return bmp / BMP_STOIC


if __name__ == "__main__":
    df = pd.read_excel("the_document.ods", engine="odf")
    xa = 50.0
    xp = 60.0
    xl = 100.0
    xc = 150.0
    fq_ges = 200.0
    ts = 50.0

    delta_xa = 5.0
    delta_xp = 200.0
    delta_xl = 7.0
    delta_xc = 8.0
    delta_fq_ges = 3.0
    delta_ts = 5.0

    xa_ufloat = ufloat(xa, delta_xa)
    xp_ufloat = ufloat(xp, delta_xp)
    xl_ufloat = ufloat(xl, delta_xl)
    xc_ufloat = ufloat(xc, delta_xc)
    ts_ufloat = ufloat(ts, delta_ts)
    fq_ges_ufloat = ufloat(fq_ges, delta_fq_ges)

    print(xc_(xa_ufloat, xp_ufloat, xl_ufloat).std_dev)
    print(delta_xc_(delta_xa, delta_xp, delta_xl))

    print(fq_ch_(xc_ufloat, fq_ges_ufloat, xa_ufloat, xp_ufloat, xl_ufloat).std_dev)
    print(
        delta_fq_ch_(
            xc, fq_ges, xa, xp, xl, delta_xc, delta_fq_ges, delta_xa, delta_xp, delta_xl
        )
    )

    print(x_ch_in_(fq_ges_ufloat, xc_ufloat, ts_ufloat, 1.0))
