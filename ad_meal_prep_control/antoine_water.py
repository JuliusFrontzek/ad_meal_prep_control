# Antoine equation; parameters applicable for water in between 1°C and 100°C
# Source: https://en.wikipedia.org/wiki/Antoine_equation
A = 8.07131
B = 1730.63
C = 233.426


def vapour_pressure_h2o(t: float):
    """
    Params:
        t: Temperature in K
    Returns:
        p: The vapour pressure of water in bar.
    """
    return 10 ** (A - B / (C + (t - 273.15))) * 133.321992076 * 1.0e-5
