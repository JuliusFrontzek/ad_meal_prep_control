# kinetic constants [1/d] (Tab. B.7 in Sörens Diss):
kchF = 0.25  # selbst gewählt
kchS = kchF / 10.0
kpr = 0.2
kli = 0.1
kdec = 0.02
fracChFast = (
    0.5  # fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)
)
mu_m_ac = 0.4  # 1/d
K_S_ac = 0.14  # kg/m^3
K_I_nh3 = 0.0306  # kg/m^3

# Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011
Kco2 = 0.025

# miscellaneous parameters (Weinrich, 2017, Tab. B.7):
R = 0.08315  # id. gas constant [bar l/mol/K]
ph2o = 0.0657  # partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0133  # atmospheric pressure [bar]
kla = 200  # mass transfer coefficient [1/d]
kp = 5e4  # friction parameter [l/bar/d]
T = 311  # operating temperature [K]
Vl = 163.0  # liquid volume, aus Sörens GitHub [m^3]
Vg = Vl / 10.0  # 26000.0  # gas volume, aus Sörens GitHub [m^3]
rho = 1000.0  # mass density of digestate [kg/m^3]
Mch4 = 16  # molar mass CH4 [kg/kmol]
Mco2 = 44  # molar mass CO2 [kg/kmol]

pHULac = 7.0
pHLLac = 6.0
nac = 3.0 / (pHULac - pHLLac)
KW = 2.078771055954360 * 10 ** (-14.0)  # mol/l
KS_IN = 0.0017  # kg/m^3
k_AB_ac = 1.0 * 10**10  # l/mol/d
k_AB_co2 = 1.0 * 10**10  # l/mol/d
k_AB_IN = 1.0 * 10**10  # l/mol/d
K_a_ac = 1.737800828749374 * 10 ** (-5.0)  # mol/l
K_a_co2 = 4.937073397534361 * 10 ** (-7.0)  # mol/l
K_a_IN = 1.110286652708067 * 10 ** (-9.0)  # mol/l

K_a_IN = 1.3809e-9
K_a_ac = 1.7378e-5
K_a_co2 = 5.0981e-7

# Parameters introduced by Julius
p_gas_storage = p0 * 1.001  # slightly larger than atmospheric pressure [bar]
p_norm = p0  # [bar]
T_norm = 273.15  # [K]
T_gas_storage = 323.15  # [K]
P_el_chp = 124.0  # kW
V_GAS_STORAGE_MAX = P_el_chp * 4.3  # m^3
