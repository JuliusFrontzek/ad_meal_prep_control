# kinetic constants [1/d] (Tab. B.7 in Sörens Diss):
kchS = 0.25
kchF = 10 * kchS  # selbst gewählt
kpr = 0.2
kli = 0.1
kdec = 0.02
fracChFast = 0.25  # fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)
mu_m_ac = 0.4  # 1/d
K_S_ac = 0.14  # g/l
K_I_nh3 = 0.0306  # g/l

# Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011
Kco2 = 0.025

# miscellaneous parameters (Weinrich, 2017, Tab. B.7):
R = 0.08315  # id. gas constant [bar l/mol/K]
ph2o = 0.0657  # partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0133  # atmospheric pressure [bar]
kla = 200  # mass transfer coefficient [1/d]
kp = 5e4  # friction parameter [l/bar/d]
T = 273.15  # operating temperature [K]
Vl = 100  # liquid volume, aus Sörens GitHub [l]
Vg = 10  # gas volume, aus Sörens GitHub [l]
rho = 1  # mass density of digestate [kg/l]
Mch4 = 16  # molar mass CH4 [kg/kmol]
Mco2 = 44  # molar mass CO2 [kg/kmol]

pHULac = 7.0
pHLLac = 6.0
nac = 3.0 / (pHULac - pHLLac)
KW = 1.0 * 10 ** (-13.7)  # mol/l
KS_IN = 0.0017  # g/l
k_AB_ac = 1.0 * 10**10  # l/mol/d
k_AB_co2 = 1.0 * 10**10  # l/mol/d
k_AB_IN = 1.0 * 10**10  # l/mol/d
K_a_ac = 1.0 * 10 ** (-4.76)  # mol/l
K_a_co2 = 1.0 * 10 ** (-6.29)  # mol/l
K_a_IN = 1.0 * 10 ** (-8.87)  # mol/l
