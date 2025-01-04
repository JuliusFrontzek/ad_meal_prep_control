%% Version
% (R2022b) Update 5
% Erstelldatum: 18.08.2023
% last modified: 12.12.2023
% Autor: Simon Hellmann

% compute inlet concentrations of ions from pH assuming no dissolved co2 in
% substrate

function [S_ac_minus,S_nh3,S_ion] = computeInletIons(pH,S_ac,S_IN)
% outputs: 
% S_ac      - mass concentration of total acetic acid
% S_ac_minus- mass concentration (MC) of dissociated acetic acid 
% S_nh3     - MC of ammonia
% S_ion     - molar concentration of left-over cations to close ion balance

% inputs: 
% pH        - pH value of substrate (assumed the same as in reactor)
% S_ac      - GC measurement of acetic acid [g/L]
% S_IN      - total inorganic nitrogen concentration [g/L_FM]

% dissociation constants: 
K_a_IN = 1.3490E-9; % from Sörens Diss, Tab. B.2 (values therein are correct, K_a_IN and K_a_co2 are false in his github!)
K_a_ac = 1.7378E-5;
K_w = 2.07877e-14; 

% compute dissociated versions of ions (assume dissociation equilibrium):
SHPlus = 10^(-pH);  
S_ac_minus = K_a_ac*S_ac/(K_a_ac + SHPlus); 
S_nh3 = K_a_IN*S_IN/(K_a_IN + SHPlus); 

% compute remaining ions (these can also be negative!): 
Phi = K_w/SHPlus - SHPlus; % ion balance

% compute remaining ion concentration: 
S_nh4 = S_IN - S_nh3;       % free ammonium nitrogen
S_ion = Phi - S_nh4/17 + S_ac_minus/60; 

end