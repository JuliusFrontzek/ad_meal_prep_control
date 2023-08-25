%% Version
% (R2022b) Update 5
% Erstelldatum: 18.08.2023
% Autor: Simon Hellmann

% compute inlet concentrations of ions from pH assuming no dissolved co2 in
% substrate

function [S_ac_minus, S_nh3, S_ion] = computeInletConcentrationsIons(pH,Ac,NH4N)
% outputs: 
% S_ac_minus- mass concentration (MC) of dissociated acetic acid 
% S_nh3     - MC of ammonia
% S_ion     - molar concentration of left-over ions to close ion balance

% inputs: 
% pH        - pH value of substrate (assumed the same as in reactor)
% Ac        - GC measurement of acetic acid [mg/L]  
% NH4N      - free ammonium mass fraction [g/kg_FM]

% define constants: 
K_a_IN = 1.3809E-9; % dissociation constants
K_a_ac = 1.7378E-5;
K_w = 2.07877105595436e-14; 
rhoFM = 1;    % assumed mass density of fresh matter [kg/L]

% transform lab measurements into concentrations: 
Sac = Ac/1000;      % unit change [mg/L] -> [g/L_FM]
SIN = NH4N*rhoFM;   % unit change [g/kgFM] -> [g/L_FM]

% compute dissociated versions of ions:
SHPlus = 10^(-pH);  
S_ac_minus = K_a_ac*Sac/(K_a_ac + SHPlus); 
S_nh3 = K_a_IN*SIN/(K_a_IN + SHPlus); 

% compute remaining ions (these can also be negative!): 
Phi = K_w/SHPlus - SHPlus; % ion balance
S_ion = Phi - NH4N/17 + S_ac_minus/60; 

end