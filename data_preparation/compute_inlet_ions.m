%% Version: Matlab R2022b, Update 6
% Author: Simon Hellmann
% first created: 2024-04-10

function [S_ac_minus,S_nh3,S_ion] = compute_inlet_ions(pH,S_ac,S_IN)

% Info: compute inlet concentrations of ionic species from pH assuming no 
% dissolved co2 in the substrate

% Outputs: 
% S_ac      - mass concentration of total acetic acid
% S_ac_minus- mass concentration (MC) of dissociated acetic acid 
% S_nh3     - MC of ammonia
% S_ion     - molar concentration of left-over cations to close ion balance

% Inputs: 
% pH        - pH value of substrate
% S_ac      - GC measurement of acetic acid [g/L]
% S_IN      - total inorganic nitrogen concentration [g/L_FM]

% get dissociation and constants and ion product of water: 
K_a_IN = 1.3490E-9; % from Sörens Diss, Tab. B.2 (values therein are correct, K_a_IN and K_a_co2 are false in his github!)
K_a_ac = 1.7378E-5;
K_w = 2.07877E-14; 

% compute dissociated versions of ions (assume dissociation equilibrium):
SHPlus = 10^(-pH);  
S_ac_minus =K_a_ac*S_ac/(K_a_ac + SHPlus); 
S_nh3 =     K_a_IN*S_IN/(K_a_IN + SHPlus); 

% compute remaining ion concentration: 
Phi = K_w/SHPlus - SHPlus; % ion balance (can also be negative)
S_nh4 = S_IN - S_nh3;       % free ammonium nitrogen
S_ion = Phi - S_nh4/17 + S_ac_minus/60; 

end