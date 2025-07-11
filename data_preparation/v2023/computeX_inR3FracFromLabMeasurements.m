%% Version
% (R2022b) Update 6
% Erstelldatum: 24.08.2023
% last modified: 12.12.2023
% Autor: Simon Hellmann

function [xIn,TS] = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich, ...
    substrate_lit_values)

% uses lab data of one sample and transforms it into inlet concentrations
% of ADM1-R3-frac model and TS

% xIn       - vector of inlet concentrations for ADM1-R3-frac 
% TS        - total solids content of sample [-]
% lab data  - results from lab analysis of one substrate sample
% BMP       - BMP of specific substrate acc. to KTBL
% BMP_stoich- max. BMP for agricultural substrates [L_N/kg_FoTS]
% substrate_lit_values - contains literature values for specific substrates
    
% compute inlet concentrations of all states except ions 
rhoFM = 1;  % assumed mass density of free matter [kg_FM/l_FM]
TS = labData(1)/100;    % unit change [%] -> [-]
   
% assume no gas release in substrate:
S_ch4_gas = 0; 
S_co2_gas = 0; 
% assume that no gasses dissolved in liquid phase: 
S_ch4 = 0; 
S_IC = 0; 
S_hco3 = 0;             % since S_IC is also zero

%% process values that are non-zero:
XA = labData(4); 
X_ash = XA*TS*rhoFM;    % unit change [g/kg_TS] -> [g/l_FM]

water = labData(2); 
S_h2o = water/100*1000*rhoFM;  % unit change [%FM] -> [g/l_FM]

% macro nutrients: 
XP = labData(5); 
X_pr = XP*TS*rhoFM;     % unit change [g/kg_TS] -> [g/l_FM]
XL = labData(6); 
X_li = XL*TS*rhoFM;     % unit change [g/kg_TS] -> [g/l_FM]
X_ch_slow = 0;          % all inflowing carbohydrates assigned to fast fraction
% compute fast carbohydrates in via BMP
FQ_ges = BMP/BMP_stoich; 
XC = labData(7); 
FQ_ch = XC^(-1)*(FQ_ges*(1000 - XA) - XP - XL); 
X_ch_fast = FQ_ch*XC*TS*rhoFM; 

% compute bacteria (estimated acc. to Félix)
X_biomass = 0.001*(XC + XP + XL);   % total biomass is 0,1% fraction of macro nutrients
X_bac = 0.95*X_biomass*TS*rhoFM;    % assign 95% of biomass to X_bac...
X_ac = 0.05*X_biomass*TS*rhoFM;     % ...and 5% to X_ac

% lab data measurements:
NH4N = labData(3);          % ammonium nitrogen [g/L]
S_ac = labData(10)/1000;    % acetic acid [g/L]
pH = 7.4;   % default inlet pH about netral

% hier Sonderwerte abrufen, sofern verfügbar:
if ~isempty(substrate_lit_values)
    pH = substrate_lit_values.pH; 
    if isfield(substrate_lit_values,'NH4N')
        NH4N = substrate_lit_values.NH4N; 
    end
    if isfield(substrate_lit_values,'ACeq')
        S_ac = substrate_lit_values.ACeq;
    end
end

S_IN = NH4N*rhoFM;      % unit change [g/kg FM] -> [g/L FM]

% compute inlet concentrations of ions from pre-defined pH:
[S_ac_minus,S_nh3,S_ion] = computeInletIons(pH,S_ac,S_IN); 

% summarize all in one vector of inlet concentrations: 
xIn = [S_ac, S_ch4, S_IC, S_IN, S_h2o, X_ch_fast, X_ch_slow, X_pr, X_li,...
       X_bac, X_ac, X_ash, S_ion, S_ac_minus, S_hco3, S_nh3, S_ch4_gas, S_co2_gas]';

end
