%% Version
% (R2022b) Update 5
% Erstelldatum: 24.08.2023
% Autor: Simon Hellmann

function xIn = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich)

% used lab data of one sample and transforms it into inlet concentrations
% of ADM1-R3-frac model

% xIn       - vector of inlet concentrations for ADM1-R3-frac 
% lab data  - results from lab analysis of one substrate sample
% BMP       - BMP of specific substrate acc. to KTBL
% BMP_stoich- max. BMP for agricultural substrates [L_N/kg_FoTS]
    
% compute inlet concentrations of all states except ions 
rhoFM = 1;  % assumed mass density of free matter [kg_FM/l_FM]
TS = labData(1)/100;        % unit change [%] -> [-]
   
% assume no gas release in substrate:
S_ch4_gas = 0; 
S_co2_gas = 0; 
% assume that no gasses dissolve in liquid phase: 
S_ch4 = 0; 
S_IC = 0; 
S_hco3 = 0;             % since S_IC is also zero

NH4N = labData(3); 
S_IN = NH4N*rhoFM;      % unit change [g/kg_FM] -> [g/l_FM]

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

% compute bacteria: 
X_bac = 0.95*(XC + XP + XL)*TS*rhoFM;   % estimated acc. to Félix
X_ac = 0.05*(XC + XP + XL)*TS*rhoFM;    % estimated acc. to Félix

% acetic acid:
Ac = labData(8); 
S_ac = Ac/1000;         % unit change [mg/l] -> [g/l]

% compute inlet concentrations of ions from pre-defined pH:
pH = 7.5;   % regular operation condition
[S_ac_minus, S_nh3, S_ion] = computeInletConcentrationsIons(pH,Ac,NH4N); 

% summarize all in one vector of inlet concentrations: 
xIn = [S_ac, S_ch4, S_IC, S_IN, S_h2o, X_ch_fast, X_ch_slow, X_pr, X_li,...
       X_bac, X_ac, X_ash, S_ion, S_ac_minus, S_hco3, S_nh3, S_ch4_gas, S_co2_gas]';

end

