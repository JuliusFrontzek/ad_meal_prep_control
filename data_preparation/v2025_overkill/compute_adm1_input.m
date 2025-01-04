%% Version: Matlab R2024a, Update 3
% Author: Simon Hellmann
% first created: 2025-01-03

function x_inlet_default = compute_adm1_input(substrate_data_means, ...
                            substrate_name,model_structure, ...
                            flag_CH_method)

% Info: computes default input vector of the ADM1 models -R3 & -R4 based on
%       mean substrate characterizations at DBFZ lab

% Outputs: 
% - x_inlet_default:    vector of default inlet concentrations (size depends
%                       on model structure)

% Inputs: 
% - substrate_data_means: table of mean DBFZ lab measurements of specific substrate
% - substrate_name:     name of substrate 
% - model_structure:    'R3', 'R3-Core', 'R3-frac', 'R4'
% - flag_CH_method:     1: compute carbs via BMP. 2: via Lignin

%% BMP values acc. to KTBL: 
literature_table = readtable('../../data/data_in/literature_values_agr_substrates.xlsx');
row_idx_sub = strcmp(literature_table.substrate, substrate_name); 
BMP = literature_table{row_idx_sub,"bmp"}; 

%% all inlet values assumed as zero
% assume no gas release in substrate:
S_ch4_gas = 0; 
S_co2_gas = 0; 
% assume that no gasses dissolved in liquid phase: 
S_ch4 = 0; 
S_IC = 0; 
S_hco3 = 0;                 % because S_IC is also zero
% assume no microbial biomass in substrate: 
X_bac = 0; 
X_ac = 0; 

%% easily computed inlet concentrations
rho_FM = 1;  % assumed mass density of free matter [kg_FM/L_FM]
TS = substrate_data_means.TS/100;       % unit change [%] -> [-]

NH4_N = substrate_data_means.NH4_N;     % ammonium nitrogen [g/L]
S_IN = NH4_N*rho_FM;            % unit change [g/kg FM] -> [g/L FM]

XA = substrate_data_means.XA;   % ash
X_ash = XA*TS*rho_FM;           % unit change [g/kg_TS] -> [g/L_FM]
S_h2o = rho_FM*(1-TS)*1000;     % water; unit change [%] -> [g/L_FM]

%% macro nutrients: 
XP = substrate_data_means.XP;   % raw protein [g/kg_TS]
XL = substrate_data_means.XL;   % raw fats [g/kg_TS]
XC = 1000 - XA - XP - XL;       % raw carbohydrates (Diss. Weinrich, Tab.3.4) [g/kg_TS]
X_pr = XP*TS*rho_FM;            % unit change [g/kg_TS] -> [g/L_FM]
X_li = XL*TS*rho_FM;            % unit change [g/kg_TS] -> [g/L_FM]

% compute carbs' fermentability via BMP (Method 1): 
BMP_stoich = 420; % max. BMP for agricultural substrates [L_N/kg_FoTS]
FQ_ges = BMP/BMP_stoich;    % total fermentability
FQ_ch = XC^(-1)*(FQ_ges*(1000 - XA) - XP - XL); % fermentability of carbs

% compute carbs' fermentability via lignin (Method 2, optional): 
if flag_CH_method == 2 && ~isnan(substrate_data_means.XLig) % only pursue if mean lignin value exist
    XLig = substrate_data_means.XLig;    % lignin content [g/kg_TS]
    FQ_ch = (XC - XLig)/XC; 
end

X_ch = FQ_ch*XC*TS*rho_FM; % unit change [g/kg_TS] -> [g/L_FM]

% % compute microbial biomasses (Félix Delory estimates it like this)
% X_mic_biomass = 0.001*(XC + XP + XL);   % total biomass is 0,1% fraction of macro nutrients
% X_bac = 0.95*X_mic_biomass*TS*rho_FM;   % assign 95% of biomass to X_bac...
% X_ac =  0.05*X_mic_biomass*TS*rho_FM;   % ...and 5% to X_ac

%% VFAs

% preprocess GC measurements (summarize i/n-BU & i/n-VA, dilution, unit change):
dilution = substrate_data_means.dilution; 
ac_eff = substrate_data_means.AC * dilution; 
pro_eff = substrate_data_means.PRO * dilution; 
i_bu_eff = substrate_data_means.i_BU * dilution; 
n_bu_eff = substrate_data_means.n_BU * dilution; 
bu_eff = i_bu_eff + n_bu_eff; % total butyric acid
i_va_eff = substrate_data_means.i_VA * dilution; 
n_va_eff = substrate_data_means.n_VA * dilution; 
va_eff = i_va_eff + n_va_eff; % total valeric acid
% apply unit change [mg/L] -> [g/L] to acids: 
S_ac = ac_eff/1000;    % acetic acid
S_pro = pro_eff/1000;  % proprionic acid
S_bu = bu_eff/1000;    % butyric acid
S_va = va_eff/1000;    % valeric acid

% aggregate VFAs to equivalent AC contrentration based on molar masses: 
% molar masses [kg/kmol]: 
M_ac = 60.05; 
M_pro = 74.08;
M_bu = 88.1;
M_va = 102.13; 

S_ac_eq = S_ac + S_pro/M_pro*M_ac + S_bu/M_bu*M_ac + S_va/M_va*M_ac;

%% ions based on substrate's pH:
pH = substrate_data_means.pH; 
if isnan(pH)
    pH = literature_table{row_idx_sub,"ph"}; 
end
[S_ac_minus,S_nh3,S_ion] = compute_inlet_ions(pH,S_ac_eq,S_IN); 

%% build vector of inlet concentrations 
% first based on ADM1-R3
x_inlet_default_r3 = [S_ac_eq,S_ch4,S_IC,S_IN,S_h2o,X_ch,X_pr,X_li,X_bac,...
            X_ac,X_ash,S_ion,S_ac_minus,S_hco3,S_nh3,S_ch4_gas,S_co2_gas]';

% correct for model structure: 
x_inlet_default = correct_adm1_input_for_model_structure(x_inlet_default_r3,model_structure); 

end
