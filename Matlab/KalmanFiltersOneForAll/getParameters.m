%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  params = getParameters(flagModel, flagFrac, system, parameters)
% returns model parameters summarized in struct params
% params - struct containing model parameters as arrays aNum, cNum and thNum
% flagModel -   3: ADM1-R3; 4: ADM1-R4
% flagFrac -    0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% system -      table containing system parameters V_liq, V_gas and p_atm
% parameters    struct containing model parameters for given model

% time-invariant parameters
s = system.Variables;   % s = systemParameters
% extract system parameter values: 
V_liq = s(1);  % gas volume
V_gas = s(2);  % liquid volume
% p_atm = s(3);  % atmospheric pressure [bar] XY: Achtung: Typo im GitHub
p_atm = 1.0133; 

switch flagModel
    case 3
        % renaming for easier understanding: 
%         systemInput = input.ADM1_R3.Variables;
        parameters = parameters.ADM1_R3.Variables; % modelParameters
        
        % extract model parameter values:  
        K_H_ch4 =  parameters(1);  % Henry parameter ch4 [mol/l/bar]
        K_H_co2 =  parameters(2);  % Henry parameter co2 [mol/l/bar]
        K_S_IN =  parameters(3);  % half-saturation constant nitrogen limitation [g/l]
        K_I_nh3 =  parameters(4); % ammonia inbibition constant [g/l]
        % KaIN =  parameters(5);  % dissociation constant ammonia [mol/l]
        % Kaac =  parameters(6);  % dissociation constant acetic acid [mol/l]
        % Kaco2 =  parameters(7); % dissociation constant carbonate [mol/l]
        % corrected values acc. to Simons computation (Sören made slight mistake):
        K_a_IN = 1.3809E-9; 
        K_a_ac = 1.7378E-5; 
        K_a_co2 = 5.0981E-7;
        K_S_ac =  parameters(8);  % half-saturation constant acetoclastic methanogenesis [g/l]
        K_w =  parameters(9);    % ion product of water [mol/l]
        R =  parameters(10);    % universal gas constant [bar l/mol/K]
        T =  parameters(11);    % operating temperature [K]
        k_AB_IN =  parameters(12);    % kin. dissociation constant ammonia [l/mol/d] 
        k_AB_ac =  parameters(13);    % kin. dissociation constant acetic acid [l/mol/d]
        k_AB_co2 =  parameters(14);   % kin. dissociation constant carbonate [l/mol/d]
        k_La =  parameters(15);  % mass transfer coefficient [1/d]
        k_ch =  parameters(16);  % hydrolysis constant carbohydrates [1/d]
        k_dec =  parameters(17); % hydrolysis constant biomass decay [1/d]
        k_li =  parameters(18);  % hydrolysis constant lipids [1/d]
        k_m_ac =  parameters(19);  % max. growth rate \mu_m [1/d]
        k_p =  parameters(20);   % friction parameter [l/bar/d]
        k_pr =  parameters(21);  % hydrolysis constant proteins [1/d]
        pK_l_ac =  parameters(22); % lower pH boundary (formerly pHLL)  
        pK_u_ac =  parameters(23); % upper pH boundary (formerly pHUL)
        p_h2o =  parameters(24); % partial pressure of water in gas phase (saturated) [bar]
        
        rho = 1000;        % mass density of digestate [g/l]
        Mch4 = 16;      % molar mass CH4 [kg/kmol]
        Mco2 = 44;      % molar mass CO2 [kg/kmol]
        n_ac = 3/(pK_u_ac - pK_l_ac); 
        
        % order model parameters in the rights structures (prepare simulation)
        c1 = 1/V_liq; 
        c2 = n_ac; 
        c3 = 10^(-(3/2)*(pK_u_ac + pK_l_ac)/(pK_u_ac - pK_l_ac)); 
        c4 = 4*K_w; 
        c5 = k_La; 
        c6 = k_La*K_H_ch4*R*T; 
        c7 = k_La*K_H_co2*R*T; 
        c8 = K_S_IN; 
        c9 = k_AB_ac;
        c10 = k_AB_co2; 
        c11 = k_AB_IN; 
        c12 = k_La*V_liq/V_gas; 
        c13 = k_p/p_atm*(R*T/Mch4)^2;
        c14 = 2*k_p/p_atm*(R*T)^2/Mch4/Mco2;
        c15 = k_p/p_atm*(R*T/Mco2)^2;
        c16 = k_p/p_atm*R*T/Mch4*(2*p_h2o - p_atm); 
        c17 = k_p/p_atm*R*T/Mco2*(2*p_h2o - p_atm); 
        c18 = k_p/p_atm*(p_h2o - p_atm)*p_h2o; 
        c19 = R*T/Mch4;
        c20 = R*T/Mco2;
        c21 = rho; 
        c22 = -k_p/V_gas/p_atm*(R*T/Mch4)^2;
        c23 = -2*k_p/V_gas/p_atm*(R*T)^2/Mch4/Mco2;
        c24 = -k_p/V_gas/p_atm*(R*T/Mco2)^2;
        c25 = -k_p/V_gas/p_atm*(R*T/Mch4)*(2*p_h2o - p_atm);
        c26 = -k_p/V_gas/p_atm*(R*T/Mco2)*(2*p_h2o - p_atm);
        c27 = -k_La*V_liq/V_gas*K_H_ch4*R*T - k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;
        c28 = -k_La*V_liq/V_gas*K_H_co2*R*T - k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;
        c29 = k_AB_ac*K_a_ac; 
        c30 = k_AB_co2*K_a_co2;
        c31 = k_AB_IN*K_a_IN; 
        c32 = V_liq/V_gas;
        c33 = -k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;

        % combine all in column vector: 
        cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,...
                c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33]';
    
    case 4
        % code-Abschnitt für R4 und R4-frac noch überprüfen!
        parameters = parameters.ADM1_R4.Variables; % modelParameters
        
        k_ch = parameters(6); 
        k_pr = parameters(10); 
        k_li = parameters(8);
        k_dec = parameters(7);  
        if flagFrac == 1
            k_ch_F = k_ch;         % war vorher der Wert für kch
            k_ch_S = 1E-1*k_ch_F;   % selbst gewählt
            fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)
        end

        % Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
        K_H_ch4 = parameters(1);      
        K_H_co2 = parameters(2); 
        
        % miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
        R = parameters(3);    % id. gas constant [bar l/mol/K]
        p_h2o = parameters(11);  % 0.0657 partial pressure of water in gas phase (saturated) [bar]
%         p_h2o = 0;              % normal conditions [bar]
        % check if pH2O is zero!
        k_La = parameters(5);       % mass transfer coefficient [1/d]
        k_p = parameters(9);       % friction parameter [l/bar/d]
        T = parameters(4);     % operating temperature [K]
        % XY: check if normal conditions!
        rho = 1;        % mass density of digestate [kg/l]
        Mch4 = 16;      % molar mass CH4 [kg/kmol]
        Mco2 = 44;      % molar mass CO2 [kg/kmol]
        
        %% order model parameters in the rights structures (prepare simulation)
        c1 = 1/V_liq; 
        c2 = k_La; 
        c3 = k_La*K_H_ch4*R*T; 
        c4 = k_La*K_H_co2*R*T; 
        c5 = k_La*V_liq/V_gas; 
        c6 = k_p/p_atm*(R*T/Mch4)^2;
        c7 = 2*k_p/p_atm*(R*T)^2/Mch4/Mco2;
        c8 = k_p/p_atm*(R*T/Mco2)^2;
        c9 = k_p/p_atm*R*T/Mch4*(2*p_h2o - p_atm); 
        c10 = k_p/p_atm*R*T/Mco2*(2*p_h2o - p_atm); 
        c11 = k_p/p_atm*(p_h2o - p_atm)*p_h2o; 
        c12 = R*T/Mch4;
        c13 = R*T/Mco2;
        c14 = rho; 
        c15 = -k_p/p_atm/V_gas*(R*T/Mch4)^2;
        c16 = -2*k_p/p_atm/V_gas*(R*T)^2/Mch4/Mco2;
        c17 = -k_p/p_atm/V_gas*(R*T/Mco2)^2;
        c18 = -k_p/p_atm/V_gas*(R*T/Mch4)*(2*p_h2o - p_atm);
        c19 = -k_p/p_atm/V_gas*(R*T/Mco2)*(2*p_h2o - p_atm);
        c20 = -k_La*V_liq/V_gas*K_H_ch4*R*T - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
        c21 = -k_La*V_liq/V_gas*K_H_co2*R*T - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
        c22 = V_liq/V_gas;
        c23 = -k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
        % combine all in column vector: 
        cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23]';

end

% time-variant parameters
switch flagModel
    case 3
        if flagFrac == 0
            thNum = [k_ch, 0, k_pr, k_li, k_dec, k_m_ac, K_S_ac, K_I_nh3, 0]';
        else 
            % XY: Fall R3-frac nachtragen
        end
    
    case 4
        if flagFrac == 0
            thNum = [k_ch,  0,      k_pr, k_li, k_dec, 0]'; 
        else
            thNum = [k_ch_F,k_ch_S, k_pr, k_li, k_dec, fracChFast]'; 
        end

end

% petersen matrix acc. to Weinrich 2021
switch flagModel
    case 3
        if flagFrac == 0
            % Note:all aij with |aij ~= 1|, only take the absolute value (see arXiv paper). 
            % Index "Num" for Numeric
            aNum = [0.6555, 0.081837, 0.2245,  0.016932, 0.057375, -1,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
                    0.9947, 0.069636, 0.10291, 0.17456,  0.47666,   0,     -1,      0,      0.13486,0, 0, 0, 0, 0, 0, 0,    0;
                    1.7651, 0.19133,  0.64716, 0.024406, 0.44695,   0,      0,     -1,      0.1621, 0, 0, 0, 0, 0, 0, 0,    0;
                    26.5447,6.7367,  18.4808,  0.15056,  0.4778,    0,      0,      0,      0,      1, 0, 0, 0, 0, 0, 0,    0; 
                    0,      0,        0,       0,        0,         0.18,   0.77,   0.05,  -1,      0, 0, 0, 0, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0.18,   0.77,   0.05,   0,     -1, 0, 0, 0, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0,-1, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0,-1, 0, 0,    0; 
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0,-1, 0,    0;
                    0,     -1,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0, 0, c32,  0; 
                    0,      0,       -1,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0, 0, 0,    c32;]';
        else
            aNum = [0.6555, 0.081837, 0.2245,  0.016932, 0.057375, -1,      0,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
                    0.6555, 0.081837, 0.2245,  0.016932, 0.057375,  0,     -1,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
                    0.9947, 0.069636, 0.10291, 0.17456,  0.47666,   0,      0,     -1,      0,      0.13486,0, 0, 0, 0, 0, 0, 0,    0;
                    1.7651, 0.19133,  0.64716, 0.024406, 0.44695,   0,      0,      0,     -1,      0.1621, 0, 0, 0, 0, 0, 0, 0,    0;
                    26.5447,6.7367,  18.4808,  0.15056,  0.4778,    0,      0,      0,      0,      0,      1, 0, 0, 0, 0, 0, 0,    0; 
                    0,      0,        0,       0,        0,         0.18,   0,      0.77,   0.05,  -1,      0, 0, 0, 0, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0.18,   0,      0.77,   0.05,   0,     -1, 0, 0, 0, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0,-1, 0, 0, 0,    0;
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0,-1, 0, 0,    0; 
                    0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0,-1, 0,    0;
                    0,     -1,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0, 0, c32,  0; 
                    0,      0,       -1,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0, 0, 0,    c32;]';
        end
    
    case 4
        if flagFrac == 0
            aNum = [0.2482,  0.6809,   0.0207,     0.0456,    -1,      0,      0,       0.1372,    0,      0,          0; 
                    0.3221,  0.7954,   0.1689,     0.4588,     0,     -1,      0,       0.1723,    0,      0,          0; 
                    0.6393,  0.5817,   0.0344,     0.4152,     0,      0,     -1,       0.2286,    0,      0,          0; 
                    0,       0,        0,          0,          0.18    0.77,   0.05,   -1,         0,      0,          0;
                   -1,       0,        0,          0,          0,      0,      0,       0,         0,      c22,        0; 
                    0,      -1,        0,          0,          0,      0,      0,       0,         0,      0,          c22]';
        else
            aNum = [0.2482,  0.6809,   0.0207,     0.0456,    -1,      0,      0,      0,       0.1372,    0,      0,          0; 
                    0.2482,  0.6809,   0.0207,     0.0456,     0,     -1,      0,      0,       0.1372,    0,      0,          0;
                    0.3221,  0.7954,   0.1689,     0.4588,     0,      0,     -1,      0,       0.1723,    0,      0,          0; 
                    0.6393,  0.5817,   0.0344,     0.4152,     0,      0,      0,     -1,       0.2286,    0,      0,          0; 
                    0,       0,        0,          0,          0.18,   0,      0.77,   0.05,   -1,         0,      0,          0;
                   -1,       0,        0,          0,          0,      0,      0,      0,       0,         0,      c22,        0; 
                    0,      -1,        0,          0,          0,      0,      0,      0,       0,         0,      0,          c22]';
        end

end


% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    
params.a = aNum; 
params.c = cNum; 
params.th = thNum; 

end