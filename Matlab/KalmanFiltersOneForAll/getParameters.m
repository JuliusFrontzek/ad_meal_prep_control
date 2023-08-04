%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  params = getParameters(flagModel, flagFrac, system, parameters)
% XY: Zweck der Funktion und Argumente beschreiben

% time-invariant parameters
s = system.Variables;   % s = systemParameters
% extract system parameter values: 
Vl = s(1);  % gas volume
Vg = s(2);  % liquid volume
% p0 = s(3);  % atmospheric pressure [bar] XY: Achtung: Typo im GitHub
p0 = 1.0133; 

switch flagModel
    case 3
        % renaming for easier understanding: 
%         systemInput = input.ADM1_R3.Variables;
        parameters = parameters.ADM1_R3.Variables; % modelParameters
        
        % extract model parameter values:  
        Kch4 =  parameters(1);  % Henry parameter ch4 [mol/l/bar]
        Kco2 =  parameters(2);  % Henry parameter co2 [mol/l/bar]
        KSIN =  parameters(3);  % half-saturation constant nitrogen limitation [g/l]
        KInh3 =  parameters(4); % ammonia inbibition constant [g/l]
        % KaIN =  parameters(5);  % dissociation constant ammonia [mol/l]
        % Kaac =  parameters(6);  % dissociation constant acetic acid [mol/l]
        % Kaco2 =  parameters(7); % dissociation constant carbonate [mol/l]
        % corrected values acc. to Simons computation (Sören made slight mistake):
        KaIN = 1.3809E-9; 
        Kaac = 1.7378E-5; 
        Kaco2 = 5.0981E-7;
        KSac =  parameters(8);  % half-saturation constant acetoclastic methanogenesis [g/l]
        Kw =  parameters(9);    % ion product of water [mol/l]
        R =  parameters(10);    % universal gas constant [bar l/mol/K]
        T =  parameters(11);    % operating temperature [K]
        kABIN =  parameters(12);    % kin. dissociation constant ammonia [l/mol/d] 
        kABac =  parameters(13);    % kin. dissociation constant acetic acid [l/mol/d]
        kABco2 =  parameters(14);   % kin. dissociation constant carbonate [l/mol/d]
        kLa =  parameters(15);  % mass transfer coefficient [1/d]
        kch =  parameters(16);  % hydrolysis constant carbohydrates [1/d]
        kdec =  parameters(17); % hydrolysis constant biomass decay [1/d]
        kli =  parameters(18);  % hydrolysis constant lipids [1/d]
        muM =  parameters(19);  % max. growth rate [1/d]
        kp =  parameters(20);   % friction parameter [l/bar/d]
        kpr =  parameters(21);  % hydrolysis constant proteins [1/d]
        pHLL =  parameters(22); % lower pH boundary  
        pHUL =  parameters(23); % upper pH boundary
        ph2o =  parameters(24); % partial pressure of water in gas phase (saturated) [bar]
        
        rho = 1000;        % mass density of digestate [g/l]
        Mch4 = 16;      % molar mass CH4 [kg/kmol]
        Mco2 = 44;      % molar mass CO2 [kg/kmol]
        nac = 3/(pHUL - pHLL); 
        
        % order model parameters in the rights structures (prepare simulation)
        c1 = 1/Vl; 
        c2 = nac; 
        c3 = 10^(-(3/2)*(pHUL + pHLL)/(pHUL - pHLL)); 
        c4 = 4*Kw; 
        c5 = kLa; 
        c6 = kLa*Kch4*R*T; 
        c7 = kLa*Kco2*R*T; 
        c8 = KSIN; 
        c9 = kABac;
        c10 = kABco2; 
        c11 = kABIN; 
        c12 = kLa*Vl/Vg; 
        c13 = kp/p0*(R*T/Mch4)^2;
        c14 = 2*kp/p0*(R*T)^2/Mch4/Mco2;
        c15 = kp/p0*(R*T/Mco2)^2;
        c16 = kp/p0*R*T/Mch4*(2*ph2o - p0); 
        c17 = kp/p0*R*T/Mco2*(2*ph2o - p0); 
        c18 = kp/p0*(ph2o - p0)*ph2o; 
        c19 = R*T/Mch4;
        c20 = R*T/Mco2;
        c21 = rho; 
        c22 = -kp/Vg/p0*(R*T/Mch4)^2;
        c23 = -2*kp/Vg/p0*(R*T)^2/Mch4/Mco2;
        c24 = -kp/Vg/p0*(R*T/Mco2)^2;
        c25 = -kp/Vg/p0*(R*T/Mch4)*(2*ph2o - p0);
        c26 = -kp/Vg/p0*(R*T/Mco2)*(2*ph2o - p0);
        c27 = -kLa*Vl/Vg*Kch4*R*T - kp/Vg/p0*(ph2o - p0)*ph2o;
        c28 = -kLa*Vl/Vg*Kco2*R*T - kp/Vg/p0*(ph2o - p0)*ph2o;
        c29 = kABac*Kaac; 
        c30 = kABco2*Kaco2;
        c31 = kABIN*KaIN; 
        c32 = Vl/Vg;
        c33 = -kp/Vg/p0*(ph2o - p0)*ph2o;

        % combine all in column vector: 
        cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,...
                c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33]';

end

% time-variant parameters
switch flagModel
    case 3
        if flagFrac == 0
            thNum = [kch, 0, kpr, kli, kdec, muM, KSac, KInh3, 0]';
        else 
            % Fall R3-frac nachtragen
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
            aNum = nan;     % XY: hier noch richtige Werte eintragen
        end
    
    case 4
        % XY: hier noch richtige Werte eintragen

end


% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    
params.a = aNum; 
params.c = cNum; 
params.th = thNum; 

end