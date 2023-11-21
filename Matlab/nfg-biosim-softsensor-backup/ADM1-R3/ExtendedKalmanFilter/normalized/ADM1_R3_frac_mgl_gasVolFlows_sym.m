%% Version
% (R2022b) Update 5
% Erstelldatum: 03.11.2023
% Autor: Simon Hellmann

function g = ADM1_R3_frac_mgl_gasVolFlows_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the 
% ADM1-R3-frac with additional terms for V_ch4_gas and V_co2_gas as
% extension

% ion balance: 
Phi = x(13) + (x(4) - x(16))/17 - x(15)/44 - x(14)/60;
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 

% measurement equations:
gasVolFlow = c(13)*x(17)^2 + c(14)*x(17)*x(18) + c(15)*x(18)^2 + c(16)*x(17) + c(17)*x(18) + c(18); 
p_ch4 = c(19)*x(17); 
p_co2 = c(20)*x(18); 
pH = -log10(SHPlus); 
S_IN = x(4);
TS = 1 - x(5)/c(21); 
VS = 1 - x(12)/(c(21)-x(5)); 
S_ac = x(1); 

% extensions: 
x_ch4 = p_ch4/(p_ch4 + p_co2); % [-]
x_co2 = p_co2/(p_ch4 + p_co2); % [-]
volFlow_ch4 = gasVolFlow*x_ch4;% [m³/d]
volFlow_co2 = gasVolFlow*x_co2;% [m³/d] 

% measurement equations
g = [gasVolFlow;
     p_ch4;
     p_co2;
     volFlow_ch4; 
     volFlow_co2;
     pH; 
     S_IN;
     TS;
     VS; 
     S_ac];

end 