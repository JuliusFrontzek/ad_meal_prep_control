%% Version
% (R2022b) Update 5
% Erstelldatum: 2.11.2023
% Autor: Simon Hellmann

function gExt = BMR4_AB_frac_mgl_gasVolFlows_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_frac with additional terms for V_ch4_gas and V_co2_gas as
% extension

% measurement equations:
gasVolFlow = c(6)*x(11)^2 + c(7)*x(11)*x(12) + c(8)*x(12)^2 + c(9)*x(11) + c(10)*x(12) + c(11); % [m³/d]
p_ch4 = c(12)*x(11); % [bar]
p_co2 = c(13)*x(12); % [bar]
S_IN = x(3);         % [kg/m³]
TS = 1 - x(4)/c(14); % [-]
VS = 1 - x(10)/(c(14)-x(4)); % [-]

% extensions:
x_ch4 = p_ch4/(p_ch4 + p_co2); % [-]
x_co2 = p_co2/(p_ch4 + p_co2); % [-]
volFlow_ch4 = gasVolFlow*x_ch4;% [m³/d]
volFlow_co2 = gasVolFlow*x_co2;% [m³/d] 

gExt =  [gasVolFlow; 
         p_ch4;
         p_co2;
         volFlow_ch4; 
         volFlow_co2;
         S_IN;
         TS;
         VS];

end 