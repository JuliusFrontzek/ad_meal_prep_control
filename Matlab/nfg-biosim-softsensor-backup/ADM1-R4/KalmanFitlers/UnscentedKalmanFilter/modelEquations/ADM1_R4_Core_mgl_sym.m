%% Version
% (R2022b) Update 5
% Erstelldatum: 06.10.2023
% Autor: Simon Hellmann

function g = ADM1_R4_Core_mgl_sym(x)
% delivers a symbolic expression (h) of measurement equation of the
% ADM1-R4-Core with bold assumptions

% assume you can measure dissolved CH4 (S_ch4)/CO2 (S_IC) and biomass (bac)

% measurement equations
g = [x(1);  % dissolved CH4 (S_ch4)
     x(2);  % dissolved CO2 (S_IC)
     x(6)]; % biomass bac
end 