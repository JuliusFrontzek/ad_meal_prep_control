%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.2023
% Autor: Simon Hellmann

function g = BMR4_AB_mgl_bac_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB

% measurement equations
g = [c(6)*x(10)^2 + c(7)*x(10)*x(11) + c(8)*x(11)^2 + c(9)*x(10) + c(10)*x(11) + c(11); % [l/d] or [m³/d], depending on parameters
     c(12)*x(10);               % pch4 [bar]
     c(13)*x(11);               % pco2 [bar]
     x(3);                      % SIN [g/l] = [kg/m³]
     1 - x(4)/c(14);            % TS [-]
     1 - x(9)/(c(14)-x(4));     % VS [-]
     x(8)];                     % biomass

end 