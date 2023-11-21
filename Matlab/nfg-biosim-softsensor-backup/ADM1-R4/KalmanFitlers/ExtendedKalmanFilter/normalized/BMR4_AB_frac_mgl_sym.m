%% Version
% (R2022b) Update 5
% Erstelldatum: 2.5.2023
% last updated: 1.11.2023
% Autor: Simon Hellmann

function g = BMR4_AB_frac_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_frac

% measurement equations
g = [c(6)*x(11)^2 + c(7)*x(11)*x(12) + c(8)*x(12)^2 + c(9)*x(11) + c(10)*x(12) + c(11); % [m³/d]
     c(12)*x(11);               % [bar]
     c(13)*x(12);               % [bar]
     x(3);                      % [kg/m³]
     1 - x(4)/c(14);            % [-]
     1 - x(10)/(c(14)-x(4))];   % [-]

end 