%% Version
% (R2022b) Update 5
% Erstelldatum: 2.5.2023
% Autor: Simon Hellmann

function g = BMR4_AB_frac_noWater_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_frac_noWater

% measurement equations
g = [c(6)*x(9)^2 + c(7)*x(9)*x(10) + c(8)*x(10)^2 + c(9)*x(9) + c(10)*x(10) + c(11); % [l/d]
     c(12)*x(9);               % [bar]
     c(13)*x(10);               % [bar]
     x(3)];                      % [g/l]
end 