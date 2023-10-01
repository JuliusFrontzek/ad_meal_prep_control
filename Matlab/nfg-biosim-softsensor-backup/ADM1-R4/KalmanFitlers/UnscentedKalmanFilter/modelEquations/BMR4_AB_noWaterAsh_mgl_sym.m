%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.2023
% Autor: Simon Hellmann

function g = BMR4_AB_noWaterAsh_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_noWaterAsh

% measurement equations
g = [c(6)*x(8)^2 + c(7)*x(8)*x(9) + c(8)*x(9)^2 + c(9)*x(8) + c(10)*x(9) + c(11); % [l/d]
     c(12)*x(8);               % [bar]
     c(13)*x(9);               % [bar]
     x(3)];                    % [g/l]
end 