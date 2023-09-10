%% Version
% (R2022b) Update 5
% Erstelldatum: 03.09.2023
% Autor: Simon Hellmann

function g = ADM1_R4_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% ADM1-R4

% measurement equations
g = [c(6)*x(10)^2 + c(7)*x(10)*x(11) + c(8)*x(11)^2 + c(9)*x(10) + c(10)*x(11) + c(11); % volFlow [l/d]
     c(12)*x(10);
     c(13)*x(11);
     x(3);
     1 - x(4)/c(14);
     1 - x(9)/(c(14)-x(4))];

end 