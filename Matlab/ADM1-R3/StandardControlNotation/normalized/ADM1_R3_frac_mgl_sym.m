%% Version
% (R2022b) Update 5
% Erstelldatum: 10.5.2023
% Autor: Simon Hellmann

function g = ADM1_R3_frac_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the 
% ADM1-R3-frac

% ion balance: 
Phi = -x(13) + (x(4) - x(16))/17 - x(15)/44 - x(14)/60;     % Achtung: VZ von x13 mit Sören klären!
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 

% measurement equations
g = [c(13)*x(17)^2 + c(14)*x(17)*x(18) + c(15)*x(18)^2 + c(16)*x(17) + c(17)*x(18) + c(18);
     c(19)*x(17);
     c(20)*x(18);
     -log10(SHPlus); 
%      x(4);                  % SIN
     x(4) - x(16);          % Snh4
     1 - x(5)/c(21);
     1 - x(12)/(c(21)-x(5)); 
     x(1)];

end 