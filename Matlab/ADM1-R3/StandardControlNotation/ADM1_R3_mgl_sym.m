%% Version
% (R2022b) Update 5
% Erstelldatum: 27.4.2023
% Autor: Simon Hellmann

function g = ADM1_R3_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the ADM1-R3

% ion balance: 
Phi = x(12) + (x(4) - x(15))/17 - x(14)/44 - x(13)/60; 
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 

% measurement equations
g = [c(13)*x(16)^2 + c(14)*x(16)*x(17) + c(15)*x(17)^2 + c(16)*x(16) + c(17)*x(17) + c(18);
     c(19)*x(16);
     c(20)*x(17);
     -log10(SHPlus); 
     x(4);
     1 - x(5)/c(21);
     1 - x(11)/(c(21)-x(5)); 
     x(1)];

end 