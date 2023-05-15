%% Version
% (R2022b) Update 5
% Erstelldatum: 10.5.2023
% Autor: Simon Hellmann

function f = ADM1_R3_frac_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the ADM1-R3-frac (incl. balancing of ash and 
% two CH-fractions) with constant values for u and xin 

% ion balance: 
Phi = -x(13) + (x(4) - x(16))/17 - x(15)/44 - x(14)/60; % Achtung: VZ von x13 mit Sören klären!
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 
% overall inhibition factor: 
Iac = c(3)/(c(3) + SHPlus^(c(2))) * x(4)/(x(4) + c(8)) * th(8)/(th(8) + x(16)); 

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(6) + a(1,2)*th(2)*x(7) + a(1,3)*th(3)*x(8) + a(1,4)*th(4)*x(9) - a(1,5)*th(6)*x(1)*x(11)/(th(7) + x(1))*Iac; 
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(6) + a(2,2)*th(2)*x(7) + a(2,3)*th(3)*x(8) + a(2,4)*th(4)*x(9) + a(2,5)*th(6)*x(1)*x(11)/(th(7) + x(1))*Iac - c(5)*x(2) + c(6)*x(17);
     c(1)*(xi(3) - x(3))*u + a(3,1)*th(1)*x(6) + a(3,2)*th(2)*x(7) + a(3,3)*th(3)*x(8) - a(3,4)*th(4)*x(9) + a(3,5)*th(6)*x(1)*x(11)/(th(7) + x(1))*Iac - c(5)*x(3) + c(5)*x(15) + c(7)*x(18);
     c(1)*(xi(4) - x(4))*u - a(4,1)*th(1)*x(6) - a(4,2)*th(2)*x(7) + a(4,3)*th(3)*x(8) - a(4,4)*th(4)*x(9) - a(4,5)*th(6)*x(1)*x(11)/(th(7) + x(1))*Iac;    
     c(1)*(xi(5) - x(5))*u - a(5,1)*th(1)*x(6) - a(5,2)*th(2)*x(7) - a(5,3)*th(3)*x(8) - a(5,4)*th(4)*x(9) + a(5,5)*th(6)*x(1)*x(11)/(th(7) + x(1))*Iac;
     c(1)*(th(9)*xi(6) - x(6))*u - th(1)*x(6) + a(6,6)*th(5)*x(10) + a(6,7)*th(5)*x(11);
     c(1)*((1-th(9))*xi(6) - x(7))*u - th(2)*x(7);
     c(1)*(xi(8) - x(8))*u - th(3)*x(8) + a(8,6)*th(5)*x(10) + a(8,7)*th(5)*x(11);
     c(1)*(xi(9) - x(9))*u - th(4)*x(9) + a(9,6)*th(5)*x(10) + a(9,7)*th(5)*x(11);
     c(1)*(xi(10) - x(10))*u + a(10,1)*th(1)*x(6) + a(10,2)*th(2)*x(7) + a(10,3)*th(3)*x(8) + a(10,4)*th(4)*x(9) - th(5)*x(10);    
     c(1)*(xi(11) - x(11))*u + th(6)*x(1)*x(11)/(th(7) + x(1))*Iac - th(5)*x(11);    
     c(1)*(xi(12) - x(12))*u;
     c(1)*(xi(13) - x(13))*u;
     c(29)*(x(1) - x(14)) - c(9)*x(14)*SHPlus;
     c(30)*(x(3) - x(15)) - c(10)*x(15)*SHPlus;
     c(31)*(x(4) - x(16)) - c(11)*x(16)*SHPlus;
     c(22)*x(17)^3 + c(23)*x(17)^2*x(18) + c(24)*x(17)*x(18)^2 + c(25)*x(17)^2 + c(26)*x(17)*x(18) + c(12)*x(2) + c(27)*x(17);    
     c(24)*x(18)^3 + c(23)*x(17)*x(18)^2 + c(22)*x(17)^2*x(18) + c(26)*x(18)^2 + c(25)*x(17)*x(18) + c(12)*x(3) - c(12)*x(15) + c(28)*x(18)];    
end 


