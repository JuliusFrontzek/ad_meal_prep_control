%% Version
% (R2022b) Update 5
% Erstelldatum: 27.4.2023
% Autor: Simon Hellmann

function f = ADM1_R3_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the ADM1-R3 (incl. balancing of ash) with 
% constant values for u and xin 

% ion balance: 
Phi = x(12) + (x(4) - x(15))/17 - x(14)/44 - x(13)/60;
% equivalent proton concentration: 
SHPlus = -Phi/2 + 0.5*sqrt(Phi^2 + c(4)); 
% overall inhibition factor: 
Iac = c(3)/(c(3) + SHPlus^(c(2))) * x(4)/(x(4) + c(8)) * th(7)/(th(7) + x(15)); 

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(6) + a(1,2)*th(2)*x(7) + a(1,3)*th(3)*x(8) - a(1,4)*th(5)*x(1)*x(10)/(th(6) + x(1))*Iac; 
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(6) + a(2,2)*th(2)*x(7) + a(2,3)*th(3)*x(8) - c(5)*x(2) + c(6)*x(16) + a(2,4)*th(5)*x(1)*x(10)/(th(6) + x(1))*Iac;
     c(1)*(xi(3) - x(3))*u + a(3,1)*th(1)*x(6) + a(3,2)*th(2)*x(7) - a(3,3)*th(3)*x(8) - c(5)*x(3) + c(5)*x(14) + c(7)*x(17) + a(3,4)*th(5)*x(1)*x(10)/(th(6) + x(1))*Iac;
     c(1)*(xi(4) - x(4))*u - a(4,1)*th(1)*x(6) + a(4,2)*th(2)*x(7) - a(4,3)*th(3)*x(8) - a(4,4)*th(5)*x(1)*x(10)/(th(6) + x(1))*Iac;    
     c(1)*(xi(5) - x(5))*u - a(5,1)*th(1)*x(6) - a(5,2)*th(2)*x(7) - a(5,3)*th(3)*x(8) + a(5,4)*th(5)*x(1)*x(10)/(th(6) + x(1))*Iac;    
     c(1)*(xi(6) - x(6))*u - th(1)*x(6) + a(6,5)*th(4)*x(9) + a(6,6)*th(4)*x(10);    
     c(1)*(xi(7) - x(7))*u - th(2)*x(7) + a(7,5)*th(4)*x(9) + a(7,6)*th(4)*x(10);
     c(1)*(xi(8) - x(8))*u - th(3)*x(8) + a(8,5)*th(4)*x(9) + a(8,6)*th(4)*x(10);
     c(1)*(xi(9) - x(9))*u + a(9,1)*th(1)*x(6) + a(9,2)*th(2)*x(7) + a(9,3)*th(3)*x(8) - th(4)*x(9);    
     c(1)*(xi(10) - x(10))*u - th(4)*x(10) + th(5)*x(1)*x(10)/(th(6) + x(1))*Iac;    
     c(1)*(xi(11) - x(11))*u;
     c(1)*(xi(12) - x(12))*u;
     c(29)*(x(1) - x(13)) - c(9)*x(13)*SHPlus;
     c(30)*(x(3) - x(14)) - c(10)*x(14)*SHPlus;
     c(31)*(x(4) - x(15)) - c(11)*x(15)*SHPlus;
     c(22)*x(16)^3 + c(23)*x(16)^2*x(17) + c(24)*x(16)*x(17)^2 + c(25)*x(16)^2 + c(26)*x(16)*x(17) + c(12)*x(2) + c(27)*x(16);    
     c(24)*x(17)^3 + c(23)*x(16)*x(17)^2 + c(22)*x(16)^2*x(17) + c(26)*x(17)^2 + c(25)*x(16)*x(17) + c(12)*x(3) - c(12)*x(14) + c(28)*x(17)];    
end 


