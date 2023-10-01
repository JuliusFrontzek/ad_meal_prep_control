%% Version
% (R2022b) Update 5
% Erstelldatum: 19.09.2023
% Autor: Simon Hellmann

function f = BMR4_AB_frac_noWaterAsh_norm_ode_sym(xNorm, uNorm, xiNorm, th, c, a, Tx, Tu)
% delivers a symbolic expression (f) of the normalized right-hand side of the 
% homogeneous ODE system of the BMR4+AB-frac-noWaterAsh (no water, no ash, 2
% CH fractions) with constant values for u and xin 

% dynamic equations
f = [c(1)*(xiNorm(1) - xNorm(1))*Tu*uNorm + a(1,1)*th(1)*Tx(4)/Tx(1)*xNorm(4) + a(1,2)*th(2)*Tx(5)/Tx(1)*xNorm(5) + a(1,3)*th(3)*Tx(6)/Tx(1)*xNorm(6) + a(1,4)*th(4)*Tx(7)/Tx(1)*xNorm(7) - c(2)*xNorm(1) + c(3)*Tx(9)/Tx(1)*xNorm(9);     
     c(1)*(xiNorm(2) - xNorm(2))*Tu*uNorm + a(2,1)*th(1)*Tx(4)/Tx(2)*xNorm(4) + a(2,2)*th(2)*Tx(5)/Tx(2)*xNorm(5) + a(2,3)*th(3)*Tx(6)/Tx(2)*xNorm(6) + a(2,4)*th(4)*Tx(7)/Tx(2)*xNorm(7) - c(2)*xNorm(2) + c(4)*Tx(10)/Tx(2)*xNorm(10);    
     c(1)*(xiNorm(3) - xNorm(3))*Tu*uNorm - a(3,1)*th(1)*Tx(4)/Tx(3)*xNorm(4) - a(3,2)*th(2)*Tx(5)/Tx(3)*xNorm(5) + a(3,3)*th(3)*Tx(6)/Tx(3)*xNorm(6) - a(3,4)*th(4)*Tx(7)/Tx(3)*xNorm(7);     
%      c(1)*(xiNorm(4) - xNorm(4))*Tu*uNorm - a(4,1)*th(1)*Tx(5)/Tx(4)*xNorm(5) - a(4,2)*th(2)*Tx(6)/Tx(4)*xNorm(6) - a(4,3)*th(3)*Tx(7)/Tx(4)*xNorm(7) - a(4,4)*th(4)*Tx(8)/Tx(4)*xNorm(8);    
     c(1)*(th(6)*xiNorm(4) - xNorm(4))*Tu*uNorm - th(1)*xNorm(4) + a(4,5)*th(5)*Tx(8)/Tx(4)*xNorm(8);     
     c(1)*((1-th(6))*Tx(4)/Tx(5)*xiNorm(4) - xNorm(5))*Tu*uNorm - th(2)*xNorm(5);    
     c(1)*(xiNorm(6) - xNorm(6))*Tu*uNorm - th(3)*xNorm(6) + a(6,5)*th(5)*Tx(8)/Tx(6)*xNorm(8);     
     c(1)*(xiNorm(7) - xNorm(7))*Tu*uNorm - th(4)*xNorm(7) + a(7,5)*th(5)*Tx(8)/Tx(7)*xNorm(8);     
     c(1)*(xiNorm(8) - xNorm(8))*Tu*uNorm + a(8,1)*th(1)*Tx(4)/Tx(8)*xNorm(4) + a(8,2)*th(2)*Tx(5)/Tx(8)*xNorm(5) + a(8,3)*th(3)*Tx(6)/Tx(8)*xNorm(6) + a(8,4)*th(4)*Tx(7)/Tx(8)*xNorm(7) - th(5)*xNorm(8);    
%      c(1)*(xiNorm(10) - xNorm(10))*Tu*uNorm;    
     c(15)*Tx(9)^2*xNorm(9)^3 + c(16)*Tx(9)*Tx(10)*xNorm(9)^2*xNorm(10) + c(17)*Tx(10)^2*xNorm(9)*xNorm(10)^2 + c(18)*Tx(9)*xNorm(9)^2 + c(19)*Tx(10)*xNorm(9)*xNorm(10) + c(20)*xNorm(9) + c(5)*Tx(1)/Tx(9)*xNorm(1);    
     c(17)*Tx(10)^2*xNorm(10)^3 + c(16)*Tx(9)*Tx(10)*xNorm(9)*xNorm(10)^2 + c(15)*Tx(9)^2*xNorm(9)^2*xNorm(10) + c(19)*Tx(10)*xNorm(10)^2 + c(18)*Tx(9)*xNorm(9)*xNorm(10) + c(21)*xNorm(10) + c(5)*Tx(2)/Tx(10)*xNorm(2)];
end 
