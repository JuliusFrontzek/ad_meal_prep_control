%% Version
% (R2022b) Update 5
% Erstelldatum: 11.5.2023
% Autor: Simon Hellmann

function f = ADM1_R3_frac_norm_ode_sym(xNorm, uNorm, xiNorm, th, c, a, Tx, Tu)
% delivers a symbolic expression (f) of the normalized right-hand side of 
% the homogeneous ODE system of the ADM1-R3-frac (incl. balancing of ash 
% and two CH-fractions) with constant values for u and xin 

% ion balance: 
PhiNorm = Tx(13)*xNorm(13) + (Tx(4)*xNorm(4) - Tx(16)*xNorm(16))/17 - Tx(15)*xNorm(15)/44 - Tx(14)*xNorm(14)/60;
% equivalent proton concentration: 
SHPlusNorm = -PhiNorm/2 + 0.5*sqrt(PhiNorm^2 + c(4)); 
% overall inhibition factor: 
IacNorm = c(3)/(c(3) + SHPlusNorm^(c(2))) * xNorm(4)/(xNorm(4) + c(8)/Tx(4)) * th(8)/(th(8) + Tx(16)*xNorm(16)); 

% dynamic equations
f = [c(1)*(xiNorm(1) - xNorm(1))*Tu*uNorm + a(1,1)*th(1)*Tx(6)/Tx(1)*xNorm(6) + a(1,2)*th(2)*Tx(7)/Tx(1)*xNorm(7) + a(1,3)*th(3)*Tx(8)/Tx(1)*xNorm(8) + a(1,4)*th(4)*Tx(9)/Tx(1)*xNorm(9) - a(1,5)*th(6)*Tx(11)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm; 
     c(1)*(xiNorm(2) - xNorm(2))*Tu*uNorm + a(2,1)*th(1)*Tx(6)/Tx(2)*xNorm(6) + a(2,2)*th(2)*Tx(7)/Tx(2)*xNorm(7) + a(2,3)*th(3)*Tx(8)/Tx(2)*xNorm(8) + a(2,4)*th(4)*Tx(9)/Tx(2)*xNorm(9) + a(2,5)*th(6)*Tx(1)*Tx(11)/Tx(2)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm - c(5)*xNorm(2) + c(6)*Tx(17)/Tx(2)*xNorm(17);
     c(1)*(xiNorm(3) - xNorm(3))*Tu*uNorm + a(3,1)*th(1)*Tx(6)/Tx(3)*xNorm(6) + a(3,2)*th(2)*Tx(7)/Tx(3)*xNorm(7) + a(3,3)*th(3)*Tx(8)/Tx(3)*xNorm(8) - a(3,4)*th(4)*Tx(9)/Tx(3)*xNorm(9) + a(3,5)*th(6)*Tx(1)*Tx(11)/Tx(3)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm - c(5)*xNorm(3) + c(5)*Tx(15)/Tx(3)*xNorm(15) + c(7)*Tx(18)/Tx(3)*xNorm(18);
     c(1)*(xiNorm(4) - xNorm(4))*Tu*uNorm - a(4,1)*th(1)*Tx(6)/Tx(4)*xNorm(6) - a(4,2)*th(2)*Tx(7)/Tx(4)*xNorm(7) + a(4,3)*th(3)*Tx(8)/Tx(4)*xNorm(8) - a(4,4)*th(4)*Tx(9)/Tx(4)*xNorm(9) - a(4,5)*th(6)*Tx(1)*Tx(11)/Tx(4)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm;    
     c(1)*(xiNorm(5) - xNorm(5))*Tu*uNorm - a(5,1)*th(1)*Tx(6)/Tx(5)*xNorm(6) - a(5,2)*th(2)*Tx(7)/Tx(5)*xNorm(7) - a(5,3)*th(3)*Tx(8)/Tx(5)*xNorm(8) - a(5,4)*th(4)*Tx(9)/Tx(5)*xNorm(9) + a(5,5)*th(6)*Tx(1)*Tx(11)/Tx(5)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm;
     c(1)*(th(9)*xiNorm(6) - xNorm(6))*Tu*uNorm - th(1)*xNorm(6) + a(6,6)*th(5)*Tx(10)/Tx(6)*xNorm(10) + a(6,7)*th(5)*Tx(11)/Tx(6)*xNorm(11);
     c(1)*((1-th(9))*Tx(6)/Tx(7)*xiNorm(6) - xNorm(7))*Tu*uNorm - th(2)*xNorm(7);
     c(1)*(xiNorm(8) - xNorm(8))*Tu*uNorm - th(3)*xNorm(8) + a(8,6)*th(5)*Tx(10)/Tx(8)*xNorm(10) + a(8,7)*th(5)*Tx(11)/Tx(8)*xNorm(11);
     c(1)*(xiNorm(9) - xNorm(9))*Tu*uNorm - th(4)*xNorm(9) + a(9,6)*th(5)*Tx(10)/Tx(9)*xNorm(10) + a(9,7)*th(5)*Tx(11)/Tx(9)*xNorm(11);
     c(1)*(xiNorm(10) - xNorm(10))*Tu*uNorm + a(10,1)*th(1)*Tx(6)/Tx(10)*xNorm(6) + a(10,2)*th(2)*Tx(7)/Tx(10)*xNorm(7) + a(10,3)*th(3)*Tx(8)/Tx(10)*xNorm(8) + a(10,4)*th(4)*Tx(9)/Tx(10)*xNorm(9) - th(5)*xNorm(10);
     c(1)*(xiNorm(11) - xNorm(11))*Tu*uNorm + th(6)*Tx(1)*xNorm(1)*xNorm(11)/(th(7) + Tx(1)*xNorm(1))*IacNorm - th(5)*xNorm(11);
     c(1)*(xiNorm(12) - xNorm(12))*Tu*uNorm;
     c(1)*(xiNorm(13) - xNorm(13))*Tu*uNorm;
     c(29)*(Tx(1)/Tx(14)*xNorm(1) - xNorm(14)) - c(9)*xNorm(14)*SHPlusNorm;
     c(30)*(Tx(3)/Tx(15)*xNorm(3) - xNorm(15)) - c(10)*xNorm(15)*SHPlusNorm;
     c(31)*(Tx(4)/Tx(16)*xNorm(4) - xNorm(16)) - c(11)*xNorm(16)*SHPlusNorm;
     c(22)*Tx(17)^2*xNorm(17)^3 + c(23)*Tx(17)*Tx(18)*xNorm(17)^2*xNorm(18) + c(24)*Tx(18)^2*xNorm(17)*xNorm(18)^2 + c(25)*Tx(17)*xNorm(17)^2 + c(26)*Tx(18)*xNorm(17)*xNorm(18) + c(12)*Tx(2)/Tx(17)*xNorm(2) + c(27)*xNorm(17);
     c(24)*Tx(18)^2*xNorm(18)^3 + c(23)*Tx(17)*Tx(18)*xNorm(17)*xNorm(18)^2 + c(22)*Tx(17)^2*xNorm(17)^2*xNorm(18) + c(26)*Tx(18)*xNorm(18)^2 + c(25)*Tx(17)*xNorm(17)*xNorm(18) + c(12)*Tx(3)/Tx(18)*xNorm(3) - c(12)*Tx(15)/Tx(18)*xNorm(15) + c(28)*xNorm(18)];    
end 


