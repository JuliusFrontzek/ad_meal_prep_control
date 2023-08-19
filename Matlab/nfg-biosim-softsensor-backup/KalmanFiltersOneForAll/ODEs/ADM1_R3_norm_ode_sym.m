%% Version
% (R2022b) Update 5
% Erstelldatum: 03.07.2023
% Autor: Simon Hellmann

function f = ADM1_R3_norm_ode_sym(xNorm, uNorm, xiNorm, th, c, a, Tx, Tu)
% delivers a symbolic expression (f) of the normalized right-hand side of 
% the homogeneous ODE system of the ADM1-R3 (incl. balancing of ash) with 
% constant values for u and xin 

% ion balance: 
PhiNorm = Tx(12)*xNorm(12) + (Tx(4)*xNorm(4) - Tx(15)*xNorm(15))/17 - Tx(14)*xNorm(14)/44 - Tx(13)*xNorm(13)/60;
% equivalent proton concentration: 
SHPlusNorm = -PhiNorm/2 + 0.5*sqrt(PhiNorm^2 + c(4)); 
% overall inhibition factor: 
IacNorm = c(3)/(c(3) + SHPlusNorm^(c(2))) * xNorm(4)/(xNorm(4) + c(8)/Tx(4)) * th(8)/(th(8) + Tx(15)*xNorm(15)); 

% dynamic equations
f = [c(1)*(xiNorm(1) - xNorm(1))*Tu*uNorm + a(1,1)*th(1)*Tx(6)/Tx(1)*xNorm(6) + a(1,2)*th(3)*Tx(7)/Tx(1)*xNorm(7) + a(1,3)*th(4)*Tx(8)/Tx(1)*xNorm(8) - a(1,4)*th(6)*Tx(10)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm; 
     c(1)*(xiNorm(2) - xNorm(2))*Tu*uNorm + a(2,1)*th(1)*Tx(6)/Tx(2)*xNorm(6) + a(2,2)*th(3)*Tx(7)/Tx(2)*xNorm(7) + a(2,3)*th(4)*Tx(8)/Tx(2)*xNorm(8) + a(2,4)*th(6)*Tx(1)*Tx(10)/Tx(2)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm - c(5)*xNorm(2) + c(6)*Tx(16)/Tx(2)*xNorm(16);
     c(1)*(xiNorm(3) - xNorm(3))*Tu*uNorm + a(3,1)*th(1)*Tx(6)/Tx(3)*xNorm(6) + a(3,2)*th(3)*Tx(7)/Tx(3)*xNorm(7) - a(3,3)*th(4)*Tx(8)/Tx(3)*xNorm(8) + a(3,4)*th(6)*Tx(1)*Tx(10)/Tx(3)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm - c(5)*xNorm(3) + c(5)*Tx(14)/Tx(3)*xNorm(14) + c(7)*Tx(17)/Tx(3)*xNorm(17);
     c(1)*(xiNorm(4) - xNorm(4))*Tu*uNorm - a(4,1)*th(1)*Tx(6)/Tx(4)*xNorm(6) + a(4,2)*th(3)*Tx(7)/Tx(4)*xNorm(7) - a(4,3)*th(4)*Tx(8)/Tx(4)*xNorm(8) - a(4,4)*th(6)*Tx(1)*Tx(10)/Tx(4)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm;    
     c(1)*(xiNorm(5) - xNorm(5))*Tu*uNorm - a(5,1)*th(1)*Tx(6)/Tx(5)*xNorm(6) - a(5,2)*th(3)*Tx(7)/Tx(5)*xNorm(7) - a(5,3)*th(4)*Tx(8)/Tx(5)*xNorm(8) + a(5,4)*th(6)*Tx(1)*Tx(10)/Tx(5)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm;
     c(1)*(xiNorm(6) - xNorm(6))*Tu*uNorm - th(1)*xNorm(6) + a(6,5)*th(5)*Tx(9)/Tx(6)*xNorm(9) + a(6,6)*th(5)*Tx(10)/Tx(6)*xNorm(10);
     c(1)*(xiNorm(7) - xNorm(7))*Tu*uNorm - th(3)*xNorm(7) + a(7,5)*th(5)*Tx(9)/Tx(7)*xNorm(9) + a(7,6)*th(5)*Tx(10)/Tx(7)*xNorm(10);
     c(1)*(xiNorm(8) - xNorm(8))*Tu*uNorm - th(4)*xNorm(8) + a(8,5)*th(5)*Tx(9)/Tx(8)*xNorm(9) + a(8,6)*th(5)*Tx(10)/Tx(8)*xNorm(10);
     c(1)*(xiNorm(9) - xNorm(9))*Tu*uNorm + a(9,1)*th(1)*Tx(6)/Tx(9)*xNorm(6) + a(9,2)*th(3)*Tx(7)/Tx(9)*xNorm(7) + a(9,3)*th(4)*Tx(8)/Tx(9)*xNorm(8) - th(5)*xNorm(9);
     c(1)*(xiNorm(10) - xNorm(10))*Tu*uNorm + th(6)*Tx(1)*xNorm(1)*xNorm(10)/(th(7) + Tx(1)*xNorm(1))*IacNorm - th(5)*xNorm(10);
     c(1)*(xiNorm(11) - xNorm(11))*Tu*uNorm;
     c(1)*(xiNorm(12) - xNorm(12))*Tu*uNorm;
     c(29)*(Tx(1)/Tx(13)*xNorm(1) - xNorm(13)) - c(9)*xNorm(13)*SHPlusNorm;
     c(30)*(Tx(3)/Tx(14)*xNorm(3) - xNorm(14)) - c(10)*xNorm(14)*SHPlusNorm;
     c(31)*(Tx(4)/Tx(15)*xNorm(4) - xNorm(15)) - c(11)*xNorm(15)*SHPlusNorm;
     c(22)*Tx(16)^2*xNorm(16)^3 + c(23)*Tx(16)*Tx(17)*xNorm(16)^2*xNorm(17) + c(24)*Tx(17)^2*xNorm(16)*xNorm(17)^2 + c(25)*Tx(16)*xNorm(16)^2 + c(26)*Tx(17)*xNorm(16)*xNorm(17) + c(12)*Tx(2)/Tx(16)*xNorm(2) + c(27)*xNorm(16);
     c(24)*Tx(17)^2*xNorm(17)^3 + c(23)*Tx(16)*Tx(17)*xNorm(16)*xNorm(17)^2 + c(22)*Tx(16)^2*xNorm(16)^2*xNorm(17) + c(26)*Tx(17)*xNorm(17)^2 + c(25)*Tx(16)*xNorm(16)*xNorm(17) + c(12)*Tx(3)/Tx(17)*xNorm(3) - c(12)*Tx(14)/Tx(17)*xNorm(14) + c(28)*xNorm(17)];    
end 
