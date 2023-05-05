%% Version
% (R2022b) Update 5
% Erstelldatum: 21.4.2023
% Autor: Simon Hellmann

function f = BMR4_AB_frac_norm_ode_sym(xNorm, uNorm, xiNorm, th, c, a, Tx, u0)
% delivers a symbolic expression (f) of the normalized right-hand side of the 
% homogeneous ODE system of the BMR4+AB-frac (incl. balancing of ash and 2
% CH fractions) with constant values for u and xin 

% % 12 states
% x = sym('x', [1,12]); 
% 
% % 1 input (feed flow rate)
% syms u;
% 
% % 6 unknown parameters (rate constants)
% th = sym('th', [1,6]); 
% 
% % known & constant parameters...: 
% c = sym('c', [1,21]); 
% % ... + 18 stoichiometric constants:
% syms a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 a55 a75 a85 a91 a92 a93 a94
% % 10 (assumed) known inlet concentrations: 
% xi = sym('xi', [1,10]);  
% 
% % define petersen matrix symbolically:
% a = [a11, a21, a31, a41, -1, 0, 0, 0, a91; 
%      a12, a22, a32, a42, 0, -1, 0, 0, a92; 
%      a13, a23, a33, a43, 0, 0, -1, 0, a93; 
%      a14, a24, a34, a44, 0, 0, 0, -1, a94;
%      0,     0,  0,  0, a55, 0 a75,a85, -1].'; % transpose is necessary for correct indexing in definition of f

% dynamic equations
f = [c(1)*(xiNorm(1) - xNorm(1))*u0*uNorm + a(1,1)*th(1)*Tx(5)/Tx(1)*xNorm(5) + a(1,2)*th(2)*Tx(6)/Tx(1)*xNorm(6) + a(1,3)*th(3)*Tx(7)/Tx(1)*xNorm(7) + a(1,4)*th(4)*Tx(8)/Tx(1)*xNorm(8) - c(2)*xNorm(1) + c(3)*Tx(11)/Tx(1)*xNorm(11);     
     c(1)*(xiNorm(2) - xNorm(2))*u0*uNorm + a(2,1)*th(1)*Tx(5)/Tx(2)*xNorm(5) + a(2,2)*th(2)*Tx(6)/Tx(2)*xNorm(6) + a(2,3)*th(3)*Tx(7)/Tx(2)*xNorm(7) + a(2,4)*th(4)*Tx(8)/Tx(2)*xNorm(8) - c(2)*xNorm(2) + c(4)*Tx(12)/Tx(2)*xNorm(12);    
     c(1)*(xiNorm(3) - xNorm(3))*u0*uNorm - a(3,1)*th(1)*Tx(5)/Tx(3)*xNorm(5) - a(3,2)*th(2)*Tx(6)/Tx(3)*xNorm(6) + a(3,3)*th(3)*Tx(7)/Tx(3)*xNorm(7) - a(3,4)*th(4)*Tx(8)/Tx(3)*xNorm(8);     
     c(1)*(xiNorm(4) - xNorm(4))*u0*uNorm - a(4,1)*th(1)*Tx(5)/Tx(4)*xNorm(5) - a(4,2)*th(2)*Tx(6)/Tx(4)*xNorm(6) - a(4,3)*th(3)*Tx(7)/Tx(4)*xNorm(7) - a(4,4)*th(4)*Tx(8)/Tx(4)*xNorm(8);    
     c(1)*(th(6)*xiNorm(5) - xNorm(5))*u0*uNorm - th(1)*xNorm(5) + a(5,5)*th(5)*Tx(9)/Tx(5)*xNorm(9);     
     c(1)*((1-th(6))*Tx(5)/Tx(6)*xiNorm(5) - xNorm(6))*u0*uNorm - th(2)*xNorm(6);    
     c(1)*(xiNorm(7) - xNorm(7))*u0*uNorm - th(3)*xNorm(7) + a(7,5)*th(5)*Tx(9)/Tx(7)*xNorm(9);     
     c(1)*(xiNorm(8) - xNorm(8))*u0*uNorm - th(4)*xNorm(8) + a(8,5)*th(5)*Tx(9)/Tx(8)*xNorm(9);     
     c(1)*(xiNorm(9) - xNorm(9))*u0*uNorm + a(9,1)*th(1)*Tx(5)/Tx(9)*xNorm(5) + a(9,2)*th(2)*Tx(6)/Tx(9)*xNorm(6) + a(9,3)*th(3)*Tx(7)/Tx(9)*xNorm(7) + a(9,4)*th(4)*Tx(8)/Tx(9)*xNorm(8) - th(5)*xNorm(9);    
     c(1)*(xiNorm(10) - xNorm(10))*u0*uNorm;    
     c(15)*Tx(11)^2*xNorm(11)^3 + c(16)*Tx(11)*Tx(12)*xNorm(11)^2*xNorm(12) + c(17)*Tx(12)^2*xNorm(11)*xNorm(12)^2 + c(18)*Tx(11)*xNorm(11)^2 + c(19)*Tx(12)*xNorm(11)*xNorm(12) + c(20)*xNorm(11) + c(5)*Tx(1)/Tx(11)*xNorm(1);    
     c(17)*Tx(12)^2*xNorm(12)^3 + c(16)*Tx(11)*Tx(12)*xNorm(11)*xNorm(12)^2 + c(15)*Tx(11)^2*xNorm(11)^2*xNorm(12) + c(19)*Tx(12)*xNorm(12)^2 + c(18)*Tx(11)*xNorm(11)*xNorm(12) + c(21)*xNorm(12) + c(5)*Tx(2)/Tx(12)*xNorm(2)];
end 
