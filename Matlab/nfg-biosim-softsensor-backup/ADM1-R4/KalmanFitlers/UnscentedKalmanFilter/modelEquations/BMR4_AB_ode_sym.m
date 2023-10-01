%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.23
% Autor: Simon Hellmann

function f = BMR4_AB_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB (incl. balancing of ash) with 
% constant values for u and xin 

% dynamic equations
f = [c( 1)*(xi( 1) - x( 1))*u + a(1,1)*th(1)*x(5) + a(1,2)*th(2)*x(6) + a(1,3)*th(3)*x(7) - c(2)*x(1) + c(3)*x(10); % S_ch4
     c( 1)*(xi( 2) - x( 2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) - c(2)*x(2) + c(4)*x(11); % S_IC
     c( 1)*(xi( 3) - x( 3))*u - a(3,1)*th(1)*x(5) + a(3,2)*th(2)*x(6) - a(3,3)*th(3)*x(7); % S_IN
     c( 1)*(xi( 4) - x( 4))*u - a(4,1)*th(1)*x(5) - a(4,2)*th(2)*x(6) - a(4,3)*th(3)*x(7); % S_h2o
     c( 1)*(xi( 5) - x( 5))*u - th(1)*x(5) + a(5,4)*th(4)*x(8);     % X_ch
     c( 1)*(xi( 6) - x( 6))*u - th(2)*x(6) + a(6,4)*th(4)*x(8);     % X_pr
     c( 1)*(xi( 7) - x( 7))*u - th(3)*x(7) + a(7,4)*th(4)*x(8);     % X_li
     c( 1)*(xi( 8) - x( 8))*u + a(8,1)*th(1)*x(5) + a(8,2)*th(2)*x(6) + a(8,3)*th(3)*x(7) - th(4)*x(8); % X_bac
     c( 1)*(xi( 9) - x( 9))*u;                                      % X_ash
     c(15)*x(10)^3 + c(16)*x(10)^2*x(11) + c(17)*x(10)*x(11)^2 + c(18)*x(10)^2 + c(19)*x(10)*x(11) + c(20)*x(10) + c(5)*x(1);  % S_gas_ch4
     c(17)*x(11)^3 + c(16)*x(10)*x(11)^2 + c(15)*x(10)^2*x(11) + c(19)*x(11)^2 + c(18)*x(10)*x(11) + c(21)*x(11) + c(5)*x(2)]; % S_gas_co2
end 
