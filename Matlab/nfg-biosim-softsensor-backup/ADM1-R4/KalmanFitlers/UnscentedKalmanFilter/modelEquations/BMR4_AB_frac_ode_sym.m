function f = BMR4_AB_frac_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB-frac (incl. balancing of ash and 2
% CH fractions) with constant values for u and xin 

% % 12 states
% x = sym('x', [1,12]); 
% ["S_ch4","S_IC","S_IN","S_h2o","X_ch_fast","X_ch_slow","X_pr","X_li","X_bac","X_ash", "S_gas_ch4","S_gas_co2"]
% % 1 input (feed flow rate) "q_in"
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
f = [c( 1)*(xi( 1)      - x( 1))*u + a(1,1)*th(1)*x(5) + a(1,2)*th(2)*x(6) + a(1,3)*th(3)*x(7) + a(1,4)*th(4)*x(8) - c(2)*x(1) + c(3)*x(11); % S_ch4
     c( 1)*(xi( 2)      - x( 2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) + a(2,4)*th(4)*x(8) - c(2)*x(2) + c(4)*x(12); % S_IC
     c( 1)*(xi( 3)      - x( 3))*u - a(3,1)*th(1)*x(5) - a(3,2)*th(2)*x(6) + a(3,3)*th(3)*x(7) - a(3,4)*th(4)*x(8);                          % S_IN
     c( 1)*(xi( 4)      - x( 4))*u - a(4,1)*th(1)*x(5) - a(4,2)*th(2)*x(6) - a(4,3)*th(3)*x(7) - a(4,4)*th(4)*x(8);                          % S_h2o
     c( 1)*(th(6)*xi(5) - x( 5))*u -        th(1)*x(5)                                                             + a(5,5)*th(5)*x(9);      % X_ch_fast
     c( 1)*((1-th(6))*xi(5) - x(6))*u                  -        th(2)*x(6);                                                                  % X_ch_slow
     c( 1)*(xi( 7)      - x( 7))*u                                         -        th(3)*x(7)                     + a(7,5)*th(5)*x(9);      % X_pr
     c( 1)*(xi( 8)      - x( 8))*u                                                             -        th(4)*x(8) + a(8,5)*th(5)*x(9);      % X_li
     c( 1)*(xi( 9)      - x( 9))*u + a(9,1)*th(1)*x(5) + a(9,2)*th(2)*x(6) + a(9,3)*th(3)*x(7) + a(9,4)*th(4)*x(8) -        th(5)*x(9);      % X_bac
     c( 1)*(xi(10)      - x(10))*u;                                                                                                          % X_ash
     c(15)*x(11)^3 + c(16)*x(11)^2*x(12) + c(17)*x(11)*x(12)^2 + c(18)*x(11)^2 + c(19)*x(11)*x(12) + c(20)*x(11) + c(5)*x(1);                % S_gas_ch4
     c(17)*x(12)^3 + c(16)*x(11)*x(12)^2 + c(15)*x(11)^2*x(12) + c(19)*x(12)^2 + c(18)*x(11)*x(12) + c(21)*x(12) + c(5)*x(2)];               % S_gas_co2
end 