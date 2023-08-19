function f = BMR4_AB_ode_sym
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB (incl. balancing of ash) with 
% constant values for u and xin 

% 11 states
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
% x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11];

% 1 input (feed flow rate)
syms u;

% 4 unknown parameters (rate constants)
syms th1 th2 th3 th4
% p = [th1; th2; th3; th4];

% 21 known & constant parameters...: 
syms c1 c3 c4 c5 c15 c16 c17 c18 c19 c20 c21
% ... + 18 stoichiometric constants:
syms a11 a13 a21 a22 a23 a31 a32 a33 a41 a42 a43 a54 a64 a74 a81 a82 a83
% 9 (assumed) known inlet concentrations: 
syms xi1 xi2 xi3 xi4 xi5 xi6 xi7 xi8 xi9

% dynamic equations
f = [c1*(xi1 - x1)*u + a11*th1*x5 + a21*th2*x6 + a13*th3*x7 - c3*x1 + c4*x10; 
     c1*(xi2 - x2)*u + a21*th1*x5 + a22*th2*x6 + a23*th3*x7 - c3*x2 + c5*x11;
     c1*(xi3 - x3)*u - a31*th1*x5 + a32*th2*x6 - a33*th3*x7; 
     c1*(xi4 - x4)*u - a41*th1*x5 - a42*th2*x6 - a43*th3*x7;
     c1*(xi5 - x5)*u - th1*x5 + a54*th4*x8; 
     c1*(xi6 - x6)*u - th2*x6 + a64*th4*x8;
     c1*(xi7 - x7)*u - th3*x7 + a74*th4*x8;
     c1*(xi8 - x8)*u + a81*th1*x5 + a82*th2*x6 + a83*th3*x7 - th4*x8;
     c1*(xi9 - x9)*u; 
     c15*x10^3 + c16*x10^2*x11 + c17*x10*x11^2 + c18*x10^2 + c19*x10*x11 + c20*x10 + c5*x1; 
     c17*x11^3 + c16*x10*x11^2 + c15*x10^2*x11 + c19*x11^2 + c18*x10*x11 + c21*x11 + c5*x2];
end 
















