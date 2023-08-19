function f = BMR4_AB_ess_ode_sym
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB-ess with constant values for u and xin 

% 8 states
syms x1 x2 x3 x4 x5 x6 x7 x8
% x = [x1; x2; x3; x4; x5; x6; x7; x8];

% 1 input (feed flow rate)
syms u;

% 6 unknown parameters (rate constants)
syms th1 th2 th3 th4
% p = [th1; th2; th3; th4];

% known & constant parameters...: 
syms c1 c2 c3 c4 c5 c15 c16 c17 c18 c19 c20 c21 
% ... + stoichiometric constants:
syms a11 a12 a13 a21 a22 a23 a54 a64 a74 a81 a82 a83
% 6 (assumed) known inlet concentrations: 
syms xi1 xi2 xi3 xi4 xi5 xi6 

% dynamic equations
f = [c1*(xi1 - x1)*u + a11*th1*x3 + a12*th2*x4  + a13*th3*x5 - c2*x1 + c3*x7; 
     c1*(xi2 - x2)*u + a21*th1*x3 + a22*th2*x4  + a23*th3*x5 - c2*x2 + c4*x8;
     c1*(xi3 - x3)*u - th1*x3 + a54*th4*x6; 
     c1*(xi4 - x4)*u - th2*x4 + a64*th4*x6;
     c1*(xi5 - x5)*u - th3*x5 + a74*th4*x6;
     c1*(xi6 - x6)*u + a81*th1*x3 + a82*th2*x4 + a83*th3*x5 - th4*x6;
     c15*x7^3 + c16*x7^2*x8 + c17*x7*x8^2 + c18*x7^2 + c19*x7*x8 + c20*x7 + c5*x1; 
     c17*x8^3 + c16*x7*x8^2 + c15*x7^2*x8 + c19*x8^2 + c18*x7*x8 + c21*x8 + c5*x2];

end 
















