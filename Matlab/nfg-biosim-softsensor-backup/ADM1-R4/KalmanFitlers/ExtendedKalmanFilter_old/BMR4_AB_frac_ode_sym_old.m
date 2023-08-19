function f = BMR4_AB_frac_ode_sym_old
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB-frac (incl. balancing of ash and 2
% CH fractions) with constant values for u and xin 

% 12 states
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
% x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

% 1 input (feed flow rate)
syms u;

% 6 unknown parameters (rate constants)
syms t1 t2 t3 t4 t5 t6
% p = [t1; t2; t3; t4];

% known & constant parameters...: 
syms c1 c2 c3 c4 c5 c15 c16 c17 c18 c19 c20 c21
% ... + 18 stoichiometric constants:
syms a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 a55 a75 a85 a91 a92 a93 a94
% 10 (assumed) known inlet concentrations: 
syms xi1 xi2 xi3 xi4 xi5 xi7 xi8 xi9 xi10

% dynamic equations
f = [c1*(xi1 - x1)*u + a11*t1*x5 + a12*t2*x6 + a13*t3*x7 + a14*t4*x8 - c2*x1 + c3*x11; 
     c1*(xi2 - x2)*u + a21*t1*x5 + a22*t2*x6 + a23*t3*x7 + a24*t4*x8 - c2*x2 + c4*x12;
     c1*(xi3 - x3)*u - a31*t1*x5 - a32*t2*x6 + a33*t3*x7 - a34*t4*x8; 
     c1*(xi4 - x4)*u - a41*t1*x5 - a42*t2*x6 - a43*t3*x7 - a44*t4*x8;
     c1*(t6*xi5 - x5)*u - t1*x5 + a55*t5*x9; 
     c1*((1-t6)*xi5 - x6)*u - t2*x6;
     c1*(xi7 - x7)*u - t3*x7 + a75*t5*x9;
     c1*(xi8 - x8)*u - t4*x8 + a85*t5*x9; 
     c1*(xi9 - x9)*u + a91*t1*x5 + a92*t2*x6 + a93*t3*x7 + a94*t4*x8 - t5*x9; 
     c1*(xi10 - x10)*u;
     c15*x11^3 + c16*x11^2*x12 + c17*x11*x12^2 + c18*x11^2 + c19*x11*x12 + c20*x11 + c5*x1; 
     c17*x12^3 + c16*x11*x12^2 + c15*x11^2*x12 + c19*x12^2 + c18*x11*x12 + c21*x12 + c5*x2];

end 
















