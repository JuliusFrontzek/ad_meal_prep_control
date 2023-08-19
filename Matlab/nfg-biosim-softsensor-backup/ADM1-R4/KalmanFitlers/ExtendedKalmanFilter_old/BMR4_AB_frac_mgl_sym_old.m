function g = BMR4_AB_frac_mgl_sym_old
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_frac

% 12 states
% x = sym('x',[1 12]);
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12

% known & constant parameters: 
syms c9 c10 c11 c12 c13 c14 c15 c16 c17
% c = sym('c',[1 17]);

% measurement equations
g = [(c6*x11^2 + c7*x11*x12 + c8*x12^2 + c9*x11 + c10*x12 + c11)/24;
     c12*x11; 
     c13*x12; 
     x3; 
     1 - x4/c14; 
     1 - x10/(c14-x4)];

end 