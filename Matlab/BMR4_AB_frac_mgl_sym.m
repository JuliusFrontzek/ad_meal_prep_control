function g = BMR4_AB_frac_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% BMR4+AB_frac

% % 12 states
% x = sym('x',[1 12]);
% 
% % known & constant parameters: 
% c = sym('c',[1 14]);

% measurement equations
g = [(c(6)*x(11)^2 + c(7)*x(11)*x(12) + c(8)*x(12)^2 + c(9)*x(11) + c(10)*x(12) + c(11))/24;
     c(12)*x(11);
     c(13)*x(12);
     x(3);
     1 - x(4)/c(14);
     1 - x(10)/(c(14)-x(4))];

end 