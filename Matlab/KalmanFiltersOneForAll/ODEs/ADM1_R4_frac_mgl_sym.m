function g = ADM1_R4_frac_mgl_sym(x,c)
% delivers a symbolic expression (h) of measurement equation of the
% ADM1_R4_frac

% measurement equations
g = [c(6)*x(11)^2 + c(7)*x(11)*x(12) + c(8)*x(12)^2 + c(9)*x(11) + c(10)*x(12) + c(11); % volFlow [l/d]
     c(12)*x(11);
     c(13)*x(12);
     x(3);
     1 - x(4)/c(14);
     1 - x(10)/(c(14)-x(4))];

end 