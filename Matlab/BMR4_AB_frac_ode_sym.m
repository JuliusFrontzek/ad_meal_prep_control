function f = BMR4_AB_frac_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the 
% homogeneous ODE system of the BMR4+AB-frac (incl. balancing of ash and 2
% CH fractions) with constant values for u and xin 

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(5) + a(1,2)*th(2)*x(6) + a(1,3)*th(3)*x(7) + a(1,4)*th(4)*x(8) - c(2)*x(1) + c(3)*x(11);     
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) + a(2,4)*th(4)*x(8) - c(2)*x(2) + c(4)*x(12);    
     c(1)*(xi(3) - x(3))*u - a(3,1)*th(1)*x(5) - a(3,2)*th(2)*x(6) + a(3,3)*th(3)*x(7) - a(3,4)*th(4)*x(8);     
     c(1)*(xi(4) - x(4))*u - a(4,1)*th(1)*x(5) - a(4,2)*th(2)*x(6) - a(4,3)*th(3)*x(7) - a(4,4)*th(4)*x(8);    
     c(1)*(th(6)*xi(5) - x(5))*u - th(1)*x(5) + a(5,5)*th(5)*x(9);     
     c(1)*((1-th(6))*xi(5) - x(6))*u - th(2)*x(6);    
     c(1)*(xi(7) - x(7))*u - th(3)*x(7) + a(7,5)*th(5)*x(9);     
     c(1)*(xi(8) - x(8))*u - th(4)*x(8) + a(8,5)*th(5)*x(9);     
     c(1)*(xi(9) - x(9))*u + a(9,1)*th(1)*x(5) + a(9,2)*th(2)*x(6) + a(9,3)*th(3)*x(7) + a(9,4)*th(4)*x(8) - th(5)*x(9);    
     c(1)*(xi(10) - x(10))*u;    
     c(15)*x(11)^3 + c(16)*x(11)^2*x(12) + c(17)*x(11)*x(12)^2 + c(18)*x(11)^2 + c(19)*x(11)*x(12) + c(20)*x(11) + c(5)*x(1);    
     c(17)*x(12)^3 + c(16)*x(11)*x(12)^2 + c(15)*x(11)^2*x(12) + c(19)*x(12)^2 + c(18)*x(11)*x(12) + c(21)*x(12) + c(5)*x(2)];
end 
