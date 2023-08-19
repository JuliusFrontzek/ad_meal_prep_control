function dxdt = f(t,x,feedVector,AC)
% dxdt - right-hand side of ODE 
% t - time vector
% x - states
% feedVector - combination of [feedVolFlow; inlet concentrations]
% AC - struct with stoichiometric coefficients and aggregated constant
% model parameters

% extract constant parameters out of struct: 
th = AC.th; 
c = AC.c; 
a = AC.a;

% hard reset negativer Konzentrationen:
x(x < 0) = 0; 

% Stellgröße und Eingangskonzentrationen:
u = feedVector(2);  % feedVolumeFlow
xi = feedVector(3:end); % Note: für x10, x11 (Gasphase) gibt es prinzipiell 
% keine Eingangskonzentrationen, daher hat xi auch nur 9 Einträge

%% DGLs für Zustände:

f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(5) + a(2,1)*th(2)*x(6) + a(1,3)*th(3)*x(7) - c(3)*x(1) + c(4)*x(10); 
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(5) + a(2,2)*th(2)*x(6) + a(2,3)*th(3)*x(7) - c(3)*x(2) + c(5)*x(11);
     c(1)*(xi(3) - x(3))*u - a(3,1)*th(1)*x(5) + a(3,2)*th(2)*x(6) - a(3,3)*th(3)*x(7); 
     c(1)*(xi(4) - x(4))*u - a(4,1)*th(1)*x(5)- a(4,2)*th(2)*x(6) - a(4,3)*th(3)*x(7);
     c(1)*(xi(5) - x(5))*u - th(1)*x(5) + a(5,4)*th(4)*x(8); 
     c(1)*(xi(6) - x(6))*u - th(2)*x(6) + a(6,4)*th(4)*x(8);
     c(1)*(xi(7) - x(7))*u - th(3)*x(7) + a(7,4)*th(4)*x(8);
     c(1)*(xi(8) - x(8))*u + a(8,1)*th(1)*x(5) + a(8,2)*th(2)*x(6) + a(8,3)*th(3)*x(7) - th(4)*x(8);
     c(1)*(xi(9) - x(9))*u; 
     c(18)*x(10)^3 + c(19)*x(10)^2*x(11) + c(20)*x(10)*x(11)^2 + c(21)*x(10)^2 + c(22)*x(10)*x(11) + c(23)*x(10) + c(6)*x(1); 
     c(20)*x(11)^3 + c(19)*x(10)*x(11)^2 + c(18)*x(10)^2*x(11) + c(22)*x(11)^2 + c(21)*x(10)*x(11) + c(24)*x(11) + c(6)*x(2)];

dxdt = f; 

end