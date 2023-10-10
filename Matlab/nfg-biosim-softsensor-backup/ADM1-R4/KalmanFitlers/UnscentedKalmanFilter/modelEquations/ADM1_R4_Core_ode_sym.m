%% Version
% (R2022b) Update 5
% Erstelldatum: 06.10.2023
% Autor: Simon Hellmann

function f = ADM1_R4_Core_ode_sym(x, u, xi, th, c, a)
% delivers a symbolic expression (f) of the right-hand side of the ODE 
% system of the core ADM1-R4 (no water, no ash, no nitrogen, no gas phase)
% with constant values for u and xin 

% dynamic equations
f = [c(1)*(xi(1) - x(1))*u + a(1,1)*th(1)*x(3) + a(1,2)*th(3)*x(4) + a(1,3)*th(4)*x(5); % S_ch4 without Henry
     c(1)*(xi(2) - x(2))*u + a(2,1)*th(1)*x(3) + a(2,2)*th(3)*x(4) + a(2,3)*th(4)*x(5); % S_IC without Henry
     c(1)*(xi(3) - x(3))*u - th(1)*x(3) + a(3,4)*th(5)*x(6);    % X_ch
     c(1)*(xi(4) - x(4))*u - th(3)*x(4) + a(4,4)*th(5)*x(6);    % X_pr
     c(1)*(xi(5) - x(5))*u - th(4)*x(5) + a(5,4)*th(5)*x(6);    % X_li
     c(1)*(xi(6) - x(6))*u + a(6,1)*th(1)*x(3) + a(6,2)*th(3)*x(4) + a(6,3)*th(4)*x(5) - th(5)*x(6)];      % X_bac;
end 
