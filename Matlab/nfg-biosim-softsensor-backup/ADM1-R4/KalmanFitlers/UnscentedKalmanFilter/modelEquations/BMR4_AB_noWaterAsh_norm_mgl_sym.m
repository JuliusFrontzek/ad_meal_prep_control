%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.2023
% Autor: Simon Hellmann

function g = BMR4_AB_noWaterAsh_norm_mgl_sym(xNorm,c,Tx,Ty)
% delivers a symbolic expression (h) of normalized measurement equation of 
% the BMR4+AB_noWaterAsh

% measurement equations
g = [(c(6)*Tx(8)^2*xNorm(8)^2 + c(7)*Tx(8)*Tx(9)*xNorm(8)*xNorm(9) + c(8)*Tx(9)^2*xNorm(9)^2 + c(9)*Tx(8)*xNorm(8) + c(10)*Tx(9)*xNorm(9) + c(11))/Ty(1); % volFlow [L/d]
     c(12)*Tx(8)*xNorm(8)/Ty(2);
     c(13)*Tx(9)*xNorm(9)/Ty(3);
     Tx(3)*xNorm(3)/Ty(4)];
end 