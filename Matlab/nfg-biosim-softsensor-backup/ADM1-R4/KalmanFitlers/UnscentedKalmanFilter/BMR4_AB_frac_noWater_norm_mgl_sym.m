%% Version
% (R2022b) Update 5
% Erstelldatum: 19.09.2023
% Autor: Simon Hellmann

function g = BMR4_AB_frac_noWater_norm_mgl_sym(xNorm,c,Tx,Ty)
% delivers a symbolic expression (h) of normalized measurement equation of 
% the BMR4+AB_frac_noWater

% measurement equations
g = [(c(6)*Tx(9)^2*xNorm(9)^2 + c(7)*Tx(9)*Tx(10)*xNorm(9)*xNorm(10) + c(8)*Tx(10)^2*xNorm(10)^2 + c(9)*Tx(9)*xNorm(9) + c(10)*Tx(10)*xNorm(10) + c(11))/Ty(1); % volFlow [L/d]
     c(12)*Tx(9)*xNorm(9)/Ty(2);
     c(13)*Tx(10)*xNorm(10)/Ty(3);
     Tx(3)*xNorm(3)/Ty(4)];
end 