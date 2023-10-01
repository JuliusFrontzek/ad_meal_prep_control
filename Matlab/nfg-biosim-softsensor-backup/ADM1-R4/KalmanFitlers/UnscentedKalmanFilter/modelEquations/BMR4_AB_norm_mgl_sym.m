%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.2023
% Autor: Simon Hellmann

function g = BMR4_AB_norm_mgl_sym(xNorm,c,Tx,Ty)
% delivers a symbolic expression (h) of normalized measurement equation of the
% BMR4+AB

% measurement equations
g = [(c(6)*Tx(10)^2*xNorm(10)^2 + c(7)*Tx(10)*Tx(11)*xNorm(10)*xNorm(11) + c(8)*Tx(11)^2*xNorm(11)^2 + c(9)*Tx(10)*xNorm(10) + c(10)*Tx(11)*xNorm(11) + c(11))/Ty(1); % volFlow [L/d]
     c(12)*Tx(10)*xNorm(10)/Ty(2);
     c(13)*Tx(11)*xNorm(11)/Ty(3);
     Tx(3)*xNorm(3)/Ty(4);
     (1 - Tx(4)*xNorm(4)/c(14))/Ty(5);
     (1 - Tx(9)*xNorm(9)/(c(14)-Tx(4)*xNorm(4)))/Ty(6)];

end 