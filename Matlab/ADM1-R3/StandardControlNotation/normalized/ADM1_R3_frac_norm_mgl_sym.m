%% Version
% (R2022b) Update 5
% Erstelldatum: 11.5.2023
% Autor: Simon Hellmann

function g = ADM1_R3_frac_norm_mgl_sym(xNorm,c,Tx,Ty)
% delivers a symbolic expression (h) of normalized measurement equation of 
% the ADM1-R3-frac

% ion balance: 
PhiNorm = Tx(13)*xNorm(13) + (Tx(4)*xNorm(4) - Tx(16)*xNorm(16))/17 - Tx(15)*xNorm(15)/44 - Tx(14)*xNorm(14)/60;
% equivalent proton concentration: 
SHPlusNorm = -PhiNorm/2 + 0.5*sqrt(PhiNorm^2 + c(4)); 

% measurement equations
g = [c(13)*Tx(17)^2/Ty(1)*xNorm(17)^2 + c(14)*Tx(17)*Tx(18)/Ty(1)*xNorm(17)*xNorm(18) + c(15)*Tx(18)^2/Ty(1)*xNorm(18)^2 + c(16)*Tx(17)/Ty(1)*xNorm(17) + c(17)*Tx(18)/Ty(1)*xNorm(18) + c(18)/Ty(1); % volFlow
     c(19)*Tx(17)/Ty(2)*xNorm(17);                  % pch4
     c(20)*Tx(18)/Ty(3)*xNorm(18);                  % pco2
     -1/Ty(4)*log10(SHPlusNorm);                    % pH
%    Tx(4)/Ty(5)*xNorm(4);                          % SIN
     Tx(4)/Ty(5)*xNorm(4) - Tx(16)/Ty(5)*xNorm(16); % Snh4
     1/Ty(6)*(1 - Tx(5)*xNorm(5)/c(21));            % TS    
     1/Ty(7)*(1 - Tx(12)*xNorm(12)/(c(21)-Tx(5)*xNorm(5))); % VS
     Tx(1)/Ty(8)*xNorm(1)];                         % Sac

end 