%% Version
% (R2022b) Update 5
% Erstelldatum: 03.07.2023
% Autor: Simon Hellmann

function g = ADM1_R3_norm_mgl_sym(xNorm,c,Tx,Ty)
% delivers a symbolic expression (h) of normalized measurement equation of 
% the ADM1-R3

% ion balance: 
PhiNorm = Tx(12)*xNorm(12) + (Tx(4)*xNorm(4) - Tx(15)*xNorm(15))/17 - Tx(14)*xNorm(14)/44 - Tx(13)*xNorm(13)/60;
% equivalent proton concentration: 
SHPlusNorm = -PhiNorm/2 + 0.5*sqrt(PhiNorm^2 + c(4)); 

% measurement equations
g = [c(13)*Tx(16)^2/Ty(1)*xNorm(16)^2 + c(14)*Tx(16)*Tx(17)/Ty(1)*xNorm(16)*xNorm(17) + c(15)*Tx(17)^2/Ty(1)*xNorm(17)^2 + c(16)*Tx(16)/Ty(1)*xNorm(16) + c(17)*Tx(17)/Ty(1)*xNorm(17) + c(18)/Ty(1); % volFlow
     c(19)*Tx(16)/Ty(2)*xNorm(16);                  % pch4
     c(20)*Tx(17)/Ty(3)*xNorm(17);                  % pco2
     -1/Ty(4)*log10(SHPlusNorm);                    % pH
     Tx(4)/Ty(5)*xNorm(4);                          % SIN
     1/Ty(6)*(1 - Tx(5)*xNorm(5)/c(21));            % TS    
     1/Ty(7)*(1 - Tx(11)*xNorm(11)/(c(21)-Tx(5)*xNorm(5))); % VS
     Tx(1)/Ty(8)*xNorm(1)];                         % Sac

end 

