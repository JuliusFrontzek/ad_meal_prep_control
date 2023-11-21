%% Version
% (R2022b) Update 6
% Erstelldatum: 19.11.2023
% Autor: Simon Hellmann

function  [NSE, eNSE] = compute_NSE(yMeas,yEst)

numerator_NSE = (yMeas - yEst)'*(yMeas - yEst);
yMean = mean(yMeas); 
denominator_NSE = (yMeas - yMean)'*(yMeas - yMean);

numerator_eNSE = sum(abs(yMeas - yEst));
denominator_eNSE = sum(abs(yMeas - yMean));

NSE = 1 - numerator_NSE/denominator_NSE; 
eNSE = 1 - numerator_eNSE/denominator_eNSE; 

end