%% Version
% (R2022b) Update 5
% Erstelldatum: 30.09.2023
% last modified: 30.09.2023
% Autor: Simon Hellmann

function NSE = computeNSE_abs_groundTruth(yMeas,yEst,yTrue)
% computes the Nash-Sutcliffe efficiency based on sum abs. values, but for
% instead of mean of measurements, use the ground truth (which is only
% known in simulation studies)

absNum = abs(yMeas - yEst); 
absDenom = abs(yMeas - yTrue); 

% compute sum of abs. values of numerator and denominator:
sumOfAbsNum = sum(absNum); 
sumOfAbsDenom = sum(absDenom); 

NSE = 1 - (sumOfAbsNum/sumOfAbsDenom); 

end