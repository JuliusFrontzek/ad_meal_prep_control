%% Version
% (R2022b) Update 5
% Erstelldatum: 30.09.2023
% last modified: 30.09.2023
% Autor: Simon Hellmann

function NSE = computeNSE_groundTruth(yMeas,yEst,yTrue)
% computes the Nash-Sutcliffe efficiency based on sum of squares, but for
% instead of mean of measurements, use the ground truth (which is only
% known in simulation studies)

% compute sum of squares of numerator and denominator:
sumOfsquaresNum = (yMeas - yEst)'*(yMeas - yEst); 
sumOfsquaresDenom = (yMeas - yTrue)'*(yMeas - yTrue); 

NSE = 1 - (sumOfsquaresNum/sumOfsquaresDenom); 

end