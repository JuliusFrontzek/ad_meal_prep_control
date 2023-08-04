%% Version
% (R2022b) Update 2
% Erstelldatum: 2.3.2023
% Autor: Simon Hellmann

function  [RMSSE] = computeRMSSE(yMeas,yEst)

N = numel(yMeas);   % number of samples

yDiff = yEst - yMeas;   
yAuto = diff(yMeas);    % naive forecast based on last measurement

squaredNum = 1/N*(yDiff'*yDiff); 
squaredDenom = 1/(N-1)*(yAuto'*yAuto); 

RMSSE = sqrt(squaredNum/squaredDenom); 

end