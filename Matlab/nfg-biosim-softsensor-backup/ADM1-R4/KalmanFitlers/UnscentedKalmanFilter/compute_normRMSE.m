%% Version
% (R2022b) Update 6
% Erstelldatum: 12.10.2023
% Autor: Simon Hellmann

% compute the normalized root mean squared error

function  nRMSE = compute_normRMSE(xTrue,xEst)

N = numel(xTrue);   % number of samples

xDiff = xEst - xTrue;

squaredNum = xDiff'*xDiff;      % sum of squares
squaredDenom = xTrue'*xTrue;    % sum of squares 

nRMSE = sqrt(1/N*squaredNum/squaredDenom); 

end