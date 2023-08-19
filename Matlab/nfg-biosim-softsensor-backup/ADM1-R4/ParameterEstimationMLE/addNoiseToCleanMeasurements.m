%% Version
% (R2022b) Update 2
% Erstelldatum: 02.03.23
% Autor: Simon Hellmann

function [yMeas] = addNoiseToCleanMeasurements(yClean, C)

% number of measurements: 
[N,q] = size(yClean);

% create mean and std.-deviation for normally distributed random noise:
yMean = zeros(N,q);      % zero mean online measurements
% meanOff = zeros(NOff,3);    % zero mean offline measurements 
sigma = sqrt(diag(C)');     % extract only std deviations 
sigmaMat = repmat(sigma,N,1); 
sigmaOff = repmat(sigma(4:6),NOff,1); 
noiseOn = normrnd(yMean,sigmaMat);
noiseOff = normrnd(meanOff,sigmaOff);

% corrupt clean measurements with random noise (normal distribution): 
yMeasOn = yCleanOn + noiseOn; 
yMeasOff = yCleanOff + noiseOff; 

end