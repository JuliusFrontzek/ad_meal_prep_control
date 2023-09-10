%% Version
% (R2022b) Update 5
% Erstelldatum: 13.07.2023
% Autor: Simon Hellmann

function  [yMeas,noiseCovMat] = addNoiseToMeasurements(flagModel,yClean)
% adds random noise to clean measurements yClean (normally distributed), 
% and returns noise covariance matrix

% yMeas -       noisy synthetic measurements
% noiseCovMat - noise covariance matrix (diagonal with squared std.
% deviations of sigma on diagonal)
% flagModel -   3: ADM1-R3; 4: ADM1-R4
% yClean -      noise-free outputs

N = size(yClean,1); 
sigmaV = 0.2*24;    % Trommelgaszähler Labor [l/h] --> [l/d]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.12;    % NH4-N [g/L]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]

switch flagModel
    case 3
        q = 8; 
        % define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
        % sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [l/h]
        
        sigmaPh = 0.1;      % pH [-]
        sigmaAc = 0.04;     % FOS [g/L]. Aber Achtung: brauchbarere Messgröße für 
                % S_ac eher die Einzelsäure aus GC, diese hat aber sigma=0,01 g/L 
        
        % combine all in sigma matrix:
        sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh, sigmaSIN, sigmaTS, sigmaVS, sigmaAc]; 
    case 4
        q = 6; 
        % ç in sigma matrix
        sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 

end % switch

% compute covariance matrix of measurement noise: 
noiseCovMat = diag(sigmas.^2);

% create normally distributed measurement noise matrix:
yMean = zeros(N,q); % zero mean for all online measurements 
rng('default');     % fix seed for random number generation (for replicable results)
sigmaMat = repmat(sigmas,N,1);
normalMeasNoise = normrnd(yMean,sigmaMat);
yMeas = yClean + normalMeasNoise; 

end % fun