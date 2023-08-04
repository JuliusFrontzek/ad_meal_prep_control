%% Version
% (R2022b) Update 5
% Erstelldatum: 13.07.2023
% Autor: Simon Hellmann

function  [yMeas,noiseCovMat] = addNoiseToMeasurements(flagModel,yClean)
% XY: Zweck der Funktion und Argumente beschreiben

N = size(yClean,1); 

switch flagModel
    case 3
        q = 8; 
        % define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
        % sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [l/h]
        sigmaV = 0.2*24;    % Trommelgaszähler Labor [l/h] --> [l/d]
        sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
        sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
        sigmaPh = 0.1;      % pH [-]
        sigmaSIN = 0.12;    % NH4-N [g/L]
        sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
        sigmaVS = 0.31/100; % [%] -> [-]
        sigmaAc = 0.04;     % FOS [g/L]. Aber Achtung: brauchbarere Messgröße für 
                % S_ac eher die Einzelsäure aus GC, diese hat aber sigma=0,01 g/L 
        
        % combine all in sigma matrix and covariance matrix:
        sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh, sigmaSIN, sigmaTS, sigmaVS, sigmaAc]; 
        sigmaMat = repmat(sigmas,N,1);
        noiseCovMat = diag(sigmas.^2);  % measurement noise covariance matrix
        
        % create normally distributed measurement noise matrix:
        yMean = zeros(N,q); % zero mean for all online measurements 
        rng('default');     % fix seed for random number generation (for replicable results)
        normalMeasNoise = normrnd(yMean,sigmaMat);
        yMeas = yClean + normalMeasNoise; 

end % switch
end % fun