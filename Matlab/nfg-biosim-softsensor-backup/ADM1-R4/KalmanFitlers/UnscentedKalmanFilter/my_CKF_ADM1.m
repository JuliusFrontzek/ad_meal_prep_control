%% Version
% (R2022b) Update 5
% Erstelldatum: 04.10.2023
% Autor: Julius Frontzek, Simon Hellmann

function [xPlus,PPlus] = my_CKF_ADM1(ckf,feedInfo,yMeas,tSpan,params,f,g) 
% XY: argumente anpassen an das, was wirklich gebraucht wird

% compute time and measurement update acc. to Cubature Kalman Filter

% XY: argumente und Beschreibungen anpassen: 

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p - parameter vector for model the UKF is working with
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise

[xMinus,PMinus] = predict(ckf,tSpan,feedInfo,params,f);

% yMeas = yMeas([1,2,4:end]);   % exclude pCO2
[xPlus,PPlus] = correct(ckf,yMeas,g,params.c);


end