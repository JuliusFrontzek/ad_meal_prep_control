%% Version
% (R2022b) Update 5
% Erstelldatum: 04.10.2023
% Autor: Julius Frontzek, Simon Hellmann

function [xPlus,PPlus] = my_UKF(ukf,u,yMeas,tSpan,p_UKF,testParam)

% compute time and measurement update acc. to Unscented Kalman Filter of 
% Matlab and saves also saves the results in ukf-object

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% ukf - object containing all ukf parameters and estimations
% u - current feed volume flow
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p_UKF - parameter vector for model the CKF is working with
% testParam - additional argument of the measurement function

[xMinus,PMinus] = predict(ukf, u, tSpan, p_UKF);
[xPlus,PPlus] = correct(ukf,yMeas,testParam);

end