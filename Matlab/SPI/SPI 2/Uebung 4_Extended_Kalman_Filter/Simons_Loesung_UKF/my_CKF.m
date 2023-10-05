%% Version
% (R2022b) Update 5
% Erstelldatum: 04.10.2023
% Autor: Julius Frontzek, Simon Hellmann

function [xPlus,PPlus] = my_CKF(ckf,u,yMeas,tSpan,p_CKF,testParam)

% compute time and measurement update acc. to Cubature Kalman Filter and
% saves also saves the results in ckf-object

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% ckf - object containing all ckf parameters and estimations
% u - current feed volume flow
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p_CKF - parameter vector for model the CKF is working with
% testParam - additional argument of the measurement function

[xMinus,PMinus] = predict(ckf, u, tSpan, p_CKF);
[xPlus,PPlus] = correct(ckf,yMeas,testParam);

end