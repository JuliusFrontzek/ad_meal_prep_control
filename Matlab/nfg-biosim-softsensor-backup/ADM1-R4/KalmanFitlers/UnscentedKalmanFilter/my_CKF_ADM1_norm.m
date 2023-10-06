%% Version
% (R2022b) Update 5
% Erstelldatum: 04.10.2023
% Autor: Julius Frontzek, Simon Hellmann

function [xPlusNorm,PPlusNorm] = my_CKF_ADM1_norm(ckf,feedInfoNorm,yMeas, ...
                                        tSpan,params,fNorm,gNorm,TxNum,TyNum,TuNum) 
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

[xMinusNorm,PMinusNorm] = predict(ckf,tSpan,feedInfoNorm,params,fNorm,TxNum,TuNum);
    % alte stateTansFun: StateTransitionFcn(xMinus, u, tSpan, p)
    % alte argumente:                      (ckf,    u, tSpan, p_CKF);
    % ----
    % neue stateTransFun: stateTransitionFcnADM1R4(xMinus,tSpan,feedInfo,params,f)
    % neue argumente:                             (ckf,   tSpan,feedInfo,params,f);
    
yMeasNorm = yMeas'./TyNum;   % normalized measurements
[xPlusNorm,PPlusNorm] = correct(ckf,yMeasNorm,gNorm,params.c,TxNum,TyNum);


end