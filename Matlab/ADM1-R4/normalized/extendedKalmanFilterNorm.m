%% Version
% (R2022b) Update 
% Erstelldatum: 5.06.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm,Kv] = extendedKalmanFilterNorm(xOldNorm,POldNorm, tSpan,feedInfoNorm, yMeas,params,QNorm,RNorm, fNorm,  gNorm,  dfdxNorm,dhdxNorm,TxNum,TyNum,TuNum)
                                                           % (xOld,     POld,       tSpan,feedInfo,     yMeas,params,Q,     R,      f,      g,      dfdx,    dhdx)
% compute time and measurement update acc. to Joseph-version of the EKF in
% normalized coordinates

% XY: evtl. noch zu viele Variablen übergeben
% XY: Variablen-Beschreibung anpassen

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise
% f - function handle of ODEs of system equations
% g - function handle of output equations 
% dfdx - function handle of partial derivatives of df/dx 
% dhdx - function handle of partial derivatives of df/dx 

global counter

% make print-statement if xOld contains any negative concentrations: 
if any(xOldNorm<0)
    % print statement: 
    disp(['neg. Konzentration bei ',num2str(tSpan(1)),' Tagen'])
    xOldNorm(xOldNorm<0)
    % apply clipping: 
    xOldNorm(xOldNorm < 0) = 0; 
    counter = counter + 1
end 

%% Time Update
xPOldNorm = [xOldNorm;reshape(POldNorm,[],1)];  % combined vector of xHat and P

tEvents = feedInfoNorm(:,1);    % feeding time points (an/aus)
% find relevant feeding events (during current measurement interval):
filterRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2); 
tRelEvents = tEvents(filterRelEvents); 

% integrate across intervals with constant feeding:
% Case a: constant feeding during measurement interval:
if isempty(tRelEvents)
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFunNorm = @(t,xPNorm) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm,params,QNorm,fNorm,dfdxNorm,TxNum,TuNum);    % XY Abhängigkeiten korrigieren!
    [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xPOldNorm);
    xPMinusNorm = xPSolNorm(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    xP0Norm = xPOldNorm;    % initial value for first interval
    % integrate each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(jj,2); 
        xInCurrNorm = feedInfoNorm(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFunNorm = @(t,xP) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm,params,QNorm,fNorm,dfdxNorm,TxNum,TuNum);    % XY Abhängigkeiten korrigieren!
        [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xP0Norm);
        % update initial value for next interval:
        xP0Norm = xPSolNorm(end,:)';
    end
    xPMinusNorm = xP0Norm;
end

% separate states and covariance matrix:
nStates = length(xOldNorm); 
xMinusNorm = xPMinusNorm(1:nStates);
PMinusNorm = reshape(xPMinusNorm(nStates+1:end),[nStates,nStates]);

%% Measurement Update
HNorm = dhdxNorm(xMinusNorm,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
SNorm = HNorm*PMinusNorm*HNorm' + RNorm;    % auxiliary matrix
SNormInv = SNorm\eye(size(R));  % numerically more efficient way to compute inverse
KNorm = PMinusNorm*HNorm'*SNormInv; % Kalman Gain matrix
hNorm = gNorm(xMinusNorm,params.c,TxNum,TyNum); % predicted model output
yMeasNorm = yMeas./TyNum;           % normalized measurements
Kv = KNorm*(yMeasNorm' - hNorm);    % effective correction of Kalman Gain on state estimate (n,1); 
% Kv = zeros(nStates,1);  % XY Rania

xPlusNorm = xMinusNorm + Kv;    % updated state estimation

% apply Joseph-Algorithm: 
leftMat = eye(nStates) - KNorm*HNorm; 
rightMat = leftMat'; 
PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated covariance of estimation error

% old: without Joseph Algorithm:
% PPlus = PMinus - K*H*PMinus; % updated covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 
end
