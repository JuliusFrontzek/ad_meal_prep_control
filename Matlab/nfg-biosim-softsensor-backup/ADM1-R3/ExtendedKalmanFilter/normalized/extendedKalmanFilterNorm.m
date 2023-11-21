%% Version
% (R2022b) Update 
% Erstelldatum: 5.06.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = extendedKalmanFilterNorm(xOldNorm,POldNorm,tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm,TxNum,TyNum,TuNum)
% compute time and measurement update acc. to Joseph-version of the EKF in
% normalized coordinates

% xPlusNorm - new normalized state estimate
% PPlusNorm - new normalized state error covariance matrix
% xOldNorm - old normalized state estimate
% POldNorm - old normalized state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfoNorm - combination of normalized feeding information [tEvents; normalized feedVolFlow; normalized inlet concentrations]
% yMeas - latest measurement vector (non normalized)
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% QNorm - normalized power spectral density matrix of process noise
% RNorm - normalized covariance matrix of measurement noise
% fNorm - function handle of normalized ODEs of system equations
% gNorm - function handle of normalized output equations 
% dfdxNorm - function handle of normalized partial derivatives of df/dx 
% dhdxNorm - function handle of normalized partial derivatives of df/dx 
% TxNum - normalization matrix of states
% TyNum - normalization matrix of outputs
% TuNum - normalization matrix of inputs

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
% combine vectors of xHat and P:
xPOldNorm = [xOldNorm;reshape(POldNorm,[],1)];  

tEvents = feedInfoNorm(:,1);    % feeding time points (on/off)
% find relevant feeding events (during current measurement interval):
filterRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2); 
tRelEvents = tEvents(filterRelEvents); 

% integrate across intervals with constant feeding:
% Case a: constant feeding during measurement interval:
if isempty(tRelEvents)
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFunNorm = @(t,xPNorm) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm,params,QNorm,fNorm,dfdxNorm,TxNum,TuNum);
    [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xPOldNorm);
    xPMinusNorm = xPSolNorm(end,:)';

% Case b: feeding changes during measurement interval:
else 
    % create time grid of currently relevant feeding events and measurements:
    tOverall = unique(sort([tSpan; tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    xP0Norm = xPOldNorm;    % initial value for first interval
    % integrate across each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(jj,2); 
        xInCurrNorm = feedInfoNorm(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFunNorm = @(t,xPNorm) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm,params,QNorm,fNorm,dfdxNorm,TxNum,TuNum); 
        [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xP0Norm);
        % update initial value for next interval:
        xP0Norm = xPSolNorm(end,:)';
    end
    xPMinusNorm = xP0Norm;
end

% separate states x and covariance matrix P again:
nStates = length(xOldNorm); 
xMinusNorm = xPMinusNorm(1:nStates);
PMinusNorm = reshape(xPMinusNorm(nStates+1:end),[nStates,nStates]);

%% Measurement Update
HNorm = dhdxNorm(xMinusNorm,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
SNorm = HNorm*PMinusNorm*HNorm' + RNorm;    % auxiliary matrix
SNormInv = SNorm\eye(size(RNorm));  % numerically more efficient way to compute inverse
KNorm = PMinusNorm*HNorm'*SNormInv; % Kalman Gain matrix
hNorm = gNorm(xMinusNorm,params.c,TxNum,TyNum); % normalized predicted model output
yMeasNorm = yMeas'./TyNum;          % normalized measurements
KvNorm = KNorm*(yMeasNorm - hNorm);     % effective correction of Kalman Gain on state estimate (n,1); 
% KvNorm = zeros(nStates,1);              % XY Rania

xPlusNorm = xMinusNorm + KvNorm;    % updated state estimation

% apply Joseph-Algorithm for measurement update of P-matrix: 
leftMat = eye(nStates) - KNorm*HNorm; 
rightMat = leftMat'; 
PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated covariance of estimation error

% old: without Joseph Algorithm:
% PPlus = PMinus - K*H*PMinus; % updated covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 
end