%% Version
% (R2022b) Update 5
% Erstelldatum: 30.08.2023
% last modified: 23.11.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormSlice(xOldNorm,POldNorm, ...
    feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm, ...
    TxNum,TyNum,TuNum,tSpan,nStates)

% compute time & measurement update with Joseph-version of EKF for
% multirate measurements WITHOUT delay, but only some offline measurements 
% available at major instances (Intensiv-Beprobung), normalized coordinates

% xPlusNorm - new norm. state estimate
% PPlusNorm - new norm. state error covariance matrix
% xOldNorm - old norm. state estimate
% POldNorm - old norm. state error covariance matrix
% feedInfoNorm - combined norm. feeding information [tEvents; norm. feedVolFlow; norm. inlet concentrations]
% yMeas - latest measurement vector (minor or major instance) in abs. coordinates
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% QNorm - power spectral density matrix of norm. process noise
% RNorm - covariance matrix of norm. measurement noise
% fNorm - function handle of norm. ODEs of system equations
% gNorm - function handle of norm. output equations 
% dfdxNorm - function handle of norm. partial derivatives of df/dx 
% dhdxNorm - function handle of norm. partial derivatives of df/dx 
% TxNum, TyNum, TuNum - normalization vectors of states, outputs and inputs
% tSpan - time interval between old and new measurement (& state estimate)
% nStates - # states (without sample state augmentation)

global counter

% make print-statement if xOld contains any negative concentrations: 
if nStates > 13 && any(xOldNorm<0)  % beim ADM1-R3 kann S_ion negativ werden
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
    odeFunNorm = @(t,xPNorm) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm, ...
                                params,QNorm,fNorm,dfdxNorm,TxNum,TuNum);
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
        odeFunNorm = @(t,xPNorm) dfP_dtNorm(xPNorm,feedVolFlowNorm,xInCurrNorm, ...
                        params,QNorm,fNorm,dfdxNorm,TxNum,TuNum); 
        [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xP0Norm);
        % update initial value for next interval:
        xP0Norm = xPSolNorm(end,:)';
    end
    xPMinusNorm = xP0Norm;
end

% separate states x and covariance matrix P:
xMinusNorm = xPMinusNorm(1:nStates);
PMinusNorm = reshape(xPMinusNorm(nStates+1:end),[nStates,nStates]);

%% Measurement Update
hNorm = gNorm(xMinusNorm,params.c,TxNum,TyNum);     % full predicted normalized model output
HNorm = dhdxNorm(xMinusNorm,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
yMeasNorm = yMeas'./TyNum; % full normalized measurement vector

% figure out which measurements are available and reduce hNorm, yMeasNorm,
% HNorm and RNorm accordingly: 
idxMeas = ~isnan(yMeas);    % 1 where measurement, 0 where NaN
hNorm = hNorm(idxMeas);         % still col vector
yMeasNorm = yMeasNorm(idxMeas); % still row vector
HNorm = HNorm(idxMeas,:);   
RNorm = RNorm(idxMeas,idxMeas); 

SNorm = HNorm*PMinusNorm*HNorm' + RNorm;    % auxiliary matrix
SNormInv = SNorm\eye(size(RNorm));  % efficient least squares (numerically more efficient alternative for inverse)
KNorm = PMinusNorm*HNorm'*SNormInv; % Kalman Gain matrix

KvNorm = KNorm*(yMeasNorm - hNorm);    % effective correction of Kalman Gain on state estimate (n,1); 
% KvNorm = zeros(nStates,1);  % XY Rania

xPlusNorm = xMinusNorm + KvNorm;    % updated state estimation
% Joseph-Algorithm for measurement update of P: 
leftMat = eye(nStates) - KNorm*HNorm; 
rightMat = leftMat'; 
PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated normalized covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 

end