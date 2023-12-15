%% Version
% (R2022b) Update 6
% Erstelldatum: 23.11.2023
% last modified: 14.12.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormMultiRateMultiDelay(xOldNorm,POldNorm, ...
    feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm,...
    TxNum,TyNum,TuNum,tSpan,nStates,nAug,actSamplesMat)

% compute time & measurement update with normalized multirate measurements
% with multiple augmentation before return of measurement for samples
% in normalized (norm.) coordinates

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
% nAug - # of augmentations
% actSamplesMat - matrix containing active (=not yet redeemed) sample times
% and return times [tSample, tArrival]

global counter

% make print-statement if xOld contains any negative concentrations: 
if any(xOldNorm<0)
    % print statement: 
    disp(['neg. Konzentration bei ',num2str(tSpan(1)),' Tagen'])
    xOldNorm(xOldNorm<0); 
    % apply clipping: 
    xOldNorm(xOldNorm < 0) = 0; 
    counter = counter + 1
end 

%% Time Update
xPOldNorm = [xOldNorm;reshape(POldNorm,[],1)];  % combined vector of xHat and P

tEvents = feedInfoNorm(:,1);    % feeding time points (an/aus)
% find relevant feeding events (during current measurement interval):
filterRelEvents = tEvents >= tSpan(1) & tEvents < tSpan(2); 
tRelEvents = tEvents(filterRelEvents); 

% integrate across intervals with constant feeding:
% Case a: constant feeding during measurement interval:
if isempty(tRelEvents)
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFunNorm = @(t,xPNorm) dfP_dtAugNormMulti(xPNorm,feedVolFlowNorm,xInCurrNorm, ...
                                params,QNorm,fNorm,dfdxNorm,TxNum,TuNum,...
                                nStates,nAug); 
    [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xPOldNorm);
    xPMinusNorm = xPSolNorm(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    % careful with unique: round to 3 decimal points after the 
    % comma to avoid floating point errors:
    tOverall = unique(round([tSpan;tRelEvents],3));
    nIntervals = length(tOverall) - 1; 
    xP0Norm = xPOldNorm;    % initial value for first interval
    % integrate each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(jj,2); 
        xInCurrNorm = feedInfoNorm(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFunNorm = @(t,xPNorm) dfP_dtAugNormMulti(xPNorm,feedVolFlowNorm,xInCurrNorm, ...
                        params,QNorm,fNorm,dfdxNorm,TxNum,TuNum,...
                        nStates,nAug); 
        [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xP0Norm);
        % update initial value for next interval:
        xP0Norm = xPSolNorm(end,:)';
    end
    xPMinusNorm = xP0Norm;
end

% separate states x and covariance matrix P from LONG vector xPMinus:
if nAug > 0 % augmented case: 
    xAugMinusNorm = xPMinusNorm(1:nStates*(1+nAug));
    PAugMinusNorm = reshape(xPMinusNorm(nStates*(1+nAug)+1:end),[nStates*(1+nAug),nStates*(1+nAug)]);
    xMinusNorm = xAugMinusNorm(1:nStates*(1+nAug)); 
    PMinusNorm = PAugMinusNorm; 
else % non-augmented case:
    xMinusNorm = xPMinusNorm(1:nStates);
    PMinusNorm = reshape(xPMinusNorm(nStates+1:end),[nStates,nStates]);
end

%% Measurement Update
xMinusNormCore = xMinusNorm(1:nStates); % definitely without augmentation
hNormFull = gNorm(xMinusNormCore,params.c,TxNum,TyNum);     % predicted normalized model output
HNormFull = dhdxNorm(xMinusNormCore,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
yMeasNormFull = yMeas'./TyNum; % full normalized measurement vector

% figure out which measurements are available and reduce hNorm, yMeasNorm,
% HNorm and RNorm accordingly: 
idxMeas = ~isnan(yMeas);            % 1 where there is measurement, 0 where NaN
hNorm = hNormFull(idxMeas);         % still col vector
yMeasNorm = yMeasNormFull(idxMeas); % still row vector
HNorm = HNormFull(idxMeas,:);   
RNorm = RNorm(idxMeas,idxMeas); 

% do augmentation if applicable: 
if nAug > 0    
    HNormAug = [HNorm,zeros(nnz(hNorm), nStates*nAug)];
    HNorm = HNormAug; 
end

SNorm = HNorm*PMinusNorm*HNorm' + RNorm;% auxiliary matrix [q,q]
SNormInv = SNorm\eye(size(RNorm));      % efficient least squares (numerically more efficient alternative for inverse)
KNorm = PMinusNorm*HNorm'*SNormInv;     % Kalman Gain matrix [n,q]
KvNorm = KNorm*(yMeasNorm - hNorm);     % effective correction of Kalman Gain on state estimate (n,1); 
% KvNorm = zeros(nStates*(1 + nAug),1);  % XY Rania

% build the indicator matrix for partial state updates appropriately:
if nAug > 0
    % default is to update only original state entries, leave copies 
    % unchanged (no measurement arrival, just samples awaiting return):
    M = blkdiag(eye(nStates),zeros(nStates*nAug)); 
    if any(~isnan(actSamplesMat(:,2))) % if there is measurement coming in
        % find earliest current return time (that one is used in measurement update): 
        [~,idxFirstCurrReturn] = min(actSamplesMat(:,2)); 
        % and switch corresponding entries in M to 1's: 
        indexForM = (idxFirstCurrReturn)*nStates + 1 : (1 + idxFirstCurrReturn)*nStates; 
        M(indexForM,indexForM) = eye(nStates); 
    end
else % non-augmented case: 
    M = eye(nStates); % update all (original) states as usual
end

xPlusNorm = xMinusNorm + M*KvNorm;    % updated normalized state estimation
% Joseph-Algorithm for measurement update of P: 
leftMat = eye(nStates*(1+nAug)) - M*KNorm*HNorm; 
rightMat = leftMat'; 
PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated normalized covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 

end