%% Version
% (R2022b) Update 6
% Erstelldatum: 23.11.2023
% last modified: 23.11.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormMultiRateMultiDelay(xOldNorm,POldNorm, ...
    feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm, ...
    TxNum,TyNum,TuNum,tSpan,nStates,qOn,qOff,nAug,actSamplesMat,...
    flagDelayPeriod,flagArrival)

% compute time & measurement update with normalized multirate measurements
% with multiple augmentation before return of measurement for samples
% in normalized (norm.) coordinates

% XY: überlege, wie du die Werte aus actSamplesMat gut umsetzt. Bedenke,
% dass die Matrix auch empty sein kann!

% XY: nomenclatur anpassen!

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
% qOn, qOff - # online and offline signals
% flagAugmented -   0: non-augmented,   1: augmented
% flagDelayPeriod - 0: NOT waiting for lab measurement to return. 1: waiting
% flagArrival -     0: minor instance,  1: major instance

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
    odeFunNorm = @(t,xPNorm) dfP_dtAugNorm(xPNorm,feedVolFlowNorm,xInCurrNorm, ...
                                params,QNorm,fNorm,dfdxNorm,TxNum,TuNum,...
                                nStates,nAug);
    [~,xPSolNorm] = ode15s(odeFunNorm,tEval,xPOldNorm);
    xPMinusNorm = xPSolNorm(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    tOverall = unique(sort([tSpan;tRelEvents]));
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
xMinusCoreNorm = xMinusNorm(1:nStates); % definitely without augmentation
hNormPre = gNorm(xMinusCoreNorm,params.c,TxNum,TyNum);     % predicted normalized model output
HNormPre = dhdxNorm(xMinusCoreNorm,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
yMeasNorm = yMeas'./TyNum; % full normalized measurement vector

% figure out which measurements are available and reduce hNorm, yMeasNorm,
% HNorm and RNorm accordingly: 
idxMeas = ~isnan(yMeas);        % 1 where there is measurement, 0 where NaN
hNorm = hNormPre(idxMeas);         % still col vector
yMeasNorm = yMeasNorm(idxMeas); % still row vector
HNorm = HNormPre(idxMeas,:);   
RNorm = RNorm(idxMeas,idxMeas); 


% XY: die folgenden 4 Zeilen Code kommen wahrscheinlich erst später!
SNorm = HNorm*PMinusNorm*HNorm' + RNorm;    % auxiliary matrix
SNormInv = SNorm\eye(size(RNorm));  % efficient least squares (numerically more efficient alternative for inverse)
KNorm = PMinusNorm*HNorm'*SNormInv; % Kalman Gain matrix
KvNorm = KNorm*(yMeasNorm - hNorm);    % effective correction of Kalman Gain on state estimate (n,1); 
% XY Ende der vier Zeilen Code


% XY: rekonstruiere die korrekte Indicator Matrix aus der active samples
% matrix!
if nAug > 0 % augmented case:    
%     HAugNorm = nan(qOn+qOff,nStates*(1+nAug));     % allocate memory
%     HAugNorm(:,1:nStates) = dhdxNorm(xMinusCoreNorm,params.c,TxNum,TyNum);  % the remaining nStates cols stay 0
%     % fill augmented columns with zeros since they correspond to sample
%     % state which is fixed, thus currently not time dependent anymore:
%     HAugNorm(:,nStates+1:end) = zeros(qOn+qOff, nStates); 
    
    % XY: checke, ob die Dimensionen korrekt sind, auch wenn vorher schon 
    % nach den laut idxMeas relevanten Einträgen gesliced wurde!
    HAugNorm = [HNorm,zeros(qOn+qOff, nStates*nAug)];
    HNorm = HAugNorm; 
% else % non-augmented case:
%     HNorm = dhdxNorm(xMinusCoreNorm,params.c,TxNum,TyNum);  % partial derivative of output, evaluated at xMinus
end

% Slicing. During minor instance...
if flagArrival == 0
    % ... reduce H and R to # of measurement signals during minor instance:
    HNorm = HNorm(1:qOn,:); 
    RNorm = RNorm(1:qOn,1:qOn); 
    % ... extract reduced # of measurement signals during minor instance:
    hNorm = hNorm(1:qOn); 
    yMeasNorm = yMeas(1:qOn)'./TyNum(1:qOn); % reduced normalized measurement vector
else 
    yMeasNorm = yMeas'./TyNum; % full normalized measurement vector
end
% during major instance, you must use all of the measurements available
% (merged online and offline measurements together)

SNorm = HNorm*PMinusNorm*HNorm' + RNorm;    % auxiliary matrix
SNormInv = SNorm\eye(size(RNorm));  % efficient least squares (numerically more efficient alternative for inverse)
KNorm = PMinusNorm*HNorm'*SNormInv; % Kalman Gain matrix

KvNorm = KNorm*(yMeasNorm - hNorm);    % effective correction of Kalman Gain on state estimate (n,1); 
% KvNorm = zeros(nStates,1);  % XY Rania

% XY: passe die indicator matrix an!
if flagDelayPeriod == 1
    M = blkdiag(eye(nStates),zeros(nStates)); % indicator matrix. XY anpassen!
    xPlusNorm = xMinusNorm + M*KvNorm;    % updated normalized state estimation
    % Joseph-Algorithm for measurement update of P: 
    leftMat = eye(nStates*(1+nAug)) - M*KNorm*HNorm; 
    rightMat = leftMat'; 
    PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated normalized covariance of estimation error
else % when not waiting for lab measurement to return, update all of x and P, 
    % regardless whether augmented or not:
    xPlusNorm = xMinusNorm + KvNorm;    % updated state estimation
    % Joseph-Algorithm for measurement update of P: 
    leftMat = eye(size(xMinusNorm)) - KNorm*HNorm; 
    rightMat = leftMat'; 
    PPlusNorm = leftMat*PMinusNorm*rightMat + KNorm*RNorm*KNorm'; % updated normalized covariance of estimation error
end

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 

end