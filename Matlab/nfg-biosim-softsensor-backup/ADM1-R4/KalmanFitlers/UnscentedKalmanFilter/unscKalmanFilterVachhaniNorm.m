%% Version
% (R2022b) Update 5
% Erstelldatum: 31.08.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = unscKalmanFilterVachhaniNorm(xOldNorm,POldNorm, ...
                        tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm, ...
                        fNorm,gNorm,TxNum,TyNum,TuNum)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Vachhani et al. (2006) but withoug optimzing measurement update (normalized)

global counterSigmaInitNorm
global counterSigmaPropNorm
global counterXNorm

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
% TxNum - normalization matrix of states
% TyNum - normalization matrix of outputs
% TuNum - normalization matrix of inputs

% extract constant parameters out of struct: 
th = params.th; 
c = params.c; 
a = params.a;

% clipping if xOldNorm contains any negative concentrations: 
if any(xOldNorm<0)
    % apply clipping: 
    xOldNorm(xOldNorm < 0) = 0; 
    counterXNorm = counterXNorm + 1;
end 

nStates = length(xOldNorm); 
q = numel(yMeas); 
nSigmaPoints = 2*nStates + 1; 

%% Time Update (TU)

% define scaling parameters and weights: 
kappa = 0.05;  % default (Vachhani 2006): 1
gamma = sqrt(nStates + kappa); % scaling parameter

% % weights acc. to Julier: 
% w0 = kappa/(nStates + kappa); 
% wi = 1/(2*(nStates + kappa)); 
% w = [w0, repmat(wi,1,2*nStates)]; % for state aggregation

%% Choose Sigma Points while respecting state constraints
sqrtPOldNorm = schol(POldNorm);  % cholesky factorization

% check if state constraints are violated. If so, modify scaling parameters: 
sijNorm = [sqrtPOldNorm, -sqrtPOldNorm];    % 2n Einträge

xMinNorm = zeros(nStates,1);    % no negative concentrations
xMaxNorm = 1E2*ones(nStates,1); % there is theoretically no upper limit on concentrations, but set a practical one
% apply clipping for water concentration:
xMaxNorm(4) = 990/TxNum(4); % g/l

theta_ik = nan(1,2*nStates);     % allocate memory (früher: 1,nStates)
% go through sij column by column to recompute the scaling paramters of all 
% sigma points: 
for col = 1:2*nStates % nStates
    currCol = sijNorm(:,col);   % current column
    posScalingIdx = currCol > 0; 
    negScalingIdx = currCol < 0; 
    
    posScale = currCol(posScalingIdx); 
    % count positive scaling parameters into theta_1k: 
    posRescaling = (xMaxNorm(posScalingIdx) - xOldNorm(posScalingIdx))./posScale; 
    % hier kann nichts negatives rauskommen, also abs nicht mehr nötig
    theta_1k = min(posRescaling);% alt: min(abs(posRescaling)); 

    negScale = currCol(negScalingIdx); 
    % count all other scaling parameters into theta_2k: 
    negRescaling = (xMinNorm(negScalingIdx) - xOldNorm(negScalingIdx))./negScale; 
    theta_2k = min(negRescaling); % alt: min(abs(negRescaling));
    theta_ik(col) = min([gamma, theta_1k, theta_2k]);
end

sigmaXInitNorm = [xOldNorm, repmat(xOldNorm,1,2*nStates) + theta_ik.*sijNorm]; 

% Apply clipping to near-negative Sigma Points: 
epsilon = 1e-10; 
if any(any(sigmaXInitNorm < -epsilon))
    sigmaXInitNorm(sigmaXInitNorm < 0) = 0; 
    counterSigmaInitNorm = counterSigmaInitNorm + 1;
end

%% Propagate Sigma Points
sigmaXPropNorm = nan(nStates, 2*nStates + 1); % allocate memory

tEvents = feedInfoNorm(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
%     feedVectorNorm = feedInfoNorm; % verwende die aktuell wirksame Fütterung
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFunNorm = @(t,XNorm) fNorm(XNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
    for k = 1:nSigmaPoints
        [~,XTUSolNorm] = ode15s(odeFunNorm,tEval,sigmaXInitNorm(:,k));
        sigmaXPropNorm(:,k) = XTUSolNorm(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtEndOfIntNorm = sigmaXInitNorm;    % Startwerte für erstes Intervall (wird überschrieben)
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(m,2);
        xInCurrNorm = feedInfoNorm(m,3:end)';   % current inlet concentrations
        tEval = [tOverall(m), tOverall(m+1)];
        odeFunNorm = @(t,XNorm) fNorm(XNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
        for kk = 1:nSigmaPoints
            [~,XTUSolNorm] = ode15s(odeFunNorm,tEval,XAtEndOfIntNorm(:,kk));
            XAtEndOfIntNorm(:,kk) = XTUSolNorm(end,:)';
        end
    end
    sigmaXPropNorm = XAtEndOfIntNorm;
end

% if any propagated sigma points violate constraints, apply clipping: 
epsilon = 1e-4; 
if any(any(sigmaXPropNorm < -epsilon))
    sigmaXPropNorm(sigmaXPropNorm < 0) = 0; 
    counterSigmaPropNorm = counterSigmaPropNorm + 1;
end

%% compute weights acc. to Vachhani 2006: 
STheta = sum(theta_ik); 
a = (2*kappa - 1)/(2*(nStates + kappa)*(STheta - nSigmaPoints*sqrt(nStates + kappa)));
b = 1/(2*(nStates + kappa)) - (2*kappa - 1)/(2*sqrt(nStates + kappa)*(STheta - nSigmaPoints*sqrt(nStates + kappa)));
WiTemp = a*theta_ik + b; 
% because numel(Wi) = nSigmaPoints - 1 and sum != 1, we need a zero'th weight: 
Wi0 = 1 - sum(WiTemp); 
Wi = [Wi0,WiTemp]; 

%% Aggregate Sigma Points to Priors for x and P
% (do NOT distinguish between weighting of states and covariances)

% aggregate state prior:
xMinusNorm = sum(Wi.*sigmaXPropNorm,2);  

% aggregate state error cov. matrix P:
diffXPriorFromSigmaNorm = sigmaXPropNorm - xMinusNorm; 
PMinusNorm = Wi.*diffXPriorFromSigmaNorm*diffXPriorFromSigmaNorm' + QNorm; 

%% Measurement Update (MU):
% omit to choose new sigma points for measurement update (acc. to Vachhani
% 2006)!

%% Derive Sigma-Measurements and aggregate them:
YNorm = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    YNorm(:,mm) = gNorm(sigmaXPropNorm(:,mm),c,TxNum,TyNum); 
end
% aggregate outputs of sigma points in overall output:
yAggregatedNorm = sum(Wi.*YNorm,2);

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputsNorm = YNorm - yAggregatedNorm; 
PyyNorm = Wi.*diffYFromSigmaOutputsNorm*diffYFromSigmaOutputsNorm' + RNorm;

% compute cross covariance matrix states/measurements:
PxyNorm = Wi.*diffXPriorFromSigmaNorm*diffYFromSigmaOutputsNorm'; 

%% 4. Kalman Gain and actual MU:
PyyNormInv = PyyNorm\eye(q); 
KNorm = PxyNorm*PyyNormInv; 
yMeasNorm = yMeas'./TyNum;   % normalized measurements
KvNorm = KNorm*(yMeasNorm - yAggregatedNorm); 

xPlusNorm = xMinusNorm + KvNorm;
PPlusTempNorm = PMinusNorm - KNorm*PyyNorm*KNorm'; 

% make sure PPlusNorm is symmetric: 
PPlusNorm = 1/2*(PPlusTempNorm + PPlusTempNorm'); 

end
