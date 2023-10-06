%% Version
% (R2022b) Update 5
% Erstelldatum: 04.10.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = unscKalmanFilterKolasFullyAugmentedNorm(xOldNorm,POldNorm, ...
                                    tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm, ...
                                    fNorm,gNorm,TxNum,TyNum,TuNum)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) with clipping wherever possible (normalized),
% but for augmented process noise (Tab. 6)

global counterSigmaInitNorm
global counterSigmaPropNorm
global counterSigmaXNorm
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

% If xOldNorm contains negative concentrations, apply clipping: 
if any(xOldNorm<0)
    xOldNorm(xOldNorm < 0) = 0; 
    counterXNorm = counterXNorm + 1;
end 

nStates = numel(xOldNorm); 
q = numel(yMeas); 

% augment x and P: 
xOldAugNorm = [xOldNorm;zeros(nStates,1);zeros(q,1)]; 
POldAugNorm = blkdiag(POldNorm,QNorm,RNorm);  % (2*nStates+q, 2*nStates+q)

nStatesAug = numel(xOldAugNorm); 
nSigmaPointsAug = 2*(nStatesAug) + 1;   % # sigma points with augmentation

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
% beta = 0; 
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
% kappa = 3 - nStatesAug;  % acc. to Julier & Uhlmann
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
gamma = sqrt(nStatesAug + lambda); % scaling parameter
% gamma = 0.2;  % XY just to check

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStatesAug + lambda); 
Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStatesAug + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPointsAug-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPointsAug-1)]; % for covariance aggregation

%% Choose Sigma Points
sqrtPOldAugNorm = schol(POldAugNorm);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInitNorm = [xOldAugNorm,  repmat(xOldAugNorm,1,nStatesAug) + gamma*sqrtPOldAugNorm, ...
                                repmat(xOldAugNorm,1,nStatesAug) - gamma*sqrtPOldAugNorm]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInitNorm < 0))
    sigmaXInitNorm(sigmaXInitNorm < 0) = 0; 
    counterSigmaInitNorm = counterSigmaInitNorm + 1;
end

%% Propagate Sigma Points
sigmaXPropNormNom = nan(nStates, nSigmaPointsAug); % allocate memory
zeroMeanXNorm = zeros(nStates,1); % zero mean for additive noise

tEvents = feedInfoNorm(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFunNorm = @(t,XNorm) fNorm(XNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
    % create zero-mean normally distributed process noise for each sigma point:
    normalNoiseMatX = mvnrnd(zeroMeanXNorm,QNorm,nSigmaPointsAug)';
    for k = 1:nSigmaPointsAug
        [~,XTUSolNorm] = ode15s(odeFunNorm,tEval,sigmaXInitNorm(1:nStates,k));
        sigmaXPropNormNom(:,k) = XTUSolNorm(end,:)';     % nominal value (without noise)
    end 
    % add normally-distributed process noise acc. to Q (zero-mean):
    sigmaXPropNorm = sigmaXPropNormNom + normalNoiseMatX;

% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1;     
    XAtBeginOfIntNorm = sigmaXInitNorm(1:nStates,:);   % Startwerte für erstes Intervall (wird nicht mehr überschrieben)
    XAtEndOfIntNormNom = nan(size(XAtBeginOfIntNorm)); % allocate memory
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(m,2);
        xInCurrNorm = feedInfoNorm(m,3:end)';   % current inlet concentrations
        tEval = [tOverall(m), tOverall(m+1)];
        odeFunNorm = @(t,XNorm) fNorm(XNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum);
        % create zero-mean normally distributed process noise for each sigma point:
        normalNoiseMatX = mvnrnd(zeroMeanXNorm,QNorm,nSigmaPointsAug)';
        for kk = 1:nSigmaPointsAug
            [~,XTUSolNorm] = ode15s(odeFunNorm,tEval,XAtBeginOfIntNorm(:,kk));
            XAtEndOfIntNormNom(:,kk) = XTUSolNorm(end,:)';  % nominal value (without system noise)
        end
        XAtBeginOfIntNorm = XAtEndOfIntNormNom; % overwrite for next interval
    end
    % add gaussian, zero-mean process noise: 
    sigmaXPropNorm = XAtEndOfIntNormNom + normalNoiseMatX;    
end

% if any propagated sigma points violate constraints, apply clipping: 
if any(any(sigmaXPropNorm < 0))
    sigmaXPropNorm(sigmaXPropNorm < 0) = 0; 
    counterSigmaPropNorm = counterSigmaPropNorm + 1;
end

%% Aggregate Sigma Points to Priors for x and P
xMinusNorm = sum(Wx.*sigmaXPropNorm,2);  % state prior

% if any state priors violate constraints, apply clipping: 
if any(any(xMinusNorm < 0))
    xMinusNorm(xMinusNorm < 0) = 0; 
    counterXNorm = counterXNorm + 1;
end

% aggregate state error cov. matrix P:
diffXPriorFromSigmaNorm = sigmaXPropNorm - xMinusNorm; 
PMinusNorm = Wc.*diffXPriorFromSigmaNorm*diffXPriorFromSigmaNorm'; 

%------------------Measurement Update (MU)---------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them:
YNom = nan(q,nSigmaPointsAug);    % allocate memory
zeroMeanYNorm = zeros(q,1);        % zero mean for additive noise
% create zero-mean gaussian measurement noise for each sigma point:
normalNoiseMatY = mvnrnd(zeroMeanYNorm,RNorm,nSigmaPointsAug)';
for mm = 1:nSigmaPointsAug
    YNom(:,mm) = gNorm(sigmaXPropNorm(:,mm),c,TxNum,TyNum); 
end
% add gaussian, zero-mean process noise: 
YNorm = YNom + normalNoiseMatY; 
% aggregate outputs of sigma points in overall output:
yAggregatedNorm = sum(Wx.*YNorm,2);

%% 4. compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputsNorm = YNom - yAggregatedNorm; 
PyyNorm = Wc.*diffYFromSigmaOutputsNorm*diffYFromSigmaOutputsNorm'; % Kolas, Tab. 7/8

% compute cross covariance matrix states/measurements:
PxyNorm = Wc.*diffXPriorFromSigmaNorm*diffYFromSigmaOutputsNorm'; 

PyyNormInv = PyyNorm\eye(q);     % efficient least squares
KNorm = PxyNorm*PyyNormInv; 

%% update propagated sigma points individually 
% use the alternative formulation of Kolas 2009, eq. (23):
yMeasNorm = yMeas'./TyNum;   % normalized measurements
sigmaXNorm = sigmaXPropNorm + KNorm*(repmat(yMeasNorm,1,nSigmaPointsAug) - YNom);

% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaXNorm < 0))
    sigmaXNorm(sigmaXNorm < 0) = 0;
    counterSigmaXNorm = counterSigmaXNorm + 1; 
end

%% compute posteriors:
xPlusNorm = sum(Wx.*sigmaXNorm,2); 

% only for comparison: 
KvNorm = KNorm*(yMeasNorm - yAggregatedNorm);
xPlusNormvdM = xMinusNorm + KvNorm; % standard formulation of vdMerwe, normalized
disp(['max. Abweichung xPlus (fully add.):', num2str(max(abs(xPlusNormvdM - xPlusNorm)))])

% mind that Kolas proved the fully augmented case. For additive noise, the
% measurement update of P-Matrix must be slightly adapted:
diffxPlusFromSigmaXNorm = sigmaXNorm - xPlusNorm; 
PPlusKolasFullyAugmentedCaseNorm = Wc.*diffxPlusFromSigmaXNorm*diffxPlusFromSigmaXNorm'; 
PPlusTempNorm = PPlusKolasFullyAugmentedCaseNorm + KNorm*RNorm*KNorm'; % adapt formula for augmented process noise case (normalized)!

% only for comparison: 
PPlusTempNormvdM = PMinusNorm - KNorm*PyyNorm*KNorm'; 
disp(['max. Abweichung PPlus (fully add.): ', ...
      num2str(max(max(abs(PPlusTempNormvdM - PPlusTempNorm))))])

PPlusNorm = 1/2*(PPlusTempNorm + PPlusTempNorm');   % make sure PPlus is symmetric!
disp(['sum of PPlus diagonal (fully aug.): ', num2str(sum(diag(PPlusNorm)))])

end