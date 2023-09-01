%% Version
% (R2022b) Update 5
% Erstelldatum: 31.08.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = unscKalmanFilterKolasAdditiveNorm(xOldNorm,POldNorm,tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,TxNum,TyNum,TuNum)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) with clipping wherever possible (normalized)

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX

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
    counterX = counterX + 1;
end 

nStates = numel(xOldNorm); 
q = numel(yMeas); 
nSigmaPoints = 2*nStates + 1;

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
gamma = sqrt(nStates + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% Choose Sigma Points
sqrtPOldNorm = schol(POldNorm);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInitNorm = [xOldNorm, repmat(xOldNorm,1,nStates) + gamma*sqrtPOldNorm, ...
                            repmat(xOldNorm,1,nStates) - gamma*sqrtPOldNorm]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInitNorm < 0))
    sigmaXInitNorm(sigmaXInitNorm < 0) = 0; 
    counterSigmaInit = counterSigmaInit + 1;
end

%% Propagate Sigma Points
sigmaXPropNorm = nan(nStates, nSigmaPoints); % allocate memory

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
if any(any(sigmaXPropNorm < 0))
    sigmaXPropNorm(sigmaXPropNorm < 0) = 0; 
    counterSigmaProp = counterSigmaProp + 1;
end

%% Aggregate Sigma Points to Priors for x and P
xMinusNorm = sum(Wx.*sigmaXPropNorm,2);  % state prior

% if any state priors violate constraints, apply clipping: 
if any(any(xMinusNorm < 0))
    xMinusNorm(xMinusNorm < 0) = 0; 
    counterX = counterX + 1;
end

% aggregate state error cov. matrix P:
diffXPriorFromSigmaNorm = sigmaXPropNorm - xMinusNorm; 
PMinusNorm = Wc.*diffXPriorFromSigmaNorm*diffXPriorFromSigmaNorm' + QNorm; 

%------------------Measurement Update (MU)---------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them:
YNorm = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    YNorm(:,mm) = gNorm(sigmaXPropNorm(:,mm),c,TxNum,TyNum); 
end
% aggregate outputs of sigma points in overall output:
yAggregatedNorm = sum(Wx.*YNorm,2);

%% 4. compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputsNorm = YNorm - yAggregatedNorm; 
PyyNorm = Wc.*diffYFromSigmaOutputsNorm*diffYFromSigmaOutputsNorm' + RNorm;

% compute cross covariance matrix states/measurements:
PxyNorm = Wc.*diffXPriorFromSigmaNorm*diffYFromSigmaOutputsNorm'; 

PyyNormInv = PyyNorm\eye(q);     % efficient least squares
KNorm = PxyNorm*PyyNormInv; 

%% update propagated sigma points individually 
% use the alternative formulation of Kolas 2009, eq. (23):
yMeasNorm = yMeas'./TyNum;   % normalized measurements
sigmaXNorm = sigmaXPropNorm + KNorm*(repmat(yMeasNorm,1,nSigmaPoints) - YNorm);

% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaXNorm < 0))
    sigmaXNorm(sigmaXNorm < 0) = 0;
    counterSigmaX = counterSigmaX + 1; 
end

%% compute posteriors:
xPlusNorm = sum(Wx.*sigmaXNorm,2); 

% only for comparison: 
KvNorm = KNorm*(yMeasNorm' - yAggregatedNorm);
xPlusNormvdM = xMinusNorm + KvNorm; % standard formulation of vdMerwe, normalized

diffxPlusFromSigmaXNorm = sigmaXNorm - xPlusNorm; 
PPlusAugmentedCaseNorm = Wc.*diffxPlusFromSigmaXNorm*diffxPlusFromSigmaXNorm'; 

PPlusTempNorm = PPlusAugmentedCaseNorm + QNorm + KNorm*RNorm*KNorm'; % different formula for additive noise case (normalized)!

% only for comparison: 
PPlusTempNormvdM = PMinusNorm - KNorm*PyyNorm*KNorm'; 

PPlusNorm = 1/2*(PPlusTempNorm + PPlusTempNorm');   % make sure PPlus is symmetric!
% sum(diag(PPlusNorm)) % show potential divergence of P-Matrix live

end