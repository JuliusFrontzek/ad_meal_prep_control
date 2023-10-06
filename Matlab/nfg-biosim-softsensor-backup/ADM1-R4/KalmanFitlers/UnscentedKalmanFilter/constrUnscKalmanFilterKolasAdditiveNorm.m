%% Version
% (R2022b) Update 5
% Erstelldatum: 03.10.2023
% Autor: Simon Hellmann

function [xPlusNorm,PPlusNorm] = constrUnscKalmanFilterKolasAdditiveNorm(xOldNorm,POldNorm, ...
                                    tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm, ...
                                    fNorm,gNorm,TxNum,TyNum,TuNum)

% compute time and measurement update of constrained UKF acc. to  Kolas et al. 
% (2009), Tab. 10 but for additive noise (normalized).
% Note: 3 differences between add. noise and fully augm. case:
% 1. augmentation, 
% 2. computation of PMinus, which includes +Q (add. noise case) or doesnt 
% (augmented system noise and fully augmented case)
% 3. computation of posterior of P: add + K*R*K' + Q compared with fully
% augmented case!

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

global counterSigmaInitNorm
global counterSigmaPropNorm
global counterSigmaXNorm
global counterXNorm
global counterSigmaXcUKF

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
nSigmaPoints = 2*nStates + 1;

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
% beta = 0;
beta = 2;   % for Gaussian prior (Diss vdM, S.56) 
% kappa = 3 - nStates;  % acc. to Julier & Uhlmann
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
gamma = sqrt(nStates + lambda); % scaling parameter
% gamma = 0.2;  % XY just to check

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
    counterSigmaInitNorm = counterSigmaInitNorm + 1;
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
    XAtBeginOfIntNorm = sigmaXInitNorm;    % Startwerte für erstes Intervall (wird nicht mehr überschrieben)
    XAtEndOfIntNorm = nan(size(sigmaXInitNorm)); % allocate memory
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(m,2);
        xInCurrNorm = feedInfoNorm(m,3:end)';   % current inlet concentrations
        tEval = [tOverall(m), tOverall(m+1)];
        odeFunNorm = @(t,XNorm) fNorm(XNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
        for kk = 1:nSigmaPoints
            [~,XTUSolNorm] = ode15s(odeFunNorm,tEval,XAtBeginOfIntNorm(:,kk));
            XAtEndOfIntNorm(:,kk) = XTUSolNorm(end,:)';
        end
        XAtBeginOfIntNorm = XAtEndOfIntNorm; % overwrite for next interval
    end
    sigmaXPropNorm = XAtEndOfIntNorm;
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
PMinusNorm = Wc.*diffXPriorFromSigmaNorm*diffXPriorFromSigmaNorm' + QNorm; % adapted for additive noise case acc. to Kolas, Tab. 5

%% ------------------Measurement Update (MU)------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them
YNorm = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    YNorm(:,mm) = gNorm(sigmaXPropNorm(:,mm),c,TxNum,TyNum); 
end

% aggregate outputs of sigma points in overall output:
yAggregatedNorm = sum(Wx.*YNorm,2);

%% run constrained optimization to determine sigmaX
% consider inequalities acc. to fmincon documentation: allow only positive 
% state values (=concentrations):
% XY Achtung: bisher noch keine Wasser-Konzentrationen > 1000 abgefangen, die
% für TS- und oTS-Berechnung Probleme machen
A = -eye(nStates); 
b = zeros(nStates,1);  
sigmaXOptNorm = nan(nStates,nSigmaPoints);    % allocate memory
yMeasNorm = yMeas'./TyNum;   % normalized measurements
% optimize all updated (normalized) sigma points: 
for k = 1:nSigmaPoints
    ukfCostFunNorm = @(sigmaXNorm) evaluateCUKFCostFunNorm(sigmaXNorm,sigmaXPropNorm(:,k), ...
            yMeasNorm,RNorm,PMinusNorm,gNorm,c,TxNum,TyNum); 
    % choose the old sigmaXPropNorm as initial value for optimization:
    sigmaXOptNorm(:,k) = fmincon(ukfCostFunNorm,sigmaXPropNorm(:,k),A,b); 
end 

% this clipping should no longer be required thanks to optimization:
% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaXOptNorm < 0))
    sigmaXOptNorm(sigmaXOptNorm < 0) = 0;
    counterSigmaXcUKF = counterSigmaXcUKF + 1;
end

%% compute posteriors:
xPlusNorm = sum(Wx.*sigmaXOptNorm,2); 

% mind: Kolas considers fully augmented case, so computation of
% posteriors must be adapted. I believe Vachhani (who also considers
% additive noise case) is wrong!
diffxPlusFromSigmaXNorm = sigmaXOptNorm - xPlusNorm; 

% I think just as for the unconstrained case, we need to adapt computation
% of the posterior of P by incorporating K, Q and R, so calculate K (normalized):
% compute cov. matrix of output Pyy:
diffYFromSigmaOutputsNorm = YNorm - yAggregatedNorm; 
PyyNorm = Wc.*diffYFromSigmaOutputsNorm*diffYFromSigmaOutputsNorm' + RNorm;

% compute cross covariance matrix states/measurements:
PxyNorm = Wc.*diffXPriorFromSigmaNorm*diffYFromSigmaOutputsNorm'; 

PyyNormInv = PyyNorm\eye(q);     % efficient least squares
KNorm = PxyNorm*PyyNormInv; 

PPlusFullyAugmentedCaseNorm = Wc.*diffxPlusFromSigmaXNorm*diffxPlusFromSigmaXNorm'; 
% adapt Kolas (2009) just like for the unconstrained case:
PPlusTempNorm = PPlusFullyAugmentedCaseNorm + KNorm*RNorm*KNorm' + QNorm; % different formula for additive noise case (normalized)!

% make sure PPlus is symmetric:
PPlusNorm = 1/2*(PPlusTempNorm + PPlusTempNorm');
disp(['sum of PPlusNorm diagonal (cUKF-add.): ', num2str(sum(diag(PPlusNorm)))])

end