%% Version
% (R2022b) Update 6
% Erstelldatum: 10.10.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = constrUnscKalmanFilterKolasQPFullyAugmentedCore(xOld,POld, ...
                                    tSpan,feedInfo,yMeas,params,Q,R,f,g)

% compute time and measurement update of constrained QP-UKF acc. to  Kolas 
% et al. (2009), Tab. 10 for fully augmented noise and QP formulation (abs. coordinates).
% assume measurements of ADM1-R4-Core

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector (non normalized)
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise
% f - function handle of ODEs of system equations
% g - function handle of output equations 

global counterSigmaInit
global counterSigmaProp
% global counterSigmaX
global counterX
global counterSigmaXcUKF

% extract constant parameters out of struct: 
th = params.th; 
c = params.c; 
a = params.a;

% % if xOld contains negative concentrations, apply clipping: 
% if any(xOld<0)
%     xOld(xOld < 0) = 0; 
%     counterX = counterX + 1;
% end 

nStates = numel(xOld); 
q = numel(yMeas); 

% augment x and P: 
xOldAug = [xOld;zeros(nStates,1);zeros(q,1)]; 
POldAug = blkdiag(POld,Q,R);   % (2*nStates+q, 2*nStates+q)

nStatesAug = numel(xOldAug); 
nSigmaPointsAug = 2*(nStatesAug) + 1;   % # sigma points with augmentation

%% 1. Time Update (TU)

% re-define scaling parameters and weights for fully augmented case: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.0;  % leichte Abweichung zu Kolas (er nimmt 0)
% lambda = alpha^2*(nStates + kappa) - nStates; 
lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
% gamma = sqrt(nStates + lambda); % scaling parameter
gamma = sqrt(nStatesAug + lambda); % scaling parameter
% gamma = 0.5;  % XY just to check

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
% Wx0 = lambda/(nStates + lambda); 
Wx0 = lambda/(nStatesAug + lambda); 
Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStatesAug + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPointsAug-1)];
Wc = [Wc0, repmat(Wi,1,nSigmaPointsAug-1)];

% ensure weights are plausble:
if (sum(Wx) < 0.9) | (sum(Wx) > 1.1)
    disp('Wx ist nicht 1 in Summe!')
end
% if (sum(Wc) < 0.9) | (sum(Wc) > 1.1)
%     disp('Wc ist nicht 1 in Summe!')
% end

%% 1.1) Choose Sigma Points
sqrtPOld = schol(POldAug);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOldAug, repmat(xOldAug,1,nStatesAug) + gamma*sqrtPOld, ...
                       repmat(xOldAug,1,nStatesAug) - gamma*sqrtPOld]; 

% % Apply clipping to negative Sigma Points: 
% if any(any(sigmaXInit < 0)) 
%     sigmaXInit(sigmaXInit < 0) = 0; 
%     counterSigmaInit = counterSigmaInit + 1;
% end

%% 1.2) Propagate all Sigma Points through system ODEs
sigmaXPropNom = nan(nStates, nSigmaPointsAug); % allocate memory
zeroMeanX = zeros(nStates,1); % zero mean for additive noise

tEvents = feedInfo(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVolFlow = feedInfo(2); 
    xInCurr = feedInfo(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFun = @(t,X) f(X,feedVolFlow,xInCurr,th,c,a); 
    for k = 1:nSigmaPointsAug
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInit(1:nStates,k));
        sigmaXPropNom(:,k) = XTUSol(end,:)';     % nominal value (without noise)
    end 
    % create zero-mean normally distributed process noise for each sigma point:
    normalNoiseMatX = mvnrnd(zeroMeanX,Q,nSigmaPointsAug)';
    % add normally-distributed process noise acc. to Q (zero-mean):
    sigmaXProp = sigmaXPropNom + normalNoiseMatX;

% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtBeginOfInt = sigmaXInit(1:nStates,:);   % Startwerte für erstes Intervall (wird nicht mehr überschrieben)
    XAtEndOfIntNom = nan(size(XAtBeginOfInt)); % allocate memory
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVolFlow = feedInfo(m,2);
        xInCurr = feedInfo(m,3:end)';   % current inlet concentrations
        tEval = [tOverall(m), tOverall(m+1)];
        odeFun = @(t,X) f(X,feedVolFlow,xInCurr,th,c,a);
        for kk = 1:nSigmaPointsAug
            [~,XTUSol] = ode15s(odeFun,tEval,XAtBeginOfInt(:,kk));
            XAtEndOfIntNom(:,kk) = XTUSol(end,:)';  % nominal value (without system noise)
        end
        XAtBeginOfInt = XAtEndOfIntNom; % overwrite for next interval
    end
    % create zero-mean normally distributed process noise for each sigma point:
    normalNoiseMatX = mvnrnd(zeroMeanX,Q,nSigmaPointsAug)';
    % add gaussian, zero-mean process noise: 
    sigmaXProp = XAtEndOfIntNom + normalNoiseMatX;    
end

% % draw noise matrix: 
% plot(normalNoiseMat(1,:),normalNoiseMat(2,:),'+');

% % if any propagated sigma points violate constraints, apply clipping: 
% if any(any(sigmaXProp < 0))
%     sigmaXProp(sigmaXProp < 0) = 0; 
%     counterSigmaProp = counterSigmaProp + 1;
% end

%% 1.3) Aggregate Sigma Points to Priors for x and P
% xMinus = sum(Wx.*sigmaXProp(:,1:nSigmaPointsNom),2);  % state prior
xMinus = sum(Wx.*sigmaXProp,2);  % state prior

% % if any state priors violate constraints, apply clipping:
% if any(any(xMinus < 0))
%     xMinus(xMinus < 0) = 0; 
%     counterX = counterX + 1; 
% end

% aggregate state error cov. matrix P:
% diffXPriorFromSigma = sigmaXProp(:,1:nSigmaPointsNom) - xMinus; 
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; % adapted for fully augmented noise case acc. to Kolas, Tab. 10

%% 2. Measurement Update (MU)

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
YNom = nan(q,nSigmaPointsAug);    % allocate memory
for mm = 1:nSigmaPointsAug
    YNom(:,mm) = g(sigmaXProp(:,mm)); 
end
% create zero-mean normally distributed measurement noise for each sigma point:
zeroMeanY = zeros(q,1);
normalNoiseMatY = mvnrnd(zeroMeanY,R,nSigmaPointsAug)';
% add normally-distributed process noise acc. to Q (zero-mean):
Y = YNom + normalNoiseMatY;

% 2.2) aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

%% run constrained optimization to determine sigmaX
% consider inequalities acc. to fmincon documentation: allow only positive 
% state values (=concentrations):
A = -eye(nStates); 
b = zeros(nStates,1);  
sigmaXOpt = nan(nStates,nSigmaPointsAug);    % allocate memory
options = optimoptions('quadprog','Display','none'); % suppress command window output 
% optimize all updated sigma points: 
for k = 1:nSigmaPointsAug
    % compute matrices H and f for QP-solver quadprog:
    PMinusX = PMinus(1:nStates,1:nStates);  % extract only x-part of fully augmented P-Matrix
    [HMat,fTranspose] = computeQPCostFunMatrices(sigmaXProp(:,k),yMeas,R,PMinusX); 
    x0QP = sigmaXProp(:,k); % initial vector for optimization
    sigmaXOpt(:,k) = quadprog(HMat,fTranspose,A,b,[],[],[],[],x0QP,options); 
end 

% this clipping should no longer be required thanks to optimization:
% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaXOpt < 0))
    sigmaXOpt(sigmaXOpt < 0) = 0;
    counterSigmaXcUKF = counterSigmaXcUKF + 1;
end

%% compute posteriors:
xPlus = sum(Wx.*sigmaXOpt,2); 

diffxPlusFromSigmaX = sigmaXOpt - xPlus; 
PPlusKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX';

PPlus = 0.5*(PPlusKolasFullyAugmented + PPlusKolasFullyAugmented'); % ensure symmetry

end