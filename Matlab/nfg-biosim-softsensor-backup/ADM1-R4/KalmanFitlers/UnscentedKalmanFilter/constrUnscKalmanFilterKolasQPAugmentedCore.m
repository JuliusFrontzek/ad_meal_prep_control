%% Version
% (R2022b) Update 6
% Erstelldatum: 15.10.2023
% Autor: Simon Hellmann

function [xPlus,PPlus,nIter] = constrUnscKalmanFilterKolasQPAugmentedCore(xOld,POld, ...
                                    tSpan,feedInfo,yMeas,params,Q,R,f,g)

% compute time and measurement update of constrained QP-UKF acc. to  Kolas 
% et al. (2009), Tab. 10 for augmented process noise and QP formulation (abs. coordinates).
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

% global counterSigmaInit
% global counterSigmaProp
% global counterX
% global counterSigmaXcUKF

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
xOldAug = [xOld;zeros(nStates,1)]; 
POldAug = blkdiag(POld,Q);   % (2*nStates, 2*nStates)

nStatesAug = numel(xOldAug); 
nSigmaPointsAug = 2*(nStatesAug) + 1;   % # sigma points with augmentation

%% 1. Time Update (TU)

% re-define scaling parameters and weights for fully augmented case: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.0;  % leichte Abweichung zu Kolas (er nimmt 0)
% this creates a false scaling:
% lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
% gamma = sqrt(nStatesAug + lambda); % scaling parameter
% this creates the correct scaling:
lambda = alpha^2*(nStates + kappa) - nStates; 
gamma = sqrt(nStates + lambda); % scaling parameter
% gamma = 1;  % XY just to check

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

tEvents = feedInfo(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% additive process noise to be added after propagation:
addNoiseOnSigmapointsXMat = sigmaXInit(nStates+1:2*nStates,:); 

% we can only perform integration when feeding is constant!
% Fall a: konst. F�tterung w�hrend gesamtem Messintervalls (keine �nderung)
if isempty(tRelEvents)
    feedVolFlow = feedInfo(2); 
    xInCurr = feedInfo(3:end)';  % current inlet concentrations
    tEval = tSpan;
    odeFun = @(t,X) f(X,feedVolFlow,xInCurr,th,c,a); 
    for k = 1:nSigmaPointsAug
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInit(1:nStates,k));
        sigmaXPropNom(:,k) = XTUSol(end,:)';     % nominal value (without noise)
    end 
    % add effect of process noise to sigma points (Kolas, Tab. 8, Line 2):
    sigmaXProp = sigmaXPropNom + addNoiseOnSigmapointsXMat;

% Fall b: ver�nderliche F�tterung w�hrend Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus F�tterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtBeginOfInt = sigmaXInit(1:nStates,:);   % Startwerte f�r erstes Intervall (wird nicht mehr �berschrieben)
    XAtEndOfIntNom = nan(size(XAtBeginOfInt)); % allocate memory
    % integriere Intervall-weise, sodass w�hrend der Intervalle konstante
    % F�tterungen herrschen:
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
    % add effect of process noise to sigma points (Kolas, Tab. 8, Line 2):
    sigmaXProp = XAtEndOfIntNom + addNoiseOnSigmapointsXMat;    
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
Y = nan(q,nSigmaPointsAug);    % allocate memory
for mm = 1:nSigmaPointsAug
    Y(:,mm) = g(sigmaXProp(:,mm)); 
end
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
    [HMat,fTranspose] = computeQPCostFunMatrices(sigmaXProp(:,k),yMeas,R,PMinus); 
    x0QP = sigmaXProp(:,k); % initial vector for optimization
    [sigmaXOpt(:,k),fval,exitflag,output] = quadprog(HMat,fTranspose,A,b,[],[],[],[],x0QP,options); 
end 

% % this clipping should no longer be required thanks to optimization:
% % if updated sigma points violate constraints, apply clipping: 
% if any(any(sigmaXOpt < 0))
%     sigmaXOpt(sigmaXOpt < 0) = 0;
%     counterSigmaXcUKF = counterSigmaXcUKF + 1;
% end

%% compute posteriors:
xPlus = sum(Wx.*sigmaXOpt,2); 

diffxPlusFromSigmaX = sigmaXOpt - xPlus; 
PPlusKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX';
PPlusKolasAugmented = PPlusKolasFullyAugmented; % acc. to Vachhani (2006)

PPlus = 0.5*(PPlusKolasAugmented + PPlusKolasAugmented'); % ensure symmetry

%% return # iterations 
% but only if someone explicitly asks for them when calling this function:
nIter = output.iterations; 
if nargout < 3 
    nIter = [];
end

end