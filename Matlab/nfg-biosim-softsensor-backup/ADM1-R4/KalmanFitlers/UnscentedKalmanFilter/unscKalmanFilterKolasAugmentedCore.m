%% Version
% (R2022b) Update 5
% Erstelldatum: 06.10.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = unscKalmanFilterKolasAugmentedCore(xOld,POld, ...
                                    tSpan,feedInfo,yMeas,params,Q,R,f,g)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) with clipping wherever possible, but for 
% augmented process noise (Tab. 6), use abs. coordinates

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX

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
POldAug = blkdiag(POld,Q);   % (nStates+q, nStates+q)

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
% gamma = 0.05;  % XY just to check

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
% Wx0 = lambda/(nStates + lambda); 
Wx0 = lambda/(nStatesAug + lambda); 
% Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
% Wi = 1/(2*(nStates + lambda)); 
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
    % create zero-mean normally distributed process noise for each sigma point:
    normalNoiseMatX = mvnrnd(zeroMeanX,Q,nSigmaPointsAug)';
    for k = 1:nSigmaPointsAug
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInit(1:nStates,k));
        sigmaXPropNom(:,k) = XTUSol(end,:)';     % nominal value (without noise)
    end 
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
        % create zero-mean normally distributed process noise for each sigma point:
        normalNoiseMatX = mvnrnd(zeroMeanX,Q,nSigmaPointsAug)';
        for kk = 1:nSigmaPointsAug
            [~,XTUSol] = ode15s(odeFun,tEval,XAtBeginOfInt(:,kk));
            XAtEndOfIntNom(:,kk) = XTUSol(end,:)';  % nominal value (without system noise)
        end
        XAtBeginOfInt = XAtEndOfIntNom; % overwrite for next interval
    end
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
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; % adapted for augmented process noise case acc. to Kolas, Tab. 6

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

%% 2.3) compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs' + R; % Kolas, Tab. 6

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% PyyInv = Pyy\eye(q);     % efficient least squares
% K = Pxy*PyyInv; 
K = Pxy/Pyy; 

%% 2.4) update propagated sigma points individually 
% use alternative formulation of Kolas 2009, eq. (23):
sigmaX = sigmaXProp + K*(repmat(yMeas',1,nSigmaPointsAug) - Y);

% % if updated sigma points violate constraints, apply clipping: 
% if any(any(sigmaX < 0))
%     sigmaX(sigmaX < 0) = 0;
%     counterSigmaX = counterSigmaX + 1;
% end

%% 2.5) compute posteriors:
xPlus = sum(Wx.*sigmaX,2); 

% only for comparison: 
Kv = K*(yMeas' - yAggregated);
xPlusvdM = xMinus + Kv; % standard formulation of vdMerwe
disp(['max. Abweichung xPlus (aug.):', num2str(max(abs(xPlusvdM - xPlus)))])

% Kolas (2009), Table 8:
diffxPlusFromSigmaX = sigmaX - xPlus; 
PPlusReformulatedKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX';
PPlusReformulatedKolasAugmented = PPlusReformulatedKolasFullyAugmented + K*R*K'; % adapted Kolas' proof in appendix for augmented process noise
PPlusVachhaniTemp = PPlusReformulatedKolasFullyAugmented; 

% only for comparison: 
PPlusTempvdM = PMinus - K*Pyy*K'; 
PPlusvdM = 1/2*(PPlusTempvdM + PPlusTempvdM');  % regularization
disp(['max. Abweichung PPlus (aug.): ', ...
      num2str(max(max(abs(PPlusvdM - PPlusReformulatedKolasFullyAugmented))))])

% make sure PPlus is symmetric:
PPlus = 1/2*(PPlusReformulatedKolasAugmented + PPlusReformulatedKolasAugmented');   
% show potential divergence/falling asleep of P-Matrix live:
disp(['sum of PPlus diagonal (aug.): ', num2str(sum(diag(PPlus)))])

end