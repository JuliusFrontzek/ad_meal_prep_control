%% Version
% (R2022b) Update 5
% Erstelldatum: 16.09.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = my_UKF_augmented(xOld,POld,u,yMeas,tSpan,p,Q,R)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) for augmented system noise (Table 6) (no clipping)

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p - parameter vector for model the UKF is working with
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise

% if xOld contains negative concentrations, apply clipping: 
if any(xOld<0)
    xOld(xOld < 0) = 0; 
    counterX = counterX + 1;
end 

nStates = numel(xOld); 
q = 2; 

% augment x and P: 
xOldAug = [xOld;zeros(nStates,1)]; 
POldAug = blkdiag(POld,Q);   % (nStates+q, nStates+q)

nStatesAug = numel(xOldAug); 
% nSigmaPointsNom = 2*nStates + 1;        % nominal # sigma points (without augmentation)
nSigmaPointsAug = 2*(nStatesAug) + 1;   % # sigma points with augmentation

%% 1. Time Update (TU)

% re-define scaling parameters and weights for fully augmented case: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
% lambda = alpha^2*(nStates + kappa) - nStates; 
lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
% gamma = sqrt(nStates + lambda); % scaling parameter
gamma = sqrt(nStatesAug + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
% Wx0 = lambda/(nStates + lambda); 
Wx0 = lambda/(nStatesAug + lambda); 
% Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
% Wi = 1/(2*(nStates + lambda)); 
Wi = 1/(2*(nStatesAug + lambda)); 
% Wx = [Wx0, repmat(Wi,1,nSigmaPointsNom-1)]; % for state aggregation
% Wc = [Wc0, repmat(Wi,1,nSigmaPointsNom-1)]; % for covariance aggregation
Wx = [Wx0, repmat(Wi,1,nSigmaPointsAug-1)];
Wc = [Wc0, repmat(Wi,1,nSigmaPointsAug-1)];

% ensure weights are plausble:
if (sum(Wx) < 0.9) | (sum(Wx) > 1.1)
    disp('Wx ist nicht 1 in Summe!')
end
if (sum(Wc) < 0.9) | (sum(Wc) > 1.1)
    disp('Wc ist nicht 1 in Summe!')
end

%% 1.1) Choose Sigma Points
sqrtPOld = schol(POldAug);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOldAug, repmat(xOldAug,1,nStatesAug) + gamma*sqrtPOld, ...
                       repmat(xOldAug,1,nStatesAug) - gamma*sqrtPOld]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInit < 0)) 
    sigmaXInit(sigmaXInit < 0) = 0; 
    counterSigmaInit = counterSigmaInit + 1;
end

%% 1.2) Propagate all Sigma Points through system ODEs
sigmaXProp = nan(nStates, nSigmaPointsAug); % allocate memory
normalNoiseMat = nan(size(sigmaXProp));     % allocate memory

% integriere Sytemverhalten für alle Sigmapunkte im Interval t_span:
odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
for k = 1:nSigmaPointsAug
    [~,XTUSol] = ode45(odeFun,tSpan,sigmaXInit(1:nStates,k));
    sigmaXPropNom = XTUSol(end,:)';
    
    % add normally-distributed process noise acc. to Q (zero-mean):
    zeroMean = zeros(nStates,1);
    normalNoise = mvnrnd(zeroMean,Q,1)';
    normalNoiseMat(:,k) = normalNoise;
    sigmaXProp(:,k) = sigmaXPropNom + normalNoise;
end 

% % draw noise matrix: 
% plot(normalNoiseMat(1,:),normalNoiseMat(2,:),'+');

% if any propagated sigma points violate constraints, apply clipping: 
if any(any(sigmaXProp < 0))
    sigmaXProp(sigmaXProp < 0) = 0; 
    counterSigmaProp = counterSigmaProp + 1;
end

%% 1.3) Aggregate Sigma Points to Priors for x and P
% xMinus = sum(Wx.*sigmaXProp(:,1:nSigmaPointsNom),2);  % state prior
xMinus = sum(Wx.*sigmaXProp,2);  % state prior

% if any state priors violate constraints, apply clipping:
if any(any(xMinus < 0))
    xMinus(xMinus < 0) = 0; 
    counterX = counterX + 1; 
end

% aggregate state error cov. matrix P:
% diffXPriorFromSigma = sigmaXProp(:,1:nSigmaPointsNom) - xMinus; 
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; % adapted for additive noise case acc. to Kolas, Tab. 5

%% 2. Measurement Update (MU)

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPointsAug);    % allocate memory
for mm = 1:nSigmaPointsAug
    % hier den richtigen Aufruf der Messgleichung hinzufügen!
    Y(:,mm) = messgleichung(sigmaXProp(:,mm)); 
end
%% 2.2) aggregate outputs of sigma points in overall output:
% yAggregated = sum(Wx.*Y(:,1:nSigmaPointsNom),2);
yAggregated = sum(Wx.*Y,2);

%% 2.3) compute Kalman Gain

% compute cov. matrix of output Pyy:
% diffYFromSigmaOutputs = Y(:,1:nSigmaPointsNom) - yAggregated; 
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs' + R; % Kolas, Tab. 6

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

PyyInv = Pyy\eye(q);     % efficient least squares
K = Pxy*PyyInv; 

%% 2.4) update propagated sigma points individually 
% use alternative formulation of Kolas 2009, eq. (23):
sigmaX = sigmaXProp + K*(repmat(yMeas,1,nSigmaPointsAug) - Y);

% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaX < 0))
    sigmaX(sigmaX < 0) = 0;
    counterSigmaX = counterSigmaX + 1;
end

%% 2.5) compute posteriors:
% xPlus = sum(Wx.*sigmaX(:,1:nSigmaPointsNom),2); 
xPlus = sum(Wx.*sigmaX,2); 

% only for comparison: 
Kv = K*(yMeas - yAggregated);
xPlusvdM = xMinus + Kv; % standard formulation of vdMerwe
disp(['max. Abweichung xPlus (aug.):', num2str(max(abs(xPlusvdM - xPlus)))])

% Kolas (2009), Table 8:
% diffxPlusFromSigmaX = sigmaX(:,1:nSigmaPointsNom) - xPlus; 
diffxPlusFromSigmaX = sigmaX - xPlus; 
PPlusReformulatedKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX'; 
PPlusReformulatedKolasAugmented = PPlusReformulatedKolasFullyAugmented + K*R*K'; % + Q in add. noise case

% % only for comparison: 
PPlusTempvdM = PMinus - K*Pyy*K'; 
PPlusvdM = 1/2*(PPlusTempvdM + PPlusTempvdM');  % regularization
disp(['max. Abweichung PPlus (aug.): ', ...
      num2str(max(max(abs(PPlusvdM - PPlusReformulatedKolasAugmented))))])

% make sure PPlus is symmetric:
PPlus = 1/2*(PPlusReformulatedKolasAugmented + PPlusReformulatedKolasAugmented');   
% show potential divergence/falling asleep of P-Matrix live:
disp(['sum of PPlus diagonal (aug.): ', num2str(sum(diag(PPlus)))])

end