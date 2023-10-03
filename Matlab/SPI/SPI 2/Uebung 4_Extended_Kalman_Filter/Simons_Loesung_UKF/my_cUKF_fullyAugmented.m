%% Version
% (R2022b) Update 5
% Erstelldatum: 03.10.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = my_cUKF_fullyAugmented(xOld,POld,u,yMeas,tSpan,p,Q,R)

% compute time and measurement update of constrained UKF acc. to 
% Kolas et al. (2009), Tab. 10 (without clipping)

% global counterSigmaInit
% global counterSigmaProp
global counterSigmaXcUKF
% global counterX

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p - parameter vector for model the UKF is working with
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise

% % if xOld contains negative concentrations, apply clipping: 
% if any(xOld<0)
%     xOld(xOld < 0) = 0; 
%     counterX = counterX + 1;
% end 

nStates = numel(xOld); 
q = 2; 

% augment x and P: 
xOldAug = [xOld;zeros(nStates,1);zeros(q,1)]; 
POldAug = blkdiag(POld,Q,R);   % (2*nStates+q, 2*nStates+q)

nStatesAug = numel(xOldAug); 
nSigmaPointsAug = 2*(nStatesAug) + 1;   % # sigma points with augmentation

%% 1. Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
gamma = sqrt(nStatesAug + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStatesAug + lambda); 
Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStatesAug + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPointsAug-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPointsAug-1)]; % for covariance aggregation

%% 1.1) Choose Sigma Points
sqrtPOldAug = schol(POldAug);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOldAug,  repmat(xOldAug,1,nStatesAug) + gamma*sqrtPOldAug, ...
                        repmat(xOldAug,1,nStatesAug) - gamma*sqrtPOldAug]; 

% % Apply clipping to negative Sigma Points: 
% if any(any(sigmaXInit < 0)) 
%     sigmaXInit(sigmaXInit < 0) = 0; 
%     counterSigmaInit = counterSigmaInit + 1;
% end

%% 1.2) Propagate all Sigma Points through system ODEs
sigmaXProp = nan(nStates, nSigmaPointsAug); % allocate memory
normalNoiseMatX = nan(size(sigmaXProp));    % allocate memory

% integriere Sytemverhalten für alle Sigmapunkte im Interval t_span:
odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
for k = 1:nSigmaPointsAug
    [~,XTUSol] = ode45(odeFun,tSpan,sigmaXInit(1:nStates,k));
    sigmaXPropNom = XTUSol(end,:)'; % nominal propagation without noise

    % add normally-distributed process noise acc. to Q (zero-mean):
    zeroMeanX = zeros(nStates,1);
    normalNoiseX = mvnrnd(zeroMeanX,Q,1)';
    normalNoiseMatX(:,k) = normalNoiseX;
    sigmaXProp(:,k) = sigmaXPropNom + normalNoiseX;
end 

% sigmaXProp = zeros(size(sigmaXProp)); % Spaß für Terrance
% % draw noise matrix: 
% plot(normalNoiseMatX(1,:),normalNoiseMatX(2,:),'+');

% % if any propagated sigma points violate constraints, apply clipping: 
% if any(any(sigmaXProp < 0))
%     sigmaXProp(sigmaXProp < 0) = 0; 
%     counterSigmaProp = counterSigmaProp + 1;
% end

%% 1.3) Aggregate Sigma Points to Priors for x and P
xMinus = sum(Wx.*sigmaXProp,2);  % state prior

% % if any state priors violate constraints, apply clipping:
% if any(any(xMinus < 0))
%     xMinus(xMinus < 0) = 0; 
%     counterX = counterX + 1; 
% end

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; % fully augmented case acc. to Kolas, Tab. 10

%% 2. Measurement Update (MU)

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
% not needed anymore in NLP-formulation of cUKF. Only needed for
% QP-forulation:
% Y = nan(q,nSigmaPointsAug);    % allocate memory
% for mm = 1:nSigmaPointsAug
%     % hier den richtigen Aufruf der Messgleichung hinzufügen!
%     Y(:,mm) = messgleichung(sigmaXProp(:,mm)); 
% end

% aggregate outputs of sigma points in overall output:
% yAggregated = sum(Wx.*Y,2);

%% run constrained optimization to determine sigmaX
% consider inequalities acc. to fmincon documentation: allow only positive 
% state values (=concentrations):
A = -eye(nStates); 
b = zeros(nStates,1);  
sigmaXOpt = nan(nStates,nSigmaPointsAug);    % allocate memory
% optimize all updated sigma points: 
for k = 1:nSigmaPointsAug
    costFun = @(sigmaX) evaluateCUKFCostFun(sigmaX,sigmaXProp(:,k),yMeas,R,PMinus); 
    % choose the old sigmaXProp as initial value for optimization:
    sigmaXOpt(:,k) = fmincon(costFun,sigmaXProp(:,k),A,b); 
end 

% this clipping should no longer be required thanks to optimization:
% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaXOpt < 0))
    sigmaXOpt(sigmaXOpt < 0) = 0;
    counterSigmaXcUKF = counterSigmaXcUKF + 1;
end

%% 2.5) compute posteriors:
xPlus = sum(Wx.*sigmaXOpt,2); 

% mind: although Kolas considers fully augmented case, computation of
% posteriors remains the same, see also Vachhani, who also considers
% additive noise case:
diffxPlusFromSigmaX = sigmaXOpt - xPlus; 
PPlusTemp = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX';

% make sure PPlus is symmetric:
PPlus = 1/2*(PPlusTemp + PPlusTemp');   
disp(['sum of PPlus diagonal (cUKF.): ', num2str(sum(diag(PPlus)))])

end