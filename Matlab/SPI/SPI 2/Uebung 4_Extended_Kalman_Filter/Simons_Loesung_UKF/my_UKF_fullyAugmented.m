%% Version
% (R2022b) Update 5
% Erstelldatum: 16.09.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = my_UKF_fullyAugmented(xOld,POld,u,yMeas,tSpan,p,Q,R)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) for non-additive noise (Table 8) (no clipping)

% global counterSigmaInit
% global counterSigmaProp
% global counterSigmaX
% global counterX

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p_KF - parameter vector for model the UKF is working with
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
nSigmaPoints = 2*(nStatesAug) + 1;

%% 1. Time Update (TU)

% re-define scaling parameters and weights for fully augmented case: 
alpha = 1;  % Kolas 2009, (18)
beta = 0;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.0;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
lambda = alpha^2*(nStatesAug + kappa) - nStatesAug; 
gamma = sqrt(nStates + lambda); % scaling parameter
gamma = sqrt(nStatesAug + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
% Wx0 = lambda/(nStatesAug + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
% Wc0 = lambda/(nStatesAug + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
% Wi = 1/(2*(nStatesAug + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% 1.1) Choose Sigma Points
sqrtPOld = schol(POldAug);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOldAug, repmat(xOldAug,1,2*nStates+q) + gamma*sqrtPOld, ...
                       repmat(xOldAug,1,2*nStates+q) - gamma*sqrtPOld]; 

% % Apply clipping to negative Sigma Points: 
% if any(any(sigmaXInit < 0)) 
%     sigmaXInit(sigmaXInit < 0) = 0; 
%     counterSigmaInit = counterSigmaInit + 1;
% end

%% 1.2) Propagate all Sigma Points through system ODEs
sigmaXProp = nan(nStates, nSigmaPoints); % allocate memory

% integriere Sytemverhalten f�r alle Sigmapunkte im Interval t_span:
odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
for k = 1:nSigmaPoints
    [~,XTUSol] = ode45(odeFun,tSpan,sigmaXInit(1:nStates,k));
    sigmaXProp(:,k) = XTUSol(end,:)';
end 

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
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; % adapted for additive noise case acc. to Kolas, Tab. 5

%% 2. Measurement Update (MU)

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    % hier den richtigen Aufruf der Messgleichung hinzuf�gen!
    Y(:,mm) = messgleichung(sigmaXProp(:,mm)); 
end
%% 2.2) aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

%% 2.3) compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs';

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

PyyInv = Pyy\eye(q);     % efficient least squares
K = Pxy*PyyInv; 

%% 2.4) update propagated sigma points individually 
% use alternative formulation of Kolas 2009, eq. (23):
sigmaX = sigmaXProp + K*(repmat(yMeas,1,nSigmaPoints) - Y);

% % if updated sigma points violate constraints, apply clipping: 
% if any(any(sigmaX < 0))
%     sigmaX(sigmaX < 0) = 0;
%     counterSigmaX = counterSigmaX + 1;
% end

%% 2.5) compute posteriors:
xPlus = sum(Wx.*sigmaX,2); 

% % only for comparison: 
% Kv = K*(yMeas - yAggregated);
% xPlusvdM = xMinus + Kv; % standard formulation of vdMerwe
% disp('max. Abweichung xPlus:')
% max(abs(xPlusvdM - xPlus))

% Kolas (2009), Table 8:
diffxPlusFromSigmaX = sigmaX - xPlus; 
PPlusReformulatedKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX'; 

% % only for comparison: 
% PPlusTempvdM = PMinus - K*Pyy*K'; 
% disp('max. Abweichung PPlus:')
% max(max(abs(PPlusTempvdM - PPlusReformulatedKolasAdditive)))

% make sure PPlus is symmetric:
PPlus = 1/2*(PPlusReformulatedKolasFullyAugmented + PPlusReformulatedKolasFullyAugmented');   
disp('sum of PPlus diagonal:')
sum(diag(PPlus)) % show potential divergence/falling asleep of P-Matrix live

end