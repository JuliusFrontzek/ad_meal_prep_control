%% Version
% (R2022b) Update 5
% Erstelldatum: 08.10.2023
% Autor: Simon Hellmann

function [xPlus,SPlus] = my_SR_UKF_additive(xOld,SOld,u,yMeas,tSpan,p,SQ,SR)

% compute time and measurement update acc. to SR-UKF of vdMerwe (2001)
% for additive noise (without clipping)

% global counterSigmaInit
% global counterSigmaProp
% global counterSigmaX
% global counterX

% xPlus - new state estimate
% SPlus - upper cholesky factor of new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% SOld - upper cholesky factor of old state error covariance matrix 
% tSpan - time interval between old and new measurement (& state estimate)
% yMeas - latest measurement vector
% p - parameter vector for model the UKF is working with
% SQ - lower chol. factor of power spectral density matrix of process noise
% SR - lower chol. factor of covariance matrix of measurement noise

% % if xOld contains negative concentrations, apply clipping: 
% if any(xOld<0)
%     xOld(xOld < 0) = 0; 
%     counterX = counterX + 1;
% end 

nStates = numel(xOld); 
q = 2; 
nSigmaPoints = 2*nStates + 1;

%% 1. Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % std. value of matlab: 1E-3; vorher 1, Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.0;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
gamma = sqrt(nStates + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

% ensure weights are plausble:
if (sum(Wx) < 0.9) | (sum(Wx) > 1.1)
    disp('Wx ist nicht 1 in Summe!')
end
% if (sum(Wc) < 0.9) | (sum(Wc) > 1.1)
%     disp('Wc ist nicht 1 in Summe!')
% end

%% 1.1) Choose Sigma Points

sigmaXInit = [xOld, repmat(xOld,1,nStates) + gamma*SOld, ...
                    repmat(xOld,1,nStates) - gamma*SOld]; 

% % Apply clipping to negative Sigma Points: 
% if any(any(sigmaXInit < 0)) 
%     sigmaXInit(sigmaXInit < 0) = 0; 
%     counterSigmaInit = counterSigmaInit + 1;
% end

%% 1.2) Propagate all Sigma Points through system ODEs
sigmaXProp = nan(nStates, nSigmaPoints); % allocate memory

% integriere Sytemverhalten für alle Sigmapunkte im Interval t_span:
odeFun = @(t,x) my_bioprocess_ode(t,x,u,p);
for k = 1:nSigmaPoints
    [~,XTUSol] = ode45(odeFun,tSpan,sigmaXInit(:,k));
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

% perform cholupdate of PMinus (20 and 21 from vdMerwe 2001):
[~,SMinusTemp] = qr([sqrt(Wi)*(sigmaXProp(:,2:end) - xMinus), SQ]','econ');  
SMinus = cholupdate(SMinusTemp, sqrt(Wc0)*(sigmaXProp(:,1) - xMinus)); % upper triangular matrix

%% 2. Measurement Update (MU)
% optional: draw new sigma point acc. to a priori estimate xMinus and thus 
% overwrite sigmaXProp: 
% sigmaXProp = [xMinus, repmat(xMinus,1,nStates) + gamma*SMinus, ...
%                       repmat(xMinus,1,nStates) - gamma*SMinus]; 

% or omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    % hier den richtigen Aufruf der Messgleichung hinzufügen!
    Y(:,mm) = messgleichung(sigmaXProp(:,mm)); 
end
%% 2.2) aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

% perform cholupdate of PMinus (20 and 21 from vdMerwe 2001): 
[~,SyTemp] = qr([sqrt(Wi)*(Y(:,2:end) - yAggregated), SR]','econ');  
Sy = cholupdate(SyTemp, sqrt(Wc0)*(Y(:,1) - yAggregated)); % upper triangular matrix

% compute cross covariance matrix states/measurements:
% PxyMat = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% compute the aggregation through a recursive sum: 
PxyTemp = zeros(nStates,q); % initialization
for k = 1:nSigmaPoints
    PxyTemp = PxyTemp + Wc(k)*((sigmaXProp(:,k) - xMinus)*(Y(:,k) - yAggregated)');
end
Pxy = PxyTemp; 

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
diffXPriorFromSigma = sigmaXProp - xMinus; 
% compute cross covariance matrix states/measurements:
PxyMat = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% K = (Pxy/Sy')/Sy; % Kalman gain acc. to 27 of vdMerwe (2001)
K = Pxy/Sy/Sy';     % Kalman gain acc. to Arrese (2016): S is always the 
% lower triangular cholesky factor in vdMerwe. However, cholupdate in
% matlab only considers upper triangular factors and returns upper
% triangular factors as well, so transposition is necessary.

xPlus = xMinus + K*(yMeas - yAggregated); 
U = K*Sy';  % transposition because of same reason as for K

nColsU = size(U,2); % # columns of U
% conduct (29) of vdMerwe (2001):
SPlusTemp = SMinus; % initialization
for k = 1:nColsU
    SPlusTemp = cholupdate(SPlusTemp, U(:,k), '-'); 
end
SPlus = SPlusTemp; 

end