%% Version
% (R2022b) Update 5
% Erstelldatum: 16.09.2023
% Autor: Simon Hellmann

function [xPlus,PPlus,Kv] = my_UKF_additive(xOld,POld,u,yMeas,tSpan,p,Q,R)

% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) for additive noise (without clipping)

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
nSigmaPoints = 2*nStates + 1;

%% 1. Time Update (TU)

% define scaling parameters and weights: 
alpha = 1E-3;  % std. value of matlab: 1E-3; vorher 1, Kolas 2009, (18)
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
sqrtPOld = chol(POld);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOld, repmat(xOld,1,nStates) + gamma*sqrtPOld, ...
                    repmat(xOld,1,nStates) - gamma*sqrtPOld]; 

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

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinusMat = Wc.*diffXPriorFromSigma*diffXPriorFromSigma' + Q; % adapted for additive noise case acc. to Kolas, Tab. 5

% compute the aggregation through a recursive sum: 
PMinusTemp = zeros(nStates); % initialization
for k = 1:nSigmaPoints
    PMinusTemp = PMinusTemp + Wc(k)*((sigmaXProp(:,k) - xMinus)*(sigmaXProp(:,k) - xMinus)');
end
PMinus = PMinusTemp + Q; % identical result as before!

%% 2. Measurement Update (MU)

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% 2.1) Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    % hier den richtigen Aufruf der Messgleichung hinzufügen!
    Y(:,mm) = messgleichung(sigmaXProp(:,mm)); 
end
%% 2.2) aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

%% 2.3) compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
PyyMat = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs' + R;

% compute the aggregation through a recursive sum: 
PyyTemp = zeros(q); % initialization
for k = 1:nSigmaPoints
    PyyTemp = PyyTemp + Wc(k)*((Y(:,k) - yAggregated)*(Y(:,k) - yAggregated)');
end
Pyy = PyyTemp + R; 

% compute cross covariance matrix states/measurements:
PxyMat = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% compute the aggregation through a recursive sum: 
PxyTemp = zeros(nStates,q); % initialization
for k = 1:nSigmaPoints
    PxyTemp = PxyTemp + Wc(k)*((sigmaXProp(:,k) - xMinus)*(Y(:,k) - yAggregated)');
end
Pxy = PxyTemp; 

% PyyInv = Pyy\eye(q);     % efficient least squares
% K = Pxy*PyyInv; 
K = Pxy/Pyy; 

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

% only for comparison: 
Kv = K*(yMeas - yAggregated);
xPlusvdM = xMinus + Kv; % standard formulation of vdMerwe
disp(['max. Abweichung xPlus (add.):', num2str(max(abs(xPlusvdM - xPlus)))])

% mind that Kolas proved the fully augmented case. For additive noise, the
% measurement update of P-Matrix must be slightly adapted:
diffxPlusFromSigmaX = sigmaX - xPlus; 
PPlusReformulatedKolasFullyAugmented = (Wc.*diffxPlusFromSigmaX)*diffxPlusFromSigmaX'; 
% PPlusGentschTemp = Wc.*sigmaX*sigmaX'; % false implementation in ABC-Code cukf_update1.m:  P = (Xs*W*Xs.');
PPlusReformulatedKolasAdditive = PPlusReformulatedKolasFullyAugmented + Q + K*R*K'; 

% compute the aggregation through a recursive sum: 
PPlusTemp = zeros(nStates); % initialization
for k = 1:nSigmaPoints
    PPlusTemp = PPlusTemp + Wc(k)*((sigmaX(:,k) - xPlus)*(sigmaX(:,k) - xPlus)');
end
PPlusTemp = PPlusTemp + Q + K*R*K'; 
PPlus = 0.5*(PPlusTemp + PPlusTemp'); 

% only for comparison: 
PPlusTempvdM = PMinusMat - K*PyyMat*K'; 
PPlusvdM = 1/2*(PPlusTempvdM + PPlusTempvdM');  % regularization
disp(['max. Abweichung PPlus (add.): ', ...
      num2str(max(max(abs(PPlusvdM - PPlus))))])

% make sure PPlus is symmetric:
PPlusMat = 1/2*(PPlusReformulatedKolasAdditive + PPlusReformulatedKolasAdditive');   
disp(['sum of PPlus diagonal (add.): ', num2str(sum(diag(PPlus)))])

end