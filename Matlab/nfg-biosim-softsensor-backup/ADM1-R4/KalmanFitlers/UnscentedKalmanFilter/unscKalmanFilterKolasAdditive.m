%% Version
% (R2022b) Update 5
% Erstelldatum: 30.08.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = unscKalmanFilterKolasAdditive(xOld,POld,...
                                tSpan,feedInfo,yMeas,params,Q,R,f,g)
% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) with clipping wherever possible

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
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
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

% if xOld contains negative concentrations, apply clipping: 
if any(xOld<0)
    xOld(xOld < 0) = 0; 
    counterX = counterX + 1;
end 

nStates = numel(xOld); 
q = numel(yMeas); 
nSigmaPoints = 2*nStates + 1;

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 0.5;  % Kolas 2009, (18)
% beta = 0; 
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
% kappa = 3 - nStates;  % acc. to Julier & Uhlmann
kappa = 0.05;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
% lambda = 1;     % Vachhani
gamma = sqrt(nStates + lambda); % scaling parameter
% gamma = 0.5;  % XY just to check: that delivers good estimations!

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% Choose Sigma Points
sqrtPOld = schol(POld);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 
% sqrtPOld = chol(POld);  % cholesky factorization acc. to Matlab

sigmaXInit = [xOld, repmat(xOld,1,nStates) + gamma*sqrtPOld, ...
                    repmat(xOld,1,nStates) - gamma*sqrtPOld]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInit < 0))
    sigmaXInit(sigmaXInit < 0) = 0; 
    counterSigmaInit = counterSigmaInit + 1;
end

%% Propagate Sigma Points
sigmaXProp = nan(nStates, nSigmaPoints); % allocate memory

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
    for k = 1:nSigmaPoints
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInit(:,k));
        sigmaXProp(:,k) = XTUSol(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtBeginOfInt = sigmaXInit;    % Startwerte für erstes Intervall (wird nicht überschrieben)
    XAtEndOfInt = nan(size(sigmaXInit)); % allocate memory
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVolFlow = feedInfo(m,2);
        xInCurr = feedInfo(m,3:end)';   % current inlet concentrations
        tEval = [tOverall(m), tOverall(m+1)];
        odeFun = @(t,X) f(X,feedVolFlow,xInCurr,th,c,a); 
        for kk = 1:nSigmaPoints
            [~,XTUSol] = ode15s(odeFun,tEval,XAtBeginOfInt(:,kk));
            XAtEndOfInt(:,kk) = XTUSol(end,:)';
        end
        XAtBeginOfInt = XAtEndOfInt; % overwrite for next interval
    end
    sigmaXProp = XAtEndOfInt;
end

% if any propagated sigma points violate constraints, apply clipping: 
if any(any(sigmaXProp < 0))
    sigmaXProp(sigmaXProp < 0) = 0; 
    counterSigmaProp = counterSigmaProp + 1;
end

%% Aggregate Sigma Points to Priors for x and P
xMinus = sum(Wx.*sigmaXProp,2);  % state prior

% if any state priors violate constraints, apply clipping:

if any(any(xMinus < 0))
    xMinus(xMinus < 0) = 0; 
    counterX = counterX + 1; 
end

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma' + Q; 

%------------------Measurement Update (MU)---------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    Y(:,mm) = g(sigmaXProp(:,mm),c); 
end
% aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

%% 4. compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs' + R;

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% PyyInv = Pyy\eye(q);     % efficient least squares
% K = Pxy*PyyInv; 
K = Pxy/Pyy; 

%% update propagated sigma points individually 
% use alternative formulation of Kolas 2009, eq. (23):
sigmaX = sigmaXProp + K*(repmat(yMeas',1,nSigmaPoints) - Y);

% if updated sigma points violate constraints, apply clipping: 
if any(any(sigmaX < 0))
    sigmaX(sigmaX < 0) = 0;
    counterSigmaX = counterSigmaX + 1;
end

%% compute posteriors:
xPlus = sum(Wx.*sigmaX,2); 

% only for comparison: 
Kv = K*(yMeas' - yAggregated);
xPlusvdM = xMinus + Kv; % standard formulation of vdMerwe
disp(['max. Abweichung xPlus (add.):', num2str(max(abs(xPlusvdM - xPlus)))])

diffxPlusFromSigmaX = sigmaX - xPlus; 
PPlusKolasFullyAugmentedCase = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX'; 
PPlusTemp = PPlusKolasFullyAugmentedCase + Q + K*R*K'; % actually different formula for additive noise case!

% only for comparison: 
PPlusTempvdM = PMinus - K*Pyy*K'; 
disp(['max. Abweichung PPlus (add.): ', ...
      num2str(max(max(abs(PPlusTempvdM - PPlusTemp))))])

PPlus = 1/2*(PPlusTemp + PPlusTemp');   % make sure PPlus is symmetric!
disp(['sum of PPlus diagonal (add.): ', num2str(sum(diag(PPlus)))])

end