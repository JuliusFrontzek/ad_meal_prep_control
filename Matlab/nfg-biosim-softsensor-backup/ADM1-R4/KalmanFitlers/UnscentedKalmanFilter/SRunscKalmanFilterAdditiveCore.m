%% Version
% (R2022b) Update 5
% Erstelldatum: 08.10.2023
% Autor: Simon Hellmann

function [xPlus,SPlus] = SRunscKalmanFilterAdditiveCore(xOld,SOld,...
                                tSpan,feedInfo,yMeas,params,SQ,SR,f,g)
% compute time and measurement update acc. to Joseph-version of the UKF
% acc. to Kolas et al. (2009) with clipping wherever possible
% assume measurements of ADM1-R4-Core

global counterSigmaInit
global counterSigmaProp
global counterX

% xPlus - new state estimate
% SPlus - upper cholesky factor of new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% SOld - upper cholesky factor of old state error covariance matrix 
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% SQ - lower chol. factor of power spectral density matrix of process noise
% SR - lower chol. factor of covariance matrix of measurement noise
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
nSigmaPoints = 2*nStates + 1;

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
% beta = 0; 
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
% kappa = 3 - nStates;  % acc. to Julier & Uhlmann
kappa = 0.0;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
% lambda = 1;     % Vachhani
gamma = sqrt(nStates + lambda); % scaling parameter
gamma = 1;  % XY just to check: that delivers good estimations!

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% Choose Sigma Points

sigmaXInit = [xOld, repmat(xOld,1,nStates) + gamma*SOld, ...
                    repmat(xOld,1,nStates) - gamma*SOld]; 

% % Apply clipping to negative Sigma Points: 
% if any(any(sigmaXInit < 0))
%     sigmaXInit(sigmaXInit < 0) = 0; 
%     counterSigmaInit = counterSigmaInit + 1;
% end

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

% % if any propagated sigma points violate constraints, apply clipping: 
% if any(any(sigmaXProp < 0))
%     sigmaXProp(sigmaXProp < 0) = 0; 
%     counterSigmaProp = counterSigmaProp + 1;
% end

%% Aggregate Sigma Points to Priors for x and P
xMinus = sum(Wx.*sigmaXProp,2);  % state prior

% if any state priors violate constraints, apply clipping:

% if any(any(xMinus < 0))
%     xMinus(xMinus < 0) = 0; 
%     counterX = counterX + 1; 
% end

% perform cholupdate of PMinus (20 and 21 from vdMerwe 2001):
[~,SMinusTemp] = qr([sqrt(Wi)*(sigmaXProp(:,2:end) - xMinus), SQ]','econ');  
SMinus = cholupdate(SMinusTemp, sqrt(Wc0)*(sigmaXProp(:,1) - xMinus)); % upper triangular matrix

%------------------Measurement Update (MU)---------------------------------

% optional: draw new sigma point acc. to a priori estimate xMinus and thus 
% overwrite sigmaXProp: 
% sigmaXProp = [xMinus, repmat(xMinus,1,nStates) + gamma*SMinus, ...
%                       repmat(xMinus,1,nStates) - gamma*SMinus]; 

% or omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them:
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    Y(:,mm) = g(sigmaXProp(:,mm)); 
end
% aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

%% 4. compute Kalman Gain

% perform cholupdate of PMinus (20 and 21 from vdMerwe 2001): 
[~,SyTemp] = qr([sqrt(Wi)*(Y(:,2:end) - yAggregated), SR]','econ');  
Sy = cholupdate(SyTemp, sqrt(Wc0)*(Y(:,1) - yAggregated)); % upper triangular matrix

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
diffXPriorFromSigma = sigmaXProp - xMinus; 
% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

K = Pxy/Sy/Sy';     % Kalman gain acc. to Arrese (2016): S is always the 
% lower triangular cholesky factor in vdMerwe. However, cholupdate in
% matlab only considers upper triangular factors and returns upper
% triangular factors as well, so transposition is necessary.

%% compute posteriors:

xPlus = xMinus + K*(yMeas' - yAggregated); 
U = K*Sy';  % transposition because of same reason as for K

nColsU = size(U,2); % # columns of U
% conduct (29) of vdMerwe (2001):
SPlusTemp = SMinus; % initialization
for k = 1:nColsU
    SPlusTemp = cholupdate(SPlusTemp, U(:,k), '-'); 
end
SPlus = SPlusTemp;

end