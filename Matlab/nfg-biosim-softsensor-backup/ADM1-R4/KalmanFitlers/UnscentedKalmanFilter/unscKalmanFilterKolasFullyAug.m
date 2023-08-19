function [xPlus,PPlus,K] = unscKalmanFilterKolasAug(xOld,POld,tSpan,feedInfo,yMeas,AC,R,Q)
% nach Kolas et al. (2009) Fully Augmented: alternative Formulierung des 
% Measurement Updates mit und Clipping an allen möglichen Stellen

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX

% time and measurement update acc. to the unscented Kalman Filter
% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% AC - struct with stoichiometric coefficients and aggregated constant
% model parameters

% tSpan(1)

% make print-statement if xOld contains any negative concentrations: 
if any(xOld<0)
%     % print statement: 
%     disp('Clipping at old estimate')
%     disp(num2str(xOld))
    % apply clipping: 
    xOld(xOld < 0) = 0; 
    counterX = counterX + 1; 
end 

nStates = length(xOld); 
nSignals = length(yMeas); 
nSigmaPoints = 2*(2*nStates+nSignals) + 1; 

%% Time Update (TU)

% define scaling parameters and weights: 
alpha = 1;  % Kolas 2009, (18)
beta = 2;   % for Gaussian prior (Diss vdM, S.56)
kappa = 0.001;  % leichte Abweichung zu Kolas (er nimmt 0)
lambda = alpha^2*(nStates + kappa) - nStates; 
gamma = sqrt(nStates + lambda); % scaling parameter

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% Choose Sigma Points
POldAug = blkdiag(POld,Q,R); 
sqrtPOld = schol(POldAug);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

xOldAug = [xOld;zeros(nStates,1);zeros(nSignals,1)]; 

sigmaXInit = [xOldAug, repmat(xOldAug,1,2*nStates+nSignals) + gamma*sqrtPOld, ...
                       repmat(xOldAug,1,2*nStates+nSignals) - gamma*sqrtPOld]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInit < 0))
%     disp('Clipping at initial sigma points:')
%     disp(num2str(sigmaXInit))
    % apply clipping to sigma points: 
    sigmaXInit(sigmaXInit < 0) = 0; 
    counterSigmaInit = counterSigmaInit + 1; 
end

%% Propagate Sigma Points
sigmaXInitx = sigmaXInit(1:nStates,:); % extract only entries belonging to x
% note: we only want first n entries of the the propagated sigma points
% which belong to x (not those of xv and xw)!
sigmaXProp = nan(nStates, nSigmaPoints); % allocate memory

tEvents = feedInfo(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
%     odeFun = @(t,X) f(t,X,feedVector,AC);
    odeFun = @(t,X) BMR4_AB_h2o_ode(t,X,feedVector,AC);
    for k = 1:nSigmaPoints
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInitx(:,k));
        sigmaXProp(:,k) = XTUSol(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan; tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtEndOfInt = sigmaXInitx;    % Startwerte für erstes Intervall (wird überschrieben)
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVector = feedInfo(m,:);
        tEval = [tOverall(m), tOverall(m+1)];
%         odeFun = @(t,X) f(t,X,feedVector,AC);
        odeFun = @(t,X) BMR4_AB_h2o_ode(t,X,feedVector,AC);
        for kk = 1:nSigmaPoints
            [~,XTUSol] = ode15s(odeFun,tEval,XAtEndOfInt(:,kk));
            XAtEndOfInt(:,kk) = XTUSol(end,:)';
        end
    end
    sigmaXProp = XAtEndOfInt;
end

% check if any propagated sigma points violate constraints: 
if any(any(sigmaXProp < 0))
%     disp('Clipping at propagated sigma points:')
%     disp(num2str(sigmaXProp))
    % apply clipping to sigma points: 
    sigmaXProp(sigmaXProp < 0) = 0; 
    counterSigmaProp = counterSigmaProp + 1; 
end

%% Aggregate Sigma Points to Priors for x and P

% aggregate state prior:
xMinus = sum(Wx.*sigmaXProp,2);  

% check if any state priors violate constraints: 
if any(any(xMinus < 0))
%     disp('Clipping at prior:')
%     disp(num2str(xMinus))
    % apply clipping to state priors: 
    xMinus(xMinus < 0) = 0; 
    counterX = counterX + 1; 
end

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma'; 

%------------------Measurement Update (MU)---------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them:

Y = BMR4_AB_mgl_h2o_mat(sigmaXProp',AC.c)'; 
yAggregated = sum(Wx.*Y,2);

%% 4. compute Kalman Gain

% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs';

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

PyyInv = Pyy\eye(nSignals);     % efficient least squares
K = Pxy*PyyInv; 

Kv = K*(yMeas' - yAggregated);

xPlus = xMinus + Kv;
PPlusTemp = PMinus - K*Pyy*K'; 

PPlus = 1/2*(PPlusTemp + PPlusTemp');   % make sure PPlus is symmetric!
sum(diag(PPlus)) % show divergence live

end