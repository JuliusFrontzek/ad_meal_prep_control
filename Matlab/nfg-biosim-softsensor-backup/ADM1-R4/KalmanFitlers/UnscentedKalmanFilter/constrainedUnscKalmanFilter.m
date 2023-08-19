function [xPlus,PPlus] = constrainedUnscKalmanFilter(xOld,POld,tSpan,feedInfo,yMeas,AC,R,Q)
% nach Vachhani et al. (2006) (full) mit Berücksichtigung von constraints

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

% print first entry of tSpan so we know where we are if an error occurs: 
tSpan(1)

nStates = length(xOld); 
% nMeas = length(yMeas); 
nSigmaPoints = 2*nStates + 1; 

%% Time Update (TU) in 4 steps

% define scaling parameters (Diss v.d.Merwe):
alpha = 0.01;  % [0,1]. should be "small" 
beta = 2;     % optimal for Gaussian distribution
kappavdM = 3; % >= 0  for pos. definiteness of P (v.d.M)
lambda = alpha^2*(nStates + kappavdM) - nStates; % Diss v.d.M., eq. (3.11)
% gamma = sqrt(nStates + lambda); % scaling parameter

% kappa = 3 - nStates;  % recommendation Julier 1995
% lambda = kappa; % v.d.M. calls it lambda, Vachhani kappa
kappa = 1; % recommendation Vachhani 2006 to avoid neg. definiteness of cov. matrices
gamma = sqrt(nStates + kappa); % default scaling parameter

%% Choose Sigma Points
try 
    sqrtPOld = chol(POld);  % cholesky factorization
catch 
    % make sure POld is symmetric pos. definite: 
    eps = 1e-5;     % all negative eigenvalues of POld are increased up to eps
    POldspd = makePSymmPosDefWhile(POld, eps); 
    sqrtPOld = chol(POldspd);  % cholesky factorization
end

% sqrtPOld = chol(POld);  % cholesky factorization

% check if state constraints are violated, and if so, modify scaling
% parameters: 
sij = [sqrtPOld, -sqrtPOld];    % 2n Einträge

xMin = zeros(nStates,1);       % no negative concentrations
xMax = 1E9*ones(nStates,1);    % there is theoretically no upper limit on concentrations, but set a practical one

theta_ik = nan(1,2*nStates);     % allocate memory (früher: 1,nStates)
% go through sij column by column to recompute the scaling paramters of all 
% sigma points: 
for col = 1:2*nStates % nStates
    currCol = sij(:,col);   % current column
    posScalingIdx = currCol > 0; 
    negScalingIdx = ~posScalingIdx; 
    
    posScale = currCol(posScalingIdx); 
    % count positive scaling parameters into theta_1k: 
    posRescaling = (xMax(posScalingIdx) - xOld(posScalingIdx))./posScale; 
    theta_1k = min(abs(posRescaling)); 

    negScale = currCol(negScalingIdx); 
    % count all other scaling parameters into theta_2k: 
    negRescaling = (xMin(negScalingIdx) - xOld(negScalingIdx))./negScale; 
    theta_2k = min(abs(negRescaling));
%     [theta_1k, theta_2k]
    theta_ik(col) = min([gamma, theta_1k, theta_2k]);
end
% XY: ohne for-Schleife möglich?

sigmaXOld = [xOld, repmat(xOld,1,2*nStates) + theta_ik.*sij]; 

% define weights acc. to v.d.M.: 
% wm0 = kappa/(nStates + kappa); 
% wc0 = lambda/(nStates + lambda) + 1-alpha^2+beta; 
% wc0 = wm0;  % Vachhani 2006
% wmi = 1/(2*(nStates + lambda)); 
% wci = wmi; 
% combine weights into vectors: 
% wm = [wm0, repmat(wmi,1,2*nStates)]; % for state aggregation
% wc = [wc0, repmat(wci,1,2*nStates)]; % for covariance aggregation

% theta_ik only has nStates entries. But we need scaling parameters for all sigma points, 
% that means the central one and all 2*nStates others: 
theta_ik_eff = [0,theta_ik,theta_ik]; 
% theta_ik_eff_P = [0,theta_ik,theta_ik]; 

% Vachhani 2006, eq. (37) corrected by zero'th summand and with distinction
% between weights for states and covariances:
STheta = sum(theta_ik_eff); 
% SThetaP = sum(theta_ik_eff_P);

% define weights acc. to Vachhani 2006:
% kappa = lambda; % kappa (from Vachhani/Julier) is theta (from v.d.Merwe)
a = (2*kappa - 1)/(2*(nStates + kappa)*(STheta - nSigmaPoints*(sqrt(nStates + kappa))));
b = 1/(2*(nStates + kappa)) - (2*kappa - 1)/(2*sqrt(nStates + kappa)*(STheta - nSigmaPoints*(sqrt(nStates + kappa))));
% aP = (2*kappa - 1)/(2*(nStates + kappa)*(SThetaP - nSigmaPoints*(sqrt(nStates + kappa))));
% bP = 1/(2*(nStates + kappa)) - (2*kappa - 1)/(2*sqrt(nStates + kappa)*(SThetaP - nSigmaPoints*(sqrt(nStates + kappa))));

% do NOT distinguish between weighting of states and covariances: 
WiTemp = a*theta_ik + b; 
% because numel(Wi) = nSigmaPoints - 1 and sum != 1, we need a zero'th weight: 
Wi0 = 1 - sum(WiTemp); 
Wi = [Wi0,WiTemp]; 

%% Propagate Sigma Points
sigmaXProp = nan(nStates, nSigmaPoints); % allocate memory

tEvents = feedInfo(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
    odeFun = @(t,X) f(t,X,feedVector,AC);
    for k = 1:nSigmaPoints
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXOld(:,k));
        sigmaXProp(:,k) = XTUSol(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan; tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtEndOfInt = sigmaXOld;    % Startwerte für erstes Intervall (wird überschrieben)
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVector = feedInfo(m,:);
        tEval = [tOverall(m), tOverall(m+1)];
        odeFun = @(t,X) f(t,X,feedVector,AC);
        for kk = 1:nSigmaPoints
            [~,XTUSol] = ode15s(odeFun,tEval,XAtEndOfInt(:,kk));
            XAtEndOfInt(:,kk) = XTUSol(end,:)';
        end
    end
    sigmaXProp = XAtEndOfInt;
end

% print if any values in sigmaXProp are still NaNs: 
anyNaNsLeftInSigmaXProp = any(any(isnan(sigmaXProp)));

%% 4. Aggregate Sigma Points to Priors for x and P

% aggregate state prior:
xMinus = sum(Wi.*sigmaXProp,2);   

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wi.*diffXPriorFromSigma*diffXPriorFromSigma' + Q; 

%% Measurement Update (MU):

% omit to choose new sigma points for measurement update (acc. to Vachhani
% 2006)!

% for every sigma-point, ensure that state constraints are satisfies:
% Neu: Berücksichtige für das measurement update, dass die Zustände positiv
% bleiben müssen. Löse also ein beschränktes Optimierungsproblem: 
A = -eye(nStates); 
b = zeros(nStates,1);  
sigmaXOpt = nan(nStates,nSigmaPoints);    % allocate memory
for k = 1:nSigmaPoints
    costFun = @(x) evaluateCUKFCostFun(x,yMeas,R,sigmaXProp(:,k),PMinus,AC.c); 
    sigmaXOpt(:,k) = fmincon(costFun,xMinus,A,b); 
end 


%% Aggregate Sigma Points to Priors for x and P

% aggregate state prior (Vachhani 2006, eq. (24)):
xPlus = sum(Wi.*sigmaXOpt,2);   

% aggregate state error cov. matrix P (Vachhani 2006, eq. (25)):
diffXPlusFromSigma = sigmaXOpt - xPlus; 
PPlus = Wi.*diffXPlusFromSigma*diffXPlusFromSigma';

end
