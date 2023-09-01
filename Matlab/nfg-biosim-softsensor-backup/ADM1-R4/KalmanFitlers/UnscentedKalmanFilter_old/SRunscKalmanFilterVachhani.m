function [xPlus,PPlus,K] = SRunscKalmanFilterVachhani(xOld,SOld,tSpan,feedInfo,yMeas,AC,Sr,Sq)
% nach Vachhani et al. (2006): passe die Gewichte für die Sigma-Punkte im
% Time Update an, unterlasse aber die Optimierung im Measurement Update.
% Nutze außerdem die numerisch stabilere SR-Erweiterungen nach v.d.M., 2001

global counterSigmaInit
global counterSigmaProp
global counterX

% time and measurement update acc. to the unscented Kalman Filter
% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% SOld - lower Cholesky factor of POld (the old state error covariance
% matrix)
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% AC - struct with stoichiometric coefficients and aggregated constant
% model parameters
% Sr = lower Cholesky factor of measurement noise cov. matrix R
% Sq = lower Cholesky factor of process noise cov. matrix Q

% make print-statement if xOld contains any negative concentrations: 
if any(xOld<0)
    % print statement: 
%     disp(['neg. Konzentration entdeckt bei ',num2str(tSpan(1)),' Tagen'])
%     xOld(xOld<0)
    % apply clipping: 
    xOld(xOld < 0) = 0; 
    counterX = counterX + 1
end 

nStates = length(xOld); 
% nMeas = length(yMeas); 
nSigmaPoints = 2*nStates + 1; 

%% Time Update (TU)

% define scaling parameters and weights: 
kappa = 0.05;  % default (Vachhani 2006): 1
gamma = sqrt(nStates + kappa); % standard scaling parameter

%% Choose Sigma Points while respecting state constraints
% sqrtPOld = chol(SOld);  % cholesky factorization

% check if state constraints are violated, and if so, modify scaling
% parameters: 
% sij = sqrtPOld;
sij = [SOld, -SOld];    % 2n Einträge

xMin = zeros(nStates,1);        % no negative concentrations
xMax = 1E9*ones(nStates,1);     % there is theoretically no upper limit on concentrations, but set a practical one

theta_ik = nan(1,2*nStates);     % allocate memory (früher: 1,nStates)
% go through sij column by column to recompute the scaling paramters of all 
% sigma points: 
for col = 1:2*nStates % nStates
    currCol = sij(:,col);   % current column
    posScalingIdx = currCol > 0; 
    negScalingIdx = currCol < 0; 
%     negScalingIdx = ~posScalingIdx; % alt
    
    posScale = currCol(posScalingIdx); 
    % count positive scaling parameters into theta_1k: 
    posRescaling = (xMax(posScalingIdx) - xOld(posScalingIdx))./posScale; 
    % hier kann nichts negatives rauskommen, also abs nicht mehr nötig
    theta_1k = min(posRescaling);% alt: min(abs(posRescaling)); 

    negScale = currCol(negScalingIdx); 
    % count all other scaling parameters into theta_2k: 
    negRescaling = (xMin(negScalingIdx) - xOld(negScalingIdx))./negScale; 
    theta_2k = min(negRescaling); % alt: min(abs(negRescaling));
    theta_ik(col) = min([gamma, theta_1k, theta_2k]);
end
% XY: ohne for-Schleife möglich (Gespräch Patrick/Michael)?

sigmaXInit = [xOld, repmat(xOld,1,2*nStates) + theta_ik.*sij]; 

% Apply clipping to negative Sigma Points: 
if any(any(sigmaXInit < 0))
    % apply clipping to sigma points: 
    sigmaXInit(sigmaXInit < 0) = 0; 
    counterSigmaInit = counterSigmaInit + 1; 
end

%% Propagate Sigma Points
sigmaXProp = nan(nStates, 2*nStates + 1); % allocate memory

tEvents = feedInfo(:,1);    % feeding time points (on/off)
idxRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2);
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
    odeFun = @(t,X) BMR4_AB_h2o_ode(t,X,feedVector,AC);
    for k = 1:2*nStates + 1
        [~,XTUSol] = ode15s(odeFun,tEval,sigmaXInit(:,k));
        sigmaXProp(:,k) = XTUSol(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan; tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtEndOfInt = sigmaXInit;    % Startwerte für erstes Intervall (wird überschrieben)
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for m = 1:nIntervals
        feedVector = feedInfo(m,:);
        tEval = [tOverall(m), tOverall(m+1)];
        odeFun = @(t,X) BMR4_AB_h2o_ode(t,X,feedVector,AC);         
        for kk = 1:2*nStates + 1
            [~,XTUSol] = ode15s(odeFun,tEval,XAtEndOfInt(:,kk));
            XAtEndOfInt(:,kk) = XTUSol(end,:)';
        end
    end
    sigmaXProp = XAtEndOfInt;
end

% % check if any propagated sigma points violate constraints: 
% if any(any(sigmaXProp < 0))
%     % apply clipping to sigma points: 
%     sigmaXProp(sigmaXProp < 0) = 0; 
%     counterSigmaProp = counterSigmaProp + 1;
% end

%% compute weights acc. to Vachhani 2006: 
STheta = sum(theta_ik); 
a = (2*kappa - 1)/(2*(nStates + kappa)*(STheta - nSigmaPoints*sqrt(nStates + kappa)));
b = 1/(2*(nStates + kappa)) - (2*kappa - 1)/(2*sqrt(nStates + kappa)*(STheta - nSigmaPoints*sqrt(nStates + kappa)));
WiExCenter = a*theta_ik + b; 
% because numel(Wi) = nSigmaPoints - 1 and sum != 1, we need a zero'th weight: 
Wi0 = 1 - sum(WiExCenter); 
Wi = [Wi0,WiExCenter]; 

%% Aggregate Sigma Points to Priors for x and P
% (do NOT distinguish between weighting of states and covariances)

% aggregate state prior:
xMinus = sum(Wi.*sigmaXProp,2);  

% aggregate sqrt of state error cov. matrix P in 2-step procedure:
% step 1: qr-decomposition:
diffXPriorFromSigmaExCenter = sigmaXProp(:,2:end) - xMinus;     % excluding the central sigma point (the first one)
augMatTU = [sqrt(WiExCenter).*diffXPriorFromSigmaExCenter, Sq]; % liegende Matrix
SMinusTemp = qr(augMatTU','econ');   % qr-decomp (20) in v.d.M. 2001, upper triangle
% step 2: cholupdate:
diffFromCentralSigmaPoint = sigmaXProp(:,nStates+1) - xMinus; 
SMinus = cholupdate(SMinusTemp, sqrt(Wi0)*diffFromCentralSigmaPoint,"+")'; % (21) in v.d.M., 2001, lower triangle
% XY: evtl. um signum-funktion mit dem Vorzeichen von Wi0 erweitern, da das
% auch negativ werden kann

%% Derive Sigma-Measurements and aggregate them:

Y = BMR4_AB_mgl_h2o_mat(sigmaXProp',AC.c)'; 
yAggregated = sum(Wi.*Y,2);

%% Measurement Update (MU):

% omit to choose new sigma points for measurement update (acc. to Vachhani
% 2006)!

% aggregate sqrt of cov. matrix Pyy in 2-step procedure:
% step 1: qr-decomposition:
diffYFromSigmaExCenter = Y(:,2:end) - yAggregated; % excluding the central sigma point (the first one)
augMatMU = [sqrt(WiExCenter).*diffYFromSigmaExCenter, Sr]; % liegende Matrix
SYTemp = qr(augMatMU','econ');   % qr-decomp (24) in v.d.M. 2001, upper triangle
% step 2: cholupdate:
diffFromCentralSigmaOutput = Y(:,nStates+1) - yAggregated; 
Sy = cholupdate(SYTemp, sqrt(Wi0)*diffFromCentralSigmaOutput,"+")'; % lower triangle
% XY: evtl. um signum-funktion mit dem Vorzeichen von Wi0 erweitern

% compute cross covariance matrix states/measurements:
diffXPriorFromSigma = sigmaXProp - xMinus; 
diffYFromSigmaOutputs = Y - yAggregated; 
Pxy = Wi.*diffXPriorFromSigma*diffYFromSigmaOutputs'; % (26) in v.d.M., 2001 

%% 4. Kalman Gain and actual MU:
% PyyInv = Pyy\eye(nMeas); 
K = (Pxy/Sy')/Sy; % Sy must be lower triangle, Sy' upper
Kv = K*(yMeas' - yAggregated);

xPlus = xMinus + Kv;

U = K*Sy; % (28) in v.d.M., 2001
[~,nCols] = size(U); 

% STemp must be upper triangle since cholupdate ignores the lower part of its argument!
STemp = SMinus';      
for col = 1:nCols
    STemp = cholupdate(STemp,U(:,col),"-")'; 
end

SPlus = STemp'; % turn back into lower triangle 
PPlusTemp = SPlus*SPlus'; 

% make sure PPlus remains symmetric
PPlus = 1/2*(PPlusTemp + PPlusTemp'); 

end
