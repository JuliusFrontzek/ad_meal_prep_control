%% Version
% (R2022b) Update 6
% Erstelldatum: 06.10.2023
% last modified: 12.10.2023
% Autor: Simon Hellmann

function [xPlus,PPlus,fCount,nIter] = constrUnscKalmanFilterKolasAdditiveCore(xOld,POld, ...
                                    tSpan,feedInfo,yMeas,params,Q,R,f,g)

% compute time and measurement update of constrained UKF acc. to  Kolas et al. 
% (2009), Tab. 10 but for additive noise (abs. coordinates).
% Note: 3 differences between add. noise and fully augm. case:
% 1. augmentation, 
% 2. computation of PMinus, which includes +Q (add. noise case) or doesnt 
% (augmented system noise and fully augmented case)
% 3. computation of posterior of P: add + K*R*K' + Q compared with fully
% augmented case!

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% fCount - # function calls of fmincon
% nIter - # iterations until convergence
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

global counterSigmaInit
global counterSigmaProp
% global counterSigmaX
global counterX
global counterSigmaXcUKF

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
gamma = sqrt(nStates + lambda); % scaling parameter
% gamma = 0.2;  % XY just to check

% weights acc. Diss vdM, (3.12) (Scaled Unscented Transformation): 
Wx0 = lambda/(nStates + lambda); 
Wc0 = lambda/(nStates + lambda) + 1 - alpha^2 + beta; 
Wi = 1/(2*(nStates + lambda)); 
Wx = [Wx0, repmat(Wi,1,nSigmaPoints-1)]; % for state aggregation
Wc = [Wc0, repmat(Wi,1,nSigmaPoints-1)]; % for covariance aggregation

%% Choose Sigma Points
sqrtPOld = schol(POld);  % cholesky factorization acc. to EKF/UKF toolbox from Finland 

sigmaXInit = [xOld, repmat(xOld,1,nStates) + gamma*sqrtPOld, ...
                        repmat(xOld,1,nStates) - gamma*sqrtPOld]; 

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
    XAtBeginOfInt = sigmaXInit;    % Startwerte für erstes Intervall (wird nicht mehr überschrieben)
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

% % if any state priors violate constraints, apply clipping: 
% if any(any(xMinus < 0))
%     xMinus(xMinus < 0) = 0; 
%     counterX = counterX + 1;
% end

% aggregate state error cov. matrix P:
diffXPriorFromSigma = sigmaXProp - xMinus; 
PMinus = Wc.*diffXPriorFromSigma*diffXPriorFromSigma' + Q; % adapted for additive noise case acc. to Kolas, Tab. 5

%% ------------------Measurement Update (MU)------------------------------

% omit to choose new sigma points for measurement update (acc. to Kolas,
% Table 4 or Vachhani 2006)!

%% Derive Sigma-Measurements and aggregate them
Y = nan(q,nSigmaPoints);    % allocate memory
for mm = 1:nSigmaPoints
    Y(:,mm) = g(sigmaXProp(:,mm)); 
end

% aggregate outputs of sigma points in overall output:
yAggregated = sum(Wx.*Y,2);

% consider inequalities acc. to fmincon documentation: allow only positive 
% state values (=concentrations):
A = -eye(nStates); 
b = zeros(nStates,1);  

sigmaXOpt = nan(nStates,nSigmaPoints);    % allocate memory

%%%%%%%%%%%%%%%%%%%%%
%% run constrained optimization to determine sigmaX without gradients/Hess
%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('fmincon',...
'Display','none');
% tic
% optimize all updated sigma points: 
for k = 1:nSigmaPoints

    ukfCostFun = @(sigmaX) evaluateCUKFCostFunCore(sigmaX,sigmaXProp(:,k), ...
                                yMeas',R,PMinus,g); 
    % choose the old sigmaXProp as initial value for optimization:
    sigmaX0 = sigmaXProp(:,k); 
    [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(ukfCostFun,sigmaX0,A,b,[],[],[],[],[],options); 
%         [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(ukfCostFun,sigmaX0,[],[],[],[],lb,ub,[],options); 

end
% output
% toc

%%%%%%%%%%%%%%%%%%%%%
%% run constrained optimization to determine sigmaX with gradients
%%%%%%%%%%%%%%%%%%%%%

% % setUp gradient for fmincon
% options = optimoptions('fmincon',...
% 'SpecifyObjectiveGradient',true,'Display','none');
% % tic
% % optimize all updated sigma points: 
% for k = 1:nSigmaPoints
% 
%     gradCostFun = @(sigmaX) evaluateGradientCUKFCostFunCore(sigmaX,sigmaXProp(:,k), ...
%                                 yMeas',R,PMinus,g); 
%     % choose the old sigmaXProp as initial value for optimization:     
%     sigmaX0 = sigmaXProp(:,k); 
%     [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(gradCostFun,sigmaX0,A,b,[],[],[],[],[],options); 
% %     [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(gradCostFun,sigmaX0,[],[],[],[],lb,ub,[],options); 
% 
% end 
% % toc
% % output

%%%%%%%%%%%%%%%%%%%%%
%% run constrained optimization to determine sigmaX with gradients & Hess
%%%%%%%%%%%%%%%%%%%%%

% % setUp gradient and Hessian for fmincon:
% myHessFcn = @(sigmaX,lambda) evaluateHessCUKFCostFunCore(sigmaX,lambda,R,PMinus); 
% options = optimoptions('fmincon',...
%     "SpecifyObjectiveGradient",true,...
%     'HessianFcn',myHessFcn,'Display','none');
% 
% % tic
% % optimize all updated sigma points: 
% for k = 1:nSigmaPoints
%     gradCostFun = @(sigmaX) evaluateGradientCUKFCostFunCore(sigmaX,sigmaXProp(:,k), ...
%                                 yMeas',R,PMinus,g); 
%     % choose the old sigmaXProp as initial value for optimization:  
%     sigmaX0 = sigmaXProp(:,k); 
%     [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(gradCostFun,sigmaX0,A,b,[],[],[],[],[],options); 
% %     [sigmaXOpt(:,k),fval,exitflag,output] = fmincon(gradCostFun,sigmaX0,[],[],[],[],lb,ub,[],options); 
% 
% end 
% % toc
% % output

% % this clipping should no longer be required thanks to optimization:
% % if updated sigma points violate constraints, apply clipping: 
% if any(any(sigmaXOpt < 0))
%     sigmaXOpt(sigmaXOpt < 0) = 0;
%     counterSigmaXcUKF = counterSigmaXcUKF + 1;
% end

%% compute posteriors:
xPlus = sum(Wx.*sigmaXOpt,2); 

% mind: Kolas considers fully augmented case, so computation of
% posteriors must be adapted. I believe Vachhani (who also considers
% additive noise case) is wrong!
diffxPlusFromSigmaX = sigmaXOpt - xPlus; 

% I think just as for the unconstrained case, we need to adapt computation
% of the posterior of P by incorporating K, Q and R, so calculate K:
% compute cov. matrix of output Pyy:
diffYFromSigmaOutputs = Y - yAggregated; 
Pyy = Wc.*diffYFromSigmaOutputs*diffYFromSigmaOutputs' + R;

% compute cross covariance matrix states/measurements:
Pxy = Wc.*diffXPriorFromSigma*diffYFromSigmaOutputs'; 

% PyyInv = Pyy\eye(q);     % efficient least squares
% K = Pxy*PyyInv;
K = Pxy/Pyy; 

PPlusKolasFullyAugmented = Wc.*diffxPlusFromSigmaX*diffxPlusFromSigmaX'; 
% adapt Kolas (2009) just like for the unconstrained case:
PPlusKolasAdditive = PPlusKolasFullyAugmented + K*R*K' + Q; % actually different formula for additive noise case!
PPlusVachhaniTemp = PPlusKolasFullyAugmented; % Vachhani (2006), (25)

% make sure PPlus is symmetric:
PPlus = 1/2*(PPlusVachhaniTemp + PPlusVachhaniTemp');
% disp(['sum of PPlus diagonal (cUKF-add.): ', num2str(sum(diag(PPlus)))])

%% return # function evaluations and # iterations 
% but only if someone explicitly asks for them when calling this function:
fCount = output.funcCount; 
nIter = output.iterations; 
if nargout < 3
    fCount = []; 
    nIter = [];
end

end