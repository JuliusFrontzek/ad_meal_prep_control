%% Version
% (R2022b) Update 6
% Erstelldatum: 06.10.2023
% Autor: Simon Hellmann

%% Unscented und Extended Kalman Filter fürs ADM1-R4-Core

addpath('generatedOutput/')
addpath('modelEquations/')

close all
clear all
clc

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX
global counterSigmaXcUKF
counterSigmaInit = 0;
counterSigmaProp = 0;
counterSigmaX = 0; 
counterX = 0; 
counterSigmaXcUKF = 0;

% Load Measurement Data:
load Messung_ADM1_R4_Core

%% Initialization and Tuning of UKF 
% initial state values:
x0Init = MESS.x0;       % intitial value from createMess-file 
x0 = x0Init;            % x0 will be replaced in every iteration later
% false initial estimate (add random number but ensure positive concentrations)
nStates = length(x0); 
rng('default');     % fix seed for random number generation (for replicable results)
xHat = x0Init.*abs(randn(nStates,1)); 

tMeas = MESS.t;
nSamples = length(tMeas);  % number of measurements taken
% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMeas); 
dt = diff_t(1);             % sampling time
t = [tMeas;tMeas(end)+dt];  % add one time interval at end

% set up raw matrices. all matrices which require initialization have 1 
% column more than the number of measurement instances:
MEAS = MESS.yMeas;
STATES = MESS.x; 
ESTIMATESUKFAdd = zeros(nSamples + 1,nStates);
ESTIMATESUKF_sysID = zeros(nSamples + 1,nStates);
ESTIMATESUKFAug = zeros(nSamples + 1,nStates);
ESTIMATESUKFFullyAug = zeros(nSamples + 1,nStates);
ESTIMATESSRUKF = zeros(nSamples + 1,nStates);
ESTIMATEScUKFNLP = zeros(nSamples + 1,nStates);
ESTIMATEScUKFQP = zeros(nSamples + 1,nStates);
ESTIMATESCKF = zeros(nSamples + 1,nStates);
ESTIMATESEKF = zeros(nSamples + 1,nStates);
COVARIANCEUKFAdd = zeros(nStates,nStates,nSamples + 1);
COVARIANCEUKF_sysID = zeros(nStates,nStates,nSamples + 1);
COVARIANCEUKFAug = zeros(nStates,nStates,nSamples + 1);
COVARIANCEUKFFullyAug = zeros(nStates,nStates,nSamples + 1);
COVARIANCEcUKFNLP = zeros(nStates,nStates,nSamples + 1);
COVARIANCEcUKFQP = zeros(nStates,nStates,nSamples + 1);
COVARIANCECKF = zeros(nStates,nStates,nSamples + 1);
COVARIANCEEKF = zeros(nStates,nStates,nSamples + 1);

buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
% Tune Kalman Filter: 
% measurement uncertainty: 
% R(1,1) = 1E-3*R(1,1);   % gasVolFlow
% R(4,4) = 1E2*R(4,4);    % SIN
% process uncertainty: 
% Q = diag([0.016, 0.555, 1, 1.263, 2.654, 0.972, 2.894, 0.374, 0.948]);
Q = eye(nStates); 

% obtain feeding information:
inputMat = MESS.inputMat;   % [tEvents,feed vol flow,inlet concentrations]
tEvents = inputMat(:,1);    % times of feeding events (on/off)

% if inputMat has first feeding entry only at t > t0, add zero row at t0:
if tEvents(1) > 0
    nColcInputMat = size(inputMat,2);   % number of columns in inputMat
    inputMat = [zeros(nColcInputMat,1),inputMat]; 
end

%% include small plant-model mismatch

trueParams = params.th; 
falseParams = trueParams*(1 + rand);
params.th = falseParams;

%% derive all function handles for EKF from symbolic model equations
% define symbolic variables ("S") (all vectors are column vectors):
xS = sym('x', [nStates,1]);     % states as col. vector
syms uS real                    % input
xiS = sym('xi', [nStates,1]);   % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));    % time-variant parameters (theta)
cS = sym('c', size(params.c));          % time-invariant parameters (known & constant)
aS = sym('a', size(params.a)); % petersen matrix with stoichiometric constants

% obtain symbolic objects of model equations:
dynamics = ADM1_R4_Core_ode_sym(xS, uS, xiS, thS, cS, aS); 
outputs = ADM1_R4_Core_mgl_sym(xS);

% derive function handles of ODEs and output equation: 
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS}); 

% partial derivatives for EKF (symbolic): 
dfdxSym = jacobian(dynamics,xS); 
dhdxSym = jacobian(outputs,xS); 

% convert to numeric function handles: 
dfdx = matlabFunction(dfdxSym, 'Vars', {xS, uS, thS, cS, aS}); 
dhdx = matlabFunction(dhdxSym, 'Vars', {xS}); 

%% create CKF object and tune it
% Set up Cubature Kalman Filter
% ckf = trackingCKF(@stateTransitionFcnADM1R4,@measurementFcnADM1R4_Core,xMinusCKF, ...
%     'HasAdditiveProcessNoise',true, 'HasAdditiveMeasurementNoise',true, ...
%     'MeasurementNoise',R, 'ProcessNoise',Q); % , 'EnableSmoothing',true

%% Tune SR-UKF van der Merwe (2001)

SQ = chol(Q,'lower'); 
SR = chol(R,'lower'); 

%% integrate across all (online) measurement intervals (like in reality):
nRuns = 20;         % # runs
runTime = nan(1,nRuns); % save run times per run here
EstimatesForRMSE1 = nan(nSamples,nStates,nRuns); % analyzed algorithm 1
EstimatesForRMSE2 = nan(nSamples,nStates,nRuns); % analyzed algorithm 2
EstimatesForRMSE3 = nan(nSamples,nStates,nRuns); % analyzed algorithm 3

for runK = 1:nRuns

% Initialize Kalman Filters:
xMinusUKFAdd = xHat;      % to be overwritten
% xMinus = x0Init;    % XY Rania
xMinusUKF_sysID = xMinusUKFAdd;
xMinusUKFAug = xMinusUKFAdd;
xMinusUKFFullyAug = xMinusUKFAdd;
xMinusSRUKF = xMinusUKFAdd;
xMinuscUKFNLP = xMinusUKFAdd; 
xMinuscUKFQP = xMinusUKFAdd; 
% xMinusCKF = xMinusUKF; 
xMinusEKF = xMinusUKFAdd; 

ESTIMATESUKFAdd(1,:) = xMinusUKFAdd;
ESTIMATESUKF_sysID(1,:) = xMinusUKF_sysID;
ESTIMATESUKFAug(1,:) = xMinusUKFAug;
ESTIMATESUKFFullyAug(1,:) = xMinusUKFAug;
ESTIMATESSRUKF(1,:) = xMinusSRUKF;
ESTIMATEScUKFNLP(1,:) = xMinuscUKFNLP;
ESTIMATEScUKFQP(1,:) = xMinuscUKFQP;
% ESTIMATESCKF(1,:) = xMinusCKF;
ESTIMATESEKF(1,:) = xMinusEKF;

P0 = diag((xHat-x0).^2);    % Schneider und Georgakis 
PMinusUKFAdd = P0;      % to overwrite
PMinusUKF_sysID = PMinusUKFAdd; 
PMinusUKFAug = PMinusUKFAdd;
PMinusUKFFullyAug = PMinusUKFAdd;
SMinusSRUKF = chol(PMinusUKFAdd,'upper'); 
PMinuscUKFNLP = PMinusUKFAdd; 
PMinuscUKFQP = PMinusUKFAdd; 
% PMinusCKF = PMinusUKF; 
PMinusEKF = PMinusUKFAdd; 

COVARIANCEUKFAdd(:,:,1) = PMinusUKFAdd; 
COVARIANCEUKF_sysID(:,:,1) = PMinusUKFAdd; 
COVARIANCEUKFAug(:,:,1) = PMinusUKFAug;
COVARIANCEUKFFullyAug(:,:,1) = PMinusUKFFullyAug;
COVARIANCEcUKFNLP(:,:,1) = PMinuscUKFNLP; 
COVARIANCEcUKFQP(:,:,1) = PMinuscUKFQP; 
% COVARIANCECKF(:,:,1) = PMinusCKF; 
COVARIANCEEKF(:,:,1) = PMinusEKF; 

% create UKF object and tune it
% Set up Unscented Kalman Filter
ukf = unscentedKalmanFilter(@stateTransitionFcnADM1R4,@measurementFcnADM1R4_Core,xMinusUKFAdd, ...
    'HasAdditiveProcessNoise',true, 'HasAdditiveMeasurementNoise',true, ...
    'MeasurementNoise',R, 'ProcessNoise',Q); % , 'EnableSmoothing',true

tic
for k = 1:nSamples 
    
    tSpan = [t(k) t(k+1)]; % measurement interval. In reality, this is the 
    % time distance between old and new measurements
        
    yMeas = MESS.yMeas(k,:);    % simulated measurement

    %% get feeding information:
    % pass only relevant feedings during the measurement interval, because
    % only those would be known in reality:
    idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
    tRelEvents = tEvents(idxRelEvents); % Auswertung anhand Index
    
    % find the critcal last feeding event before current measurement interval:
    idxLastEvent = find(tEvents < tSpan(1),1,'last');

    % Case a: constant feeding during measurement interval:
    if isempty(tRelEvents) 
        feedInfo = inputMat(idxLastEvent,:);
    % Case b: feeding changes during measurement interval:
    else
        % use all relevant feeding events during measurement interval and 
        % the critical last one before:
        feedInfo = inputMat([idxLastEvent,idxRelEvents],:);
    end
    %% execute KFs:
    % EKF: 
    [xPlusEKF,PPlusEKF,~] = extendedKalmanFilterCore(xMinusEKF,PMinusEKF, ...
                    tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx);

    % ---UKFs------------------------------
    [xPlusUKF_sysID,PPlusUKF_sysID] = my_UKF_ADM1_Core(ukf,feedInfo,yMeas,tSpan,params,f,g);    
%     [xPlusUKFAdd,PPlusUKFAdd] = unscKalmanFilterKolasAdditiveCore(xMinusUKFAdd,PMinusUKFAdd,...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);
    [xPlusUKFAug,PPlusUKFAug] = unscKalmanFilterKolasAugmentedCore(xMinusUKFAug,PMinusUKFAug,...
                            tSpan,feedInfo,yMeas,params,Q,R,f,g);
    [xPlusUKFFullyAug,PPlusUKFFullyAug] = unscKalmanFilterKolasFullyAugmentedCore(xMinusUKFFullyAug,PMinusUKFFullyAug,...
                            tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPlusSRUKF,SPlusSRUKF] = SRunscKalmanFilterAdditiveCore(xMinusSRUKF,SMinusSRUKF,...
%                             tSpan,feedInfo,yMeas,params,SQ,SR,f,g);
%     [xPluscUKFNLP,PPluscUKFNLP] = constrUnscKalmanFilterKolasAdditiveCore(xMinuscUKFNLP,PMinuscUKFNLP, ...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPluscUKFQP,PPluscUKFQP] = constrUnscKalmanFilterKolasQPAdditiveCore(xMinuscUKFQP,PMinuscUKFQP, ...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPluscUKFQP,PPluscUKFQP] = constrUnscKalmanFilterKolasQPFullyAugmentedCore(xMinuscUKFQP,PMinuscUKFQP, ...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);

%     % ---CKF-----------------------------
%     [xPlusCKF,PPlusCKF] = my_CKF_ADM1_Core(ckf,feedInfo,yMeas,tSpan,params,f,g);
    
    % save results:
    EstimatesForRMSE1(k,:,runK) = xPlusUKF_sysID; % to compute RMSE later
    EstimatesForRMSE2(k,:,runK) = xPlusUKFAug; 
    EstimatesForRMSE3(k,:,runK) = xPlusUKFFullyAug; 

    ESTIMATESEKF(k+1,:) = xPlusEKF';
    ESTIMATESUKF_sysID(k+1,:) = xPlusUKF_sysID';
%     ESTIMATESUKFAdd(k+1,:) = xPlusUKFAdd';
%     ESTIMATESSRUKF(k+1,:) = xPlusSRUKF; 
    ESTIMATESUKFAug(k+1,:) = xPlusUKFAug';
    ESTIMATESUKFFullyAug(k+1,:) = xPlusUKFFullyAug';
%     ESTIMATEScUKFNLP(k+1,:) = xPluscUKFNLP';
%     ESTIMATEScUKFQP(k+1,:) = xPluscUKFQP';
%     ESTIMATESCKF(k+1,:) = xPlusCKF';
    
    COVARIANCEEKF(:,:,k+1) = PPlusEKF; 
    COVARIANCEUKF_sysID(:,:,k+1) = PPlusUKF_sysID; 
%     COVARIANCEUKFAdd(:,:,k+1) = PPlusUKFAdd; 
    COVARIANCEUKFAug(:,:,k+1) = PPlusUKFAug; 
    COVARIANCEUKFFullyAug(:,:,k+1) = PPlusUKFFullyAug; 
%     COVARIANCEcUKFNLP(:,:,k+1) = PPluscUKFNLP;
%     COVARIANCEcUKFQP(:,:,k+1) = PPluscUKFQP;
%     COVARIANCECKF(:,:,k+1) = PPlusCKF; 

    % Update for next iteration...  
    % ... estimated state from Kalman Filter:
    xMinusEKF = xPlusEKF;
    xMinusUKF_sysID = xPlusUKF_sysID; 
%     xMinusUKFAdd = xPlusUKFAdd;
%     xMinusSRUKF = xPlusSRUKF;
    xMinusUKFAug = xPlusUKFAug;
    xMinusUKFFullyAug = xPlusUKFFullyAug; 
%     xMinuscUKFNLP = xPluscUKFNLP;
%     xMinuscUKFQP = xPluscUKFQP;
%     xMinusCKF = xPlusCKF; 

    % ... state error covariance matrices:
    PMinusEKF = PPlusEKF;
    PMinusUKF_sysID = PPlusUKF_sysID; 
%     PMinusUKFAdd = PPlusUKFAdd;
%     SMinusSRUKF = SPlusSRUKF; 
    PMinusUKFAug = PPlusUKFAug;
    PMinusUKFFullyAug = PPlusUKFFullyAug;
%     PMinuscUKFNLP = PPluscUKFNLP;
%     PMinuscUKFQP = PPluscUKFQP;
%     PMinusCKF = PPlusCKF; 

end
runTime(runK) = toc
end

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
EKFOutput = nan(nSamples + 1,q);            % allocate memory
UKFAddOutput = nan(size(EKFOutput));        % allocate memory
UKFOutput_sysID = nan(size(EKFOutput));     % allocate memory
UKFAugOutput = nan(size(EKFOutput));        % allocate memory
UKFFullyAugOutput = nan(size(EKFOutput));   % allocate memory
SRUKFOutput = nan(size(EKFOutput));         % allocate memory
cUKFNLPOutput = nan(size(EKFOutput));       % allocate memory
cUKFQPOutput = nan(size(EKFOutput));        % allocate memory
% CKFOutput = nan(size(EKFOutput));       % allocate memory
for k = 1:nSamples + 1
    EKFOutput(k,:) = g(ESTIMATESEKF(k,:)');
    UKFAddOutput(k,:) = g(ESTIMATESUKFAdd(k,:)');
    UKFOutput_sysID(k,:) = g(ESTIMATESUKF_sysID(k,:)');
    UKFAugOutput(k,:) = g(ESTIMATESUKFAug(k,:)');
    UKFFullyAugOutput(k,:) = g(ESTIMATESUKFFullyAug(k,:)');
    SRUKFOutput(k,:) = g(ESTIMATESSRUKF(k,:)');
    cUKFNLPOutput(k,:) = g(ESTIMATEScUKFNLP(k,:)');
    cUKFQPOutput(k,:) = g(ESTIMATEScUKFQP(k,:)');
%     CKFOutput(k,:) = g(ESTIMATESCKF(k,:)');
end
yClean = MESS.yClean; 
feedVolFlow = inputMat(:,2);    % [l/d]

%% compute goodness of fit for all measurements
RMSSE_UKF_SRUKF = zeros(q,1); % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurementsUKFMatlab = UKFOutput_sysID(2:end,kk); 
    estimatedMeasurementsSRUKF = SRUKFOutput(2:end,kk); 
    
    % get RMSSE between UKF and CKF: 
    RMSSE_UKF_SRUKF(kk) = computeRMSSE(estimatedMeasurementsUKFMatlab,estimatedMeasurementsSRUKF);     
end

%% compute normalized RMSE for all states and all runs: 
nRMSE1 = nan(nStates,nRuns);  % allocate memory
nRMSE2 = nan(nStates,nRuns);
nRMSE3 = nan(nStates,nRuns);
trueStates = MESS.x; 

for runK = 1:nRuns
    for m = 1:nStates
        nRMSE1(m,runK) = compute_normRMSE(trueStates(:,m),EstimatesForRMSE1(:,m,runK)); 
        nRMSE2(m,runK) = compute_normRMSE(trueStates(:,m),EstimatesForRMSE2(:,m,runK)); 
        nRMSE3(m,runK) = compute_normRMSE(trueStates(:,m),EstimatesForRMSE3(:,m,runK)); 
    end
end

% compute mean of nRMSE for measurable and non-measurable states for all algorithms: 
mean_nRMSE1y = mean(nRMSE1([1,2,6],:),2);   % algorithm 1
mean_nRMSE1x = mean(nRMSE1([3,4,5],:),2); 
mean_nRMSE2y = mean(nRMSE2([1,2,6],:),2);   % algorithm 2
mean_nRMSE2x = mean(nRMSE2([3,4,5],:),2); 
mean_nRMSE3y = mean(nRMSE3([1,2,6],:),2);   % algorithm 3
mean_nRMSE3x = mean(nRMSE3([3,4,5],:),2); 

% aggregate all means to one number each: 
meanOfMean_nRMSE1y = mean(mean_nRMSE1y);    % algorithm 1
meanOfMean_nRMSE1x = mean(mean_nRMSE1x); 
meanOfMean_nRMSE2y = mean(mean_nRMSE2y);    % algorithm 2
meanOfMean_nRMSE2x = mean(mean_nRMSE2x); 
meanOfMean_nRMSE3y = mean(mean_nRMSE3y);    % algorithm 3
meanOfMean_nRMSE3x = mean(mean_nRMSE3x); 

%% Plot results

% plot model output based on EKF estimation and compare with real
% measurements:
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figure

% S_ch4: 
subplot(3,1,1)
scatter(tMeas, MESS.yMeas(:,1),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on
plot(tMeas,yClean(:,1),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
% plot(t,EKFOutput(:,1),'DisplayName','EKF-Output',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.8); 
plot(t,UKFOutput_sysID(:,1),'DisplayName','UKF-sysID-Output',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2); 
plot(t,UKFAddOutput(:,1),'DisplayName','UKF-Add-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',2); 
plot(t,SRUKFOutput(:,1),'DisplayName','SR-UKF-Output',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1); 
% plot(t,UKFAugOutput(:,1),'DisplayName','UKF-Aug-Output',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% plot(t,UKFFullyAugOutput(:,1),'DisplayName','UKF-Fully-Aug-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1); 
% plot(t,cUKFNLPOutput(:,1),'DisplayName','cUKF-NLP-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1);
% plot(t,cUKFQPOutput(:,1),'DisplayName','cUKF-QP-Output',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,CKFOutput(:,1),'DisplayName','CKF-Output',...
%      'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
% ylim([0.4,0.85])
ylabel('S_{ch4} [kg/m3]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% S_co2:
subplot(3,1,2)
scatter(tMeas, MESS.yMeas(:,2),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,2),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% plot(t,EKFOutput(:,2),'DisplayName','EKF-Output',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.8); 
plot(t,UKFOutput_sysID(:,2),'DisplayName','UKF-sysID-Output',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2);
plot(t,UKFAddOutput(:,2),'DisplayName','UKF-Add-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',2)
plot(t,SRUKFOutput(:,2),'DisplayName','SR-UKF-Output',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1); 
% plot(t,UKFAugOutput(:,2),'DisplayName','UKF-Aug-Output',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% plot(t,UKFFullyAugOutput(:,2),'DisplayName','UKF-Fully-Aug-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1); 
% plot(t,cUKFNLPOutput(:,2),'DisplayName','cUKF-NLP-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1); 
% plot(t,cUKFQPOutput(:,2),'DisplayName','cUKF-QP-Output',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,CKFOutput(:,2),'DisplayName','CKF-Output',...
%      'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2)
% ylim([0.4,0.85])
ylabel('S_{co2} [kg/m3]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_bac:  
subplot(3,1,3)
scatter(tMeas, MESS.yMeas(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,3),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% plot(t,EKFOutput(:,3),'DisplayName','EKF-Output',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.8); 
plot(t,UKFOutput_sysID(:,3),'DisplayName','UKF-sysID-Output',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2);
plot(t,UKFAddOutput(:,3),'DisplayName','UKF-Add-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',2)
plot(t,SRUKFOutput(:,3),'DisplayName','SR-UKF-Output',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1);
% plot(t,UKFAugOutput(:,3),'DisplayName','UKF-Aug-Output',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% plot(t,UKFFullyAugOutput(:,3),'DisplayName','UKF-Fully-Aug-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1);  
% plot(t,cUKFNLPOutput(:,3),'DisplayName','cUKF-NLP-Output',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1);
% plot(t,cUKFQPOutput(:,3),'DisplayName','cUKF-QP-Output',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,CKFOutput(:,3),'DisplayName','CKF-Output',...
%      'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2);  
% ylim([0.4,0.85])
ylabel('microbial biomass [kg/m3]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of UKF and clean model output')

%% Biomass X_bac:
% sigmaBacArray = sqrt(COVARIANCEUKF(6,6,:)); 
% sigmaBac = reshape(sigmaBacArray,nSamples + 1,1);
% 
% figure
% plot(tMeas,STATES(:,6)','DisplayName','true',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% hold on
% plot(t,ESTIMATESUKF(:,6)','DisplayName','estimate',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
% plot(t,ESTIMATESUKF(:,6)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
% plot(t,ESTIMATESUKF(:,6)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
% hold off
% ylabel('biomass X_{bac} [kg/m3]')
% title('Estimated and true biomass concentration')
% legend()

%% Plot trajectories of relevant states: 
figure()

% X_ch:
subplot(3,1,1)
plot(tMeas,trueStates(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
% plot(t,ESTIMATESEKF(:,3),'DisplayName','EKF',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
plot(t,ESTIMATESUKF_sysID(:,3),'DisplayName','UKF-sysID',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2);
plot(t,ESTIMATESUKFAdd(:,3),'DisplayName','UKF-Add',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
plot(t,ESTIMATESSRUKF(:,3),'DisplayName','SR-UKF',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1);
% plot(t,ESTIMATESUKFAug(:,3),'DisplayName','UKF-Aug',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
% plot(t,ESTIMATESUKFFullyAug(:,3),'DisplayName','UKF-Fully-Aug',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1);
% plot(t,ESTIMATEScUKFNLP(:,3),'DisplayName','cUKF-NLP',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1);
% plot(t,ESTIMATEScUKFQP(:,3),'DisplayName','cUKF-QP',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,ESTIMATESCKF(:,3),'DisplayName','CKF',...
%      'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6); 
% ylim([0.4,0.7])
ylabel('X_{ch} [kg/m3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
legend('Location','NorthEast'); 

% X_pr:
subplot(3,1,2)
plot(tMeas,trueStates(:,4),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
% plot(t,ESTIMATESEKF(:,4),'DisplayName','EKF',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
plot(t,ESTIMATESUKF_sysID(:,4),'DisplayName','UKF-sysID',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2);
plot(t,ESTIMATESUKFAdd(:,4),'DisplayName','UKF-Add',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
plot(t,ESTIMATESSRUKF(:,4),'DisplayName','SR-UKF',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1);
% plot(t,ESTIMATESUKFAug(:,4),'DisplayName','UKF-Aug',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
% plot(t,ESTIMATESUKFFullyAug(:,4),'DisplayName','UKF-Fully-Aug',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1);
% plot(t,ESTIMATEScUKFNLP(:,4),'DisplayName','cUKF-NLP',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1);
% plot(t,ESTIMATEScUKFQP(:,4),'DisplayName','cUKF-QP',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,ESTIMATESCKF(:,4),'DisplayName','CKF',...
%      'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6);
% ylim([0.4,0.7])
ylabel('X_{pr} [kg/m3]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_li:
subplot(3,1,3)
plot(tMeas,trueStates(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
% plot(t,ESTIMATESEKF(:,5),'DisplayName','EKF',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
plot(t,ESTIMATESUKF_sysID(:,5),'DisplayName','UKF-sysID',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',1.2);
plot(t,ESTIMATESUKFAdd(:,5),'DisplayName','UKF-Add',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
plot(t,ESTIMATESSRUKF(:,5),'DisplayName','SR-UKF',...
     'LineStyle','--', 'Color', colorPaletteHex(5), 'LineWidth',1);
% plot(t,ESTIMATESUKFAug(:,5),'DisplayName','UKF-Aug',...
%      'LineStyle',':', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
% plot(t,ESTIMATESUKFFullyAug(:,5),'DisplayName','UKF-Fully-Aug',...
%      'LineStyle','--', 'Color', colorPaletteHex(1), 'LineWidth',1);
% plot(t,ESTIMATEScUKFNLP(:,5),'DisplayName','cUKF-NLP',...
%      'LineStyle','--', 'Color', colorPaletteHex(4), 'LineWidth',1);
% plot(t,ESTIMATEScUKFQP(:,5),'DisplayName','cUKF-QP',...
%      'LineStyle','-.', 'Color', colorPaletteHex(1), 'LineWidth',1.4); 
% plot(t,ESTIMATESCKF(:,5),'DisplayName','CKF',...
%      'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6); 
% ylim([0.4,0.7])
ylabel('X_{li} [kg/m3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of true and UKF-estimates of states')

%% create plots for ECC paper
% one representative measurable (S_co2), one non-measurable state (X_ch): 

% set up color palette for 6 plots in magma style, where the first one is 
% omitted since it is too light. Sorted from light to dark.
eccColorPalette = ["#fcfdbf","#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"]; 
colorFeeding = "#02C39A"; % light green 

ECCPlot = figure;

% S_co2:
subplot(2,1,1)
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorFeeding, 'LineWidth',2); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
yyaxis left % now plot all the curves
plot(tMeas,yClean(:,2),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', eccColorPalette(3), 'LineWidth',1.5)
hold on; 
scatter(tMeas, MESS.yMeas(:,2),70,'DisplayName','noisy measurements',...
        'Marker','.', 'Color', eccColorPalette(6)); 
% plot(t,EKFOutput(:,2),'DisplayName','EKF',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.8); 
plot(t,UKFOutput_sysID(:,2),'DisplayName','UKF-sysID',...
     'LineStyle',':', 'Color', eccColorPalette(4), 'LineWidth',1) % 2
% plot(t,UKFAddOutput(:,2),'DisplayName','UKF-Add',...
%      'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1.5)
% plot(t,SRUKFOutput(:,2),'DisplayName','SR-UKF',...
%      'LineStyle','--', 'Marker','none', 'Color', eccColorPalette(5), 'LineWidth',1); 
% plot(t,UKFAugOutput(:,2),'DisplayName','UKF-Aug',...
%      'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1); 
% plot(t,UKFFullyAugOutput(:,2),'DisplayName','UKF-Fully-Aug','Marker','none',...
%      'LineStyle','-.', 'Color', eccColorPalette(5), 'LineWidth',1.5); 
plot(t,cUKFNLPOutput(:,2),'DisplayName','cUKF-NLP',...
     'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1); 
plot(t,cUKFQPOutput(:,2),'DisplayName','cUKF-QP','Marker','none',...
     'LineStyle','-.', 'Color', eccColorPalette(5), 'LineWidth',1.5); 
% ylim([0.4,0.85])
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('S_{co2} [kg/m^3]')
legend('Location','NorthEast'); 

% X_ch:
subplot(2,1,2)
yyaxis right    % plot feeding first
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorFeeding, 'LineWidth',2); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
yyaxis left     % now plot all the curves
plot(tMeas,trueStates(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', eccColorPalette(3), 'LineWidth',1.5); 
hold on; 
% plot(t,ESTIMATESEKF(:,3),'DisplayName','EKF',...
%      'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
plot(t,ESTIMATESUKF_sysID(:,3),'DisplayName','UKF-sysID',...
     'LineStyle',':', 'Color', eccColorPalette(4), 'LineWidth',1); % 2
% plot(t,ESTIMATESUKFAdd(:,3),'DisplayName','UKF-Add',...
%      'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1.5);
% plot(t,ESTIMATESSRUKF(:,3),'DisplayName','SR-UKF',...
%      'LineStyle','--', 'Color', eccColorPalette(5), 'LineWidth',1);
% plot(t,ESTIMATESUKFAug(:,3),'DisplayName','UKF-Aug',...
%      'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1);
% plot(t,ESTIMATESUKFFullyAug(:,3),'DisplayName','UKF-Fully-Aug',...
%      'LineStyle','-.', 'Color', eccColorPalette(5), 'LineWidth',1.5);
plot(t,ESTIMATEScUKFNLP(:,3),'DisplayName','cUKF-NLP',...
     'LineStyle','-', 'Color', eccColorPalette(2), 'LineWidth',1);
plot(t,ESTIMATEScUKFQP(:,3),'DisplayName','cUKF-QP',...
     'LineStyle','-.', 'Color', eccColorPalette(5), 'LineWidth',1.5); 
% ylim([0.4,0.7])
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('X_{ch} [kg/m^3]')
xlabel('time [d]')

%% create a plot of nRMSE over tCalc

% tCalc1 = mean(runTime); 
% nRMSE1 = mean([mean_nRMSE1y;mean_nRMSE1x]); 
% 
% % alg1 = [tCalc1, nRMSE1]; 
% 
% figure_t_RMSE = figure;
% scatter()

