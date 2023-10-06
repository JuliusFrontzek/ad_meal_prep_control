%% Version
% (R2022b) Update 5
% Erstelldatum: 05.10.2023
% Autor: Simon Hellmann

%% Unscented und Cubature Kalman Filter für ADM1-R4-m3
% mit Konzentrationen in kg/m³ und Volumina in m³

addpath('generatedOutput/')
addpath('modelEquations/')

% close all
clear all
clc

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX
counterSigmaInit = 0;
counterSigmaProp = 0;
counterSigmaX = 0; 
counterX = 0; 

% Load Measurement Data:
load Messung_ADM1_R4_m3

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
ESTIMATESUKF = zeros(nSamples + 1,nStates);
COVARIANCEUKF = zeros(nStates,nStates,nSamples + 1);
% GAIN = zeros(nSamples,nStates);
ESTIMATESCKF = zeros(nSamples + 1,nStates);
ESTIMATESEKF = zeros(nSamples + 1,nStates);
COVARIANCECKF = zeros(nStates,nStates,nSamples + 1);
COVARIANCEEKF = zeros(nStates,nStates,nSamples + 1);

% Initialize Kalman Filters:
xMinusUKF = xHat;         % to be overwritten
xMinusCKF = xMinusUKF; 
xMinusEKF = xMinusUKF; 
% P0 = eye(nStates); % XY: sicher besseres Tuning möglich
P0 = diag((xHat-x0).^2);    % Schneider und Georgakis 
PMinusUKF = P0;             % to overwrite
PMinusCKF = PMinusUKF; 
PMinusEKF = PMinusUKF; 
ESTIMATESUKF(1,:) = xMinusUKF;
COVARIANCEUKF(:,:,1) = PMinusUKF; 
ESTIMATESCKF(1,:) = xMinusCKF;
COVARIANCECKF(:,:,1) = PMinusUKF; 
ESTIMATESEKF(1,:) = xMinusCKF;
COVARIANCEEKF(:,:,1) = PMinusUKF; 
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
% R = eye(size(MESS.C)); 
% Tune Kalman Filter: 
% measurement uncertainty: 
% R(1,1) = 1E2*R(1,1);    % gasVolFlow
% R(4,4) = 1E2*R(4,4);    % SIN
% process uncertainty: 
% Q = diag([0.016, 0.555, 0.563, 958.4, 1.263, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
Q = eye(nStates); 
% increase trust in S_IN, S_h2o and X_ash from model:
% Q([3,4,9],[3,4,9]) = 1E-3*eye(3); 
% Q(8,8) = 1E-1; 

% obtain feeding information:
inputMat = MESS.inputMat;   % [tEvents,feed vol flow,inlet concentrations]
tEvents = inputMat(:,1);    % times of feeding events (on/off)

% if inputMat has first feeding entry only at t > t0, add zero row at t0:
if tEvents(1) > 0
    nColcInputMat = size(inputMat,2);   % number of columns in inputMat
    inputMat = [zeros(nColcInputMat,1),inputMat]; 
end

%% derive all function handles for EKF from symbolic model equations
% define symbolic variables ("S") (all vectors are column vectors):
xS = sym('x', [nStates,1]);     % states as col. vector
syms uS real                    % input
xiS = sym('xi', [nStates,1]);   % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));    % time-variant parameters (theta)
cS = sym('c', size(params.c));          % time-invariant parameters (known & constant)
aS = sym('a', size(params.a)); % petersen matrix with stoichiometric constants

% obtain symbolic objects of model equations:
dynamics = BMR4_AB_ode_sym(xS, uS, xiS, thS, cS, aS); 
outputs = BMR4_AB_mgl_sym(xS,cS);

% derive function handles of ODEs and output equation: 
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% partial derivatives for EKF (symbolic): 
dfdxSym = jacobian(dynamics,xS); 
dhdxSym = jacobian(outputs,xS); 

% convert to numeric function handles: 
dfdx = matlabFunction(dfdxSym, 'Vars', {xS, uS, thS, cS, aS}); 
dhdx = matlabFunction(dhdxSym, 'Vars', {xS, cS}); 

%% create CKF object and tune it
% Set up Cubature Kalman Filter
% @measurementFcnADM1R4
% @measurementFcnADM1R4_reduced
Q_CKF = eye(nStates);
R_CKF = buffer * MESS.C;  
% R_CKF = R(2:6,2:6); 
ckf = trackingCKF(@stateTransitionFcnADM1R4,@measurementFcnADM1R4,xMinusCKF, ...
    'HasAdditiveProcessNoise',true, 'HasAdditiveMeasurementNoise',true, ...
    'MeasurementNoise',R_CKF, 'ProcessNoise',Q_CKF); % , 'EnableSmoothing',true

%% integrate across all (online) measurement intervals (like in reality):
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
    %% execute KFs
    
    % call EKF:
    [xPlusEKF,PPlusEKF,~] = extendedKalmanFilter(xMinusEKF,PMinusEKF, ...
                        tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx);
    
    % call UKF:
    [xPlusUKF,PPlusUKF] = unscKalmanFilterKolasAdditive(xMinusUKF,PMinusUKF,...
                            tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPlusUKF,PPlusUKF] = unscKalmanFilterKolasAugmented(xMinusUKF,PMinusUKF,...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPlusUKF,PPlusUKF] = unscKalmanFilterKolasFullyAugmented(xMinusUKF,PMinusUKF, ...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g); 

    % call CKF:
    [xPlusCKF,PPlusCKF] = my_CKF_ADM1(ckf,feedInfo,yMeas,tSpan,params,f,g);
    
    % save results:
    ESTIMATESUKF(k+1,:) = xPlusUKF';
    COVARIANCEUKF(:,:,k+1) = PPlusUKF; 
%     GAIN(k,:) = Kv;     % Kalman Gain * Innovation
    ESTIMATESCKF(k+1,:) = xPlusCKF';
    COVARIANCECKF(:,:,k+1) = PPlusCKF; 
    ESTIMATESEKF(k+1,:) = xPlusEKF';
    COVARIANCEEKF(:,:,k+1) = PPlusEKF; 

    % Update for next iteration:  
    xMinusUKF = xPlusUKF;      % estimated state from Kalman Filter
    PMinusUKF = PPlusUKF;      % state error covariance matrix
    % same for CKF: 
    xMinusCKF = xPlusCKF; 
    PMinusCKF = PPlusCKF; 
    % same for EKF: 
    xMinusEKF = xPlusEKF; 
    PMinusEKF = PPlusEKF; 

end
toc

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
UKFOutput = nan(nSamples + 1,q);    % allocate memory
CKFOutput = nan(size(UKFOutput));   % allocate memory
EKFOutput = nan(size(UKFOutput));   % allocate memory
for k = 1:nSamples + 1
    UKFOutput(k,:) = g(ESTIMATESUKF(k,:)',params.c);
    CKFOutput(k,:) = g(ESTIMATESCKF(k,:)',params.c);
    EKFOutput(k,:) = g(ESTIMATESEKF(k,:)',params.c);
end
yClean = MESS.yClean; 
feedVolFlow = inputMat(:,2);    % [l/d]

%% compute goodness of fit for all measurements
RMSSE = zeros(q,1);         % allocate memory
RMSSE_UKF_CKF_y = zeros(q,1); % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurementsUKF = UKFOutput(2:end,kk);    
    estimatedMeasurementsCKF = CKFOutput(2:end,kk);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk)] = computeRMSSE(measurements,estimatedMeasurementsUKF); 
    
    % get same between UKF and CKF: 
    RMSSE_UKF_CKF_y(kk) = computeRMSSE(estimatedMeasurementsUKF,estimatedMeasurementsCKF); 
end

RMSSE_mean = mean(RMSSE); 

%% compute RMSSE of states between ...
% ... CKF and UKF: 
RMSSE_UKF_CKF_x = zeros(nStates,1); 
for ll = 1:nStates
    RMSSE_UKF_CKF_x(ll) = computeRMSSE(ESTIMATESUKF(:,ll),ESTIMATESCKF(:,ll)); 
end
% ... UKF and true states: 
RMSSE_UKF_true_x = zeros(nStates,1); 
for jj = 1:nStates
    RMSSE_UKF_true_x(jj) = computeRMSSE(ESTIMATESUKF(2:end,jj),STATES(:,jj)); 
end
% ... CKF and true states: 
RMSSE_CKF_true_x = zeros(nStates,1); 
for lll = 1:nStates
    RMSSE_CKF_true_x(lll) = computeRMSSE(ESTIMATESCKF(2:end,lll),STATES(:,lll)); 
end
... EKF and true states: 
RMSSE_EKF_true_x = zeros(nStates,1); 
for nn = 1:nStates
    RMSSE_EKF_true_x(nn) = computeRMSSE(ESTIMATESEKF(2:end,nn),STATES(:,nn)); 
end

%% Plot results

% plot model output based on EKF estimation and compare with real
% measurements:
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figure

% gas volume flow: 
subplot(3,2,1)
scatter(tMeas, MESS.yMeas(:,1),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,1),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,UKFOutput(:,1),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% plot(t,CKFOutput(:,1),'DisplayName','CKF-Output',...
%      'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
plot(t,EKFOutput(:,1),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylim([0.1,0.2])
ylabel('gas vol flow [m3/d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tMeas, MESS.yMeas(:,2),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on
plot(tMeas,yClean(:,2),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,UKFOutput(:,2),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% plot(t,CKFOutput(:,2),'DisplayName','CKF-Output',...
%      'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
plot(t,EKFOutput(:,2),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylim([0.5,0.6])
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% pco2:
subplot(3,2,3)
scatter(tMeas, MESS.yMeas(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,3),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,3),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,CKFOutput(:,3),'DisplayName','CKF-Output',...
     'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2)
plot(t,EKFOutput(:,3),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylim([0.45,0.55])
ylabel('p_{co2} [bar]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% SIN:  
subplot(3,2,4)
scatter(tMeas, MESS.yMeas(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,4),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,4),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,CKFOutput(:,4),'DisplayName','CKF-Output',...
     'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
plot(t,EKFOutput(:,4),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylabel('inorg. nitrogen [kg/m³]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% TS:  
subplot(3,2,5)
scatter(tMeas, MESS.yMeas(:,5),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,5),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,5),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,CKFOutput(:,5),'DisplayName','CKF-Output',...
     'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
plot(t,EKFOutput(:,5),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylim([0,0.1])
ylabel('total solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% VS:  
subplot(3,2,6)
scatter(tMeas, MESS.yMeas(:,6),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,6),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,6),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,CKFOutput(:,6),'DisplayName','CKF-Output',...
     'LineStyle','-.', 'Color', 'magenta', 'LineWidth',1.2); 
plot(t,EKFOutput(:,6),'DisplayName','EKF-Output',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',0.6); 
ylim([0,1])
ylabel('volatile solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of UKF and clean model output')

%% Biomass X_bac:
sigmaBacArray = sqrt(COVARIANCEUKF(7,7,:)); 
sigmaBac = reshape(sigmaBacArray,nSamples + 1,1);

figure
plot(tMeas,STATES(:,8)','DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
hold on
plot(t,ESTIMATESUKF(:,8)','DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,ESTIMATESUKF(:,8)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
plot(t,ESTIMATESUKF(:,8)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
hold off
ylabel('biomass X_{bac} [kg/m³]')
title('Estimated and true biomass concentration')
legend()

%% Plot trajectories of relevant states: 
trueStates = MESS.x; 

figure()

% S_IN:
subplot(3,1,1)
scatter(tMeas, MESS.yMeas(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on;
plot(tMeas,trueStates(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);  
plot(t,ESTIMATESUKF(:,3),'DisplayName','UKF',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
plot(t,ESTIMATESCKF(:,3),'DisplayName','CKF',...
     'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6); 
plot(t,ESTIMATESEKF(:,3),'DisplayName','EKF',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',1); 
% ylim([0.4,0.7])
ylabel('S_{IN} [g/L]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% X_ch:
subplot(3,1,2)
plot(tMeas,trueStates(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATESUKF(:,5),'DisplayName','UKF',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
plot(t,ESTIMATESCKF(:,5),'DisplayName','CKF',...
     'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6);
plot(t,ESTIMATESEKF(:,5),'DisplayName','EKF',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',1); 
% ylim([0.4,0.7])
ylabel('X_{ch} [kg/m³]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_bac:
subplot(3,1,3)
plot(tMeas,trueStates(:,8),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATESUKF(:,8),'DisplayName','UKF',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
plot(t,ESTIMATESCKF(:,8),'DisplayName','CKF',...
     'LineStyle','--', 'Color', 'magenta', 'LineWidth',0.6); 
plot(t,ESTIMATESEKF(:,8),'DisplayName','EKF',...
     'LineStyle',':', 'Color', 'red', 'LineWidth',1); 
% ylim([0.4,0.7])
ylabel('X_{bac} [kg/m³]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of true and UKF-estimates of states')

