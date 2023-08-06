%% Version
% (R2022b) Update 5
% Erstelldatum: 23.6.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R3-frac

close all
clear all
clc

global counter
counter = 0; 

% Load Measurement Data:
load('generatedOutput\Messung_ADM1_R3_frac.mat')

%% Initialization and Tuning of EKF 
% initial state values:
x0Init = MESS.x0;       % intitial value from createMess-file 
x0 = x0Init;            % x0 will be replaced in every iteration later
% false initial estimate (add random number but ensure positive concentrations)
nStates = length(x0); 
nTheta = length(params.th); % number of time-variant parameters
rng('default');     % fix seed for random number generation (for replicable results)
xHat = x0Init.*abs(randn(nStates,1)); 
xMinus = xHat;      % to be overwritten
% xMinus = x0Init;    % XY Rania

tMeas = MESS.t;
nSamples = length(tMeas);  % number of measurements taken
% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMeas); 
dt = diff_t(1);             % sampling time
t = [tMeas;tMeas(end)+dt];  % add one time interval at end

% allocate memory for raw matrices. all matrices requiring initialization 
% have 1 column more than the number of measurement instances:
MEAS = MESS.yMeas;
STATES = MESS.x; 
ESTIMATES = zeros(nSamples + 1,nStates);
COVARIANCE = zeros(nStates,nStates,nSamples + 1);
GAIN = zeros(nSamples,nStates);

% Initialize Kalman Filter:
ESTIMATES(1,:) = xHat;
P0 = eye(nStates); % XY: sicher besseres Tuning möglich
COVARIANCE(:,:,1) = P0; 
POld = P0;      % to overwrite
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
% Tune Kalman Filter: measurement uncertainty: 
R = buffer * MESS.C; 
% Tune Kalman Filter: process uncertainty: 
Q = diag(x0Init);

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
thS = sym('th', [nTheta,1]);    % time-variant parameters (theta)
cS = sym('c', [size(params.c)]);% time-invariant parameters (known & constant)
aS = sym('a', [size(params.a)]); % petersen matrix with stoichiometric constants

% obtain symbolic objects of model equations:
dynamics = ADM1_R3_frac_ode_sym(xS, uS, xiS, thS, cS, aS); 
outputs = ADM1_R3_frac_mgl_sym(xS,cS);

% derive function handles of ODEs and output equation: 
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% partial derivatives for EKF (symbolic): 
dfdxSym = jacobian(dynamics,xS); 
dhdxSym = jacobian(outputs,xS); 

% convert to numeric function handles: 
dfdx = matlabFunction(dfdxSym, 'Vars', {xS, uS, thS, cS, aS}); 
dhdx = matlabFunction(dhdxSym, 'Vars', {xS, cS}); 

tic
% integrate across all (online) measurement intervals (like in reality):
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
    %% execute EKF

    [xPlus,PPlus,Kv] = extendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx); % standard EKF
%     [xPlus,PPlus,Kv] = constrainedExtendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC,R); % constrained EKF
    
    % save results:
    ESTIMATES(k+1,:) = xPlus';
    COVARIANCE(:,:,k+1) = PPlus; 
    GAIN(k,:) = Kv;     % Kalman Gain * Innovation
    
    % Update for next iteration:  
    xMinus = xPlus;     % estimated state from Kalman Filter
    POld = PPlus;       % state error covariance matrix
end
toc

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
EKFOutput = nan(nSamples + 1,q);   % allocate memory
for k = 1:nSamples + 1
    EKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
end
yClean = MESS.yClean; 
feedVolFlow = inputMat(:,2);    % [l/d]

%% compute goodness of fit for all measurements
RMSSE = zeros(q,1);         % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = EKFOutput(2:end,kk);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

RMSSE_mean = mean(RMSSE); 

%% Plot results

% plot model output based on EKF estimation and compare with real
% measurements:
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figure

% gas volume flow: 
subplot(4,2,1)
scatter(tMeas, MESS.yMeas(:,1)/24,'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,1)/24,'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,EKFOutput(:,1)/24,'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([-10,50])
ylabel('gas vol flow [l/h]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% pch4: 
subplot(4,2,2)
scatter(tMeas, MESS.yMeas(:,2),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on
plot(tMeas,yClean(:,2),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,EKFOutput(:,2),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.85])
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% pco2:
subplot(4,2,3)
scatter(tMeas, MESS.yMeas(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,3),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,3),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
% ylim([0.2,0.6])
ylabel('p_{co2} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% pH:  
subplot(4,2,4)
scatter(tMeas, MESS.yMeas(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,4),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,4),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('pH [-]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% SIN:
subplot(4,2,5)
scatter(tMeas, MESS.yMeas(:,5),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,5),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,5),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('inorg. nitrogen [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% TS:  
subplot(4,2,6)
scatter(tMeas, MESS.yMeas(:,6),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,6),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,6),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0,0.5])
ylabel('total solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% VS:  
subplot(4,2,7)
scatter(tMeas, MESS.yMeas(:,7),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,7),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,7),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0.4,0.8])
ylabel('volatile solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% Sac:  
subplot(4,2,8)
scatter(tMeas, MESS.yMeas(:,8),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,8),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,EKFOutput(:,8),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0,2.5])
ylabel('acetic acid [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of EKF and clean model output')

%% save results in struct EKF: 
EKF.t = t; 
EKF.tMeas = tMeas; 
EKF.xHat = ESTIMATES; 
EKF.yHat = EKFOutput;
EKF.x = STATES;
EKF.y = yClean;  
EKF.inputMat = inputMat;    % u in [l/d]
EKF.R = R;
EKF.Q = Q; 

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
mkdir(pathToResults);   % create subfolder (gives warning if exists yet)
fileName = 'EKF_R3_frac.mat'; 
save(fullfile(pathToResults,fileName), 'EKF', 'params')

%% Plot trajectories of relevant states: 
trueStates = MESS.x; 

figure()

% S_IN:
subplot(3,2,1)
plot(tMeas,trueStates(:,4),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,4),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('S_{IN} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% S_h2o:
subplot(3,2,2)
plot(tMeas,trueStates(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,5),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.8,1.1])
ylabel('S_{h2o} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% X_ch_fast:
subplot(3,2,3)
plot(tMeas,trueStates(:,6),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,6),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,fast} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_ch_slow:
subplot(3,2,4)
plot(tMeas,trueStates(:,7),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,7),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,slow} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_bac:
subplot(3,2,5)
plot(tMeas,trueStates(:,10),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,11),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{bac} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_ac:
subplot(3,2,6)
plot(tMeas,trueStates(:,11),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,11),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ac} [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of true and EKF-estimates of states')