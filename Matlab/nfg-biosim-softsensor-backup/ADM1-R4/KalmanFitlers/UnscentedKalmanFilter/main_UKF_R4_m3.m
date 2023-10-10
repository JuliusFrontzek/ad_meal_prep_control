%% Version
% (R2022b) Update 5
% Erstelldatum: 05.10.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-m3
% mit Konzentrationen in kg/m³ und Volumina in m³

addpath('generatedOutput/')
addpath('modelEquations/')

close all
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
xHat = x0Init.*(ones(nStates,1)) + 0.1*abs(randn(nStates,1)); 
xMinus = xHat;      % to be overwritten
% xMinus = x0Init;    % XY Rania

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
ESTIMATES = zeros(nSamples + 1,nStates);
COVARIANCE = zeros(nStates,nStates,nSamples + 1);
% GAIN = zeros(nSamples,nStates);

% Initialize Kalman Filter:
ESTIMATES(1,:) = xHat;
% P0 = eye(nStates); % XY: sicher besseres Tuning möglich
P0 = diag((xHat-x0).^2);    % Schneider und Georgakis 
COVARIANCE(:,:,1) = P0; 
PMinus = P0;      % to overwrite
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
% Tune Kalman Filter: measurement uncertainty: 
R(4,4) = 1E-2*R(4,4);    % SIN
% R(5,5) = 5E2*R(5,5);    % TS
% Tune Kalman Filter: process uncertainty: 
Q = diag([0.016, 0.555, 5.563, 958.4, 1.263, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
% Q = eye(nStates); 

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
    %% execute UKF

%     [xPlus,PPlus,Kv] = extendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdxNew,dhdxNew); % standard EKF
    [xPlus,PPlus] = unscKalmanFilterKolasAdditive(xMinus,PMinus,...
                            tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPlus,PPlus] = unscKalmanFilterKolasAugmented(xMinus,PMinus,...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);
%     [xPlus,PPlus] = constrUnscKalmanFilterKolasAdditive(xMinus,PMinus, ...
%                             tSpan,feedInfo,yMeas,params,Q,R,f,g);

    % save results:
    ESTIMATES(k+1,:) = xPlus';
    COVARIANCE(:,:,k+1) = PPlus; 
%     GAIN(k,:) = Kv;     % Kalman Gain * Innovation
    
    % Update for next iteration:  
    xMinus = xPlus;     % estimated state from Kalman Filter
    PMinus = PPlus;       % state error covariance matrix
end
toc

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
UKFOutput = nan(nSamples + 1,q);   % allocate memory
for k = 1:nSamples + 1
    UKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
end
yClean = MESS.yClean; 
feedVolFlow = inputMat(:,2);    % [m³/d]

%% compute goodness of fit for all measurements
RMSSE = zeros(q,1);         % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = UKFOutput(2:end,kk);    
    
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
subplot(2,2,1)
scatter(tMeas, MESS.yMeas(:,1),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,1),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,UKFOutput(:,1),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0,0.3])
ylabel('gas vol flow [m³/d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% pch4: 
subplot(2,2,2)
scatter(tMeas, MESS.yMeas(:,2),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on
plot(tMeas,yClean(:,2),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,UKFOutput(:,2),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0.4,0.85])
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% pco2:
subplot(2,2,3)
scatter(tMeas, MESS.yMeas(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,3),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,3),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0.2,0.6])
ylabel('p_{co2} [bar]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% SIN:  
subplot(2,2,4)
scatter(tMeas, MESS.yMeas(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,4),'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(t,UKFOutput(:,4),'DisplayName','UKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('inorg. nitrogen [g/L]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [m³/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of UKF and clean model output')

%% Biomass X_bac:
sigmaBacArray = sqrt(COVARIANCE(7,7,:)); 
sigmaBac = reshape(sigmaBacArray,nSamples + 1,1);

figure
plot(tMeas,STATES(:,7)','DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
hold on
plot(t,ESTIMATES(:,7)','DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,ESTIMATES(:,7)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
plot(t,ESTIMATES(:,7)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
hold off
ylabel('biomass X_{bac} [g/L]')
title('Estimated and true biomass concentration')
legend()

%% Plot trajectories of relevant states: 
trueStates = MESS.x; 

figure()

% S_IN:
subplot(3,1,1)
scatter(tMeas, MESS.yMeas(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,trueStates(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,ESTIMATES(:,3),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('S_{IN} [g/L]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_ch:
subplot(3,1,2)
plot(tMeas,trueStates(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,5),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_bac:
subplot(3,1,3)
plot(tMeas,trueStates(:,8),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,8),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{bac} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 

sgtitle('Comparison of true and UKF-estimates of states')

