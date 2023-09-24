%% Version
% (R2022b) Update 5
% Erstelldatum: 22.09.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-frac mit Online/Offline Messwerten und 
% multiple augmentation

close all
clear all
clc

global counter
counter = 0; 

% load steady-state values from Sören's implementation: 
load('generatedOutput/Messung_ADM1_R4_fracIntensiveSampling_withCampaigns.mat')

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

tMinor = MESS.tOnline;          % time grid of minor measurements 
nSamplesMinor = length(tMinor);  % number of online measurements taken
% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMinor); 
dt = diff_t(1);                 % sampling time
tKF = [tMinor;tMinor(end)+dt];  % time vector for Kalman Filter

% create unified measurement matrix in which all measurements go at their
% respective arrival times:
MEASOn = MESS.yMeasOn;
MEASStd = MESS.yMeasStdArrival; 
MEASCamp = MESS.yMeasCampArrival; 
q = size(MESS.yMeas,2);     % # signals all together 
qOn = size(MEASOn,2);       % # secondary (online) measurement signals
qStd = size(MEASStd,2);     % # offline measurement signals (not measured in campaigns)
qCamp = size(MEASCamp,2);   % # offline measurement signals (measured in campaigns)
MEASUnite = nan(nSamplesMinor,q);     % allocate memory 
% insert minor instances in first qOn columns of MEASUnite: 
MEASUnite(:,1:qOn) = MEASOn; 
% determine correct row positioning of major instances and insert respective
% measurements:
tStdArrivalShift = MESS.tStdArrivalShift; 
tCampArrivalShift = MESS.tCampArrivalShift;
[~,idxTStdArrival] = ismember(tStdArrivalShift,tMinor); 
[~,idxTCampArrival] = ismember(tCampArrivalShift,tMinor); 
MEASUnite(idxTCampArrival,qOn+1:qOn+qCamp) = MEASCamp;      % Nitrogen measurements 
MEASUnite(idxTStdArrival,qOn+qCamp+1:end) = MEASStd;        % other std. measurements

% create united time vectors of whenever offline measurements are sampled:
tStdSampleShift = MESS.tStdSampleShift;
tCampSampleShift = MESS.tCampSampleShift; 
tOfflineSampleShift = unique([tStdSampleShift;tCampSampleShift]);
% ... or returned:
tOfflineArrivalShift = unique([tStdArrivalShift;tCampArrivalShift]);

%% set up raw matrices and structure for Kalman Filterin
% Note: all matrices which require initialization have 1 column more than the 
% number of measurement instances.
STATES = MESS.xSolOn; 
ESTIMATES = zeros(nSamplesMinor + 1,nStates);
COVARIANCE = zeros(nStates,nStates,nSamplesMinor + 1);

% Initialize Kalman Filter:
ESTIMATES(1,:) = xHat;
P0 = eye(nStates); % XY: sicher besseres Tuning möglich
COVARIANCE(:,:,1) = P0; 
PMinus = P0;      % to overwrite
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
% Hand-Tune Kalman Filter: measurement uncertainty: 
R(4,4) = 5E3*R(4,4);    % SIN
R(5,5) = 5E2*R(5,5);    % TS
% Tune Kalman Filter: process uncertainty: 
Q = diag([0.016, 0.555, 0.563, 958.4, 1.263, 2.816, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
% Q(4,4) = 1E-3*Q(4,4);       % unit change for h2o [g/l] --> [kg/l]
% Q(10,10) = 1E-3*Q(10,10);   % unit change for ash [g/l] --> [kg/l]  

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
thS = sym('th', size(params.th));   % time-variant parameters (theta)
cS = sym('c', size(params.c));      % time-invariant parameters (known & constant)
aS = sym('a', size(params.a));  % petersen matrix with stoichiometric constants

% obtain symbolic objects of model equations:
dynamics = BMR4_AB_frac_ode_sym(xS, uS, xiS, thS, cS, aS); 
outputs = BMR4_AB_frac_mgl_sym(xS,cS);

% derive function handles of ODEs and output equation: 
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% partial derivatives for EKF (symbolic): 
dfdxSym = jacobian(dynamics,xS); 
dhdxSym = jacobian(outputs,xS); 

% convert to numeric function handles: 
dfdx = matlabFunction(dfdxSym, 'Vars', {xS, uS, thS, cS, aS}); 
dhdx = matlabFunction(dhdxSym, 'Vars', {xS, cS}); 

%% call EKF repeatedly
tic
% flagAugmented = 0;      % 0: non-augmented; 1: augmented
% flagDelayPeriod = 0;    % 0: not in delay period; 1: in delay period

nAug = 0;           % # state augmentations
flagArrival = 0;    % 0: no lab measurement returned in current time step. 1: lab measurement returned

% integrate across all (online) measurement intervals (like in reality):
for k = 1:nSamplesMinor 
    
    tSpan = [tKF(k) tKF(k+1)]; % measurement interval. In reality, this is the 
    tk = tSpan(1);      % t_k
    tkp1 = tSpan(2);    % t_k+1
        
    % check if you're AT any primary sampling (k==s): 
    if ismember(tk,tOfflineSampleShift)
        disp('sample got sent to lab'); 
        nAug = nAug + 1; 
        % augment and initialize state x & state err. cov. matrix P.
        % for this purpose, extract topmost core of x and P, then replicate
        xMinusCore = xMinus(1:nStates); 
        PMinusCore = PMinus(1:nStates,1:nStates); 
        xAugMinus = repmat(xMinusCore,1+nAug,1); 
        PAugMinus = repmat(PMinusCore,1+nAug,1+nAug); 
        % overwrite old (unaugmented) entities: 
        xMinus = xAugMinus; 
        PMinus = PAugMinus; 
%         flagAugmented = 1; 
%         flagDelayPeriod = 1; 
    end
    
    if nAug >= 1
        disp('waiting for some lab measurements...')
    end

    disp(['# augmentations currently: ',num2str(nAug)]);

    % get most recent measurement:
    yMeas = MEASUnite(k,:);    % simulated measurement
%     % check if you're at minor or major instance: 
%     flagMajor = all(~isnan(yMeas)); % 0: minor instance, 1: major instance
    % check if any lab measurement was returned:
    flagArrival = any(~isnan(yMeas(qOn+1:end))); % 0: minor instance, 1: major instance
    % XY: check, ob das stimmt

    if flagArrival == 1
        disp('just received lab measurements!')
    end    

    %% get feeding information:
    % pass only relevant feedings during the measurement interval, because
    % only those would be known in reality:
    idxRelEvents = find(tEvents >= tk & tEvents <= tkp1);
    tRelEvents = tEvents(idxRelEvents); % Auswertung anhand Index
    
    % find the critcal last feeding event before current measurement interval:
    idxLastEvent = find(tEvents < tk,1,'last');

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

%     [xPlus,PPlus,Kv] = extendedKalmanFilter(xMinus,PMinus,tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx); % standard EKF
%     [xPlus,PPlus,Kv] = constrainedExtendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC,R); % constrained EKF
    [xPlus,PPlus] = extendedKalmanFilterMultiRateCampaigns(xMinus,PMinus,feedInfo, ...
                        yMeas,params,Q,R,f,g,dfdx,dhdx, ...
                        tSpan,nStates,qOn,qStd,qCamp,nAug,flagArrival); % mulirate EKF with multiple augmentation
    
    % save results:
    ESTIMATES(k+1,:) = xPlus(1:nStates)';   % only keep the "real" states, not the sample state in case of augmentation
    COVARIANCE(:,:,k+1) = PPlus(1:nStates,1:nStates);   % same for P-matrix
    
    % Update for next iteration:  
    xMinus = xPlus;     % estimated state from Kalman Filter
    PMinus = PPlus;     % state error covariance matrix

    % remove state augmentation once primary measurement (& delay period) is over:
%     if ismember(tk,tOfflineArrivalShift)
    if flagArrival == 1
        disp(['primary measurement over', newline, ...
              'remove 1 augmentation']); 
        % remove augmentation of state and state err. cov. matrix P:
        nAug = nAug - 1; 
        xMinusOneAugLess = xMinus(1:(1+nAug)*nStates); 
        PMinusOneAugLess = PMinus(1:(1+nAug)*nStates,1:(1+nAug)*nStates); 
        % overwrite old (augmented) enteties: 
        xMinus = xMinusOneAugLess; 
        PMinus = PMinusOneAugLess; 
%         flagAugmented = 0; 
        flagArrival = 0; 
    end

    % XY: check, ob nAug und flagArrival Sinn machen! 
end
toc

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
EKFOutput = nan(nSamplesMinor + 1,q);   % allocate memory
for k = 1:nSamplesMinor + 1
    EKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
end

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
% plot model output based on EKF estimates, compare with real measurements:
yCleanOn = MESS.yCleanOn; 
yCleanStd = MESS.yCleanStd; % clean model outputs of std. measurements
yCleanCamp = MESS.yCleanCamp; % clean model outputs of std. measurements
feedVolFlow = inputMat(:,2);% [l/d]

colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figOutputs = figure;

% gas volume flow: 
subplot(3,2,1)
scatter(tMinor, MEASOn(:,1)/24,'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMinor,yCleanOn(:,1)/24,'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(tKF,EKFOutput(:,1)/24,'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0,15])
ylabel('gas vol flow [l/h]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight  % draw axis only as long as time vectors
% xlabel('time [d]')
legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tMinor, MEASOn(:,2),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on
plot(tMinor,yCleanOn(:,2),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(tKF,EKFOutput(:,2),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0.4,0.85])
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% pco2:
subplot(3,2,3)
scatter(tMinor, MEASOn(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMinor,yCleanOn(:,3),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,EKFOutput(:,3),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0.2,0.6])
ylabel('p_{co2} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight % draw axis only as long as time vectors
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% SIN:  
subplot(3,2,4)
scatter(tOfflineArrivalShift,MEASOff(:,1),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSampleShift; tMinor(end)],[yCleanOff(:,1);yCleanOff(end,1)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,EKFOutput(:,4),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('inorg. nitrogen [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight % draw axis only as long as time vectors
% xlabel('time [d]')
legend('Location','NorthEast'); 

% TS:  
subplot(3,2,5)
scatter(tOfflineArrivalShift,MEASOff(:,2),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSampleShift; tMinor(end)],[yCleanOff(:,2);yCleanOff(end,2)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,EKFOutput(:,5),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0,0.1])
ylabel('total solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
axis tight % draw axis only as long as time vectors
% legend('Location','NorthEast'); 

% VS:  
subplot(3,2,6)
scatter(tOfflineArrivalShift,MEASOff(:,3),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSampleShift; tMinor(end)],[yCleanOff(:,3);yCleanOff(end,3)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,EKFOutput(:,6),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylim([0,1])
ylabel('volatile solids [-]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
axis tight % draw axis only as long as time vectors
% legend('Location','NorthEast'); 

sgtitle('Comparison of EKF and clean model output')

%% Biomass X_bac:
sigmaBacArray = sqrt(COVARIANCE(9,9,:)); 
sigmaBac = reshape(sigmaBacArray,nSamplesMinor + 1,1);

figure
plot(tMinor,STATES(:,9)','DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
hold on
plot(tKF,ESTIMATES(:,9)','DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(tKF,ESTIMATES(:,9)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
plot(tKF,ESTIMATES(:,9)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
hold off
ylabel('biomass X_{bac} [g/L]')
title('Estimated and true biomass concentration')
legend()

%% Plot trajectories of relevant states: 

figure()

% S_IN:
subplot(3,2,1)
plot(tMinor,STATES(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,3),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('S_{IN} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% S_h2o:
subplot(3,2,2)
plot(tMinor,STATES(:,4),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,4),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0.8,1.1])
ylabel('S_{h2o} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_ch_fast:
subplot(3,2,3)
plot(tMinor,STATES(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,5),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,fast} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_ch_slow:
subplot(3,2,4)
plot(tMinor,STATES(:,6),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,6),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,slow} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_bac:
subplot(3,2,5)
plot(tMinor,STATES(:,9),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,9),'DisplayName','estimate',...
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
axis tight % draw axis only as long as time vectors

% X_ash:
subplot(3,2,6)
plot(tMinor,STATES(:,10),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,ESTIMATES(:,10),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0.015,0.025])
ylabel('X_{ash} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

sgtitle('Comparison of true and EKF-estimates of states')