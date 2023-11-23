%% Version
% (R2022b) Update 5
% Erstelldatum: 29.08.2023
% last modified: 23.11.2023
% Autor: Simon Hellmann

%% DAS EKF fürs ADM1-R4-frac-norm mit Online/Offline Messwerten gemäß Intensivbeprobung
% but without delay and augmentation

close all
clear all
clc

global counter
counter = 0; 

% Load Measurement Data:
load Messung_ADM1_R4_frac_norm_IntensiveSampling.mat

% get numeric values of normalization matrices: 
TxNum = TNum.Tx; % states
TyNum = TNum.Ty; % outputs
TuNum = TNum.Tu; % inputs

%% set up time vectors required for multirate measurements:
tMinor = MESS.tOnline;          % time grid of minor measurements 
nSamplesMinor = length(tMinor); % number of online measurements taken
% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMinor); 
dt = diff_t(1);                 % sampling time
tKF = [tMinor;tMinor(end)+dt];  % time vector for Kalman Filter

% get true time instances:
tStdSample = MESS.tStdSample;    

%% create united measurement array with NaNs where no (offline) measurement was returned
MEASOn = MESS.yMeasOn;
MEASStd = MESS.yMeasStd;
qOn = size(MEASOn,2);   % # secondary (online) measurements
qStd = size(MEASStd,2); % # primary (standard/offline) measurements
MEASUnite = nan(nSamplesMinor,qOn+qStd);     % allocate memory
% insert minor instances in first qOn columns: 
MEASUnite(:,1:qOn) = MEASOn; 
% insert minor instances (at right timing) in last qOff cols of MEASUnite:
[Lia,idxTStd] = ismember(tStdSample,tMinor); % determine correct row positioning of major instances in MEASUnite (right timing)
MEASUnite(idxTStd,qOn+1:qOn+qStd) = MEASStd; 

%% Initialization of EKF 
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

% initial state in normalized coordinates: 
x0InitNorm = MESS.x0Norm'; 
x0Norm = x0InitNorm; 
xHatNorm = xHat./TxNum; 
xMinusNorm = xHatNorm;      % to be overwritten

% set up raw matrices. all matrices which require initialization have 1 
% column more than the number of measurement instances:
STATES = MESS.xSolOn; 
ESTIMATES = zeros(nSamplesMinor + 1,nStates);
COVARIANCE = zeros(nStates,nStates,nSamplesMinor + 1);
% same for normalized coordinates: 
STATESNorm = MESS.xSolOnNorm; 
ESTIMATESNorm = zeros(nSamplesMinor + 1,nStates);
COVARIANCENorm = zeros(nStates,nStates,nSamplesMinor + 1); 

% Initialize Kalman Filter:
ESTIMATES(1,:) = xHat;
ESTIMATESNorm(1,:) = xHatNorm; 

%% tuning of EKF: 
P0 = diag((xHat-x0).^2);    % Schneider und Georgakis
PMinus = P0;                % to overwrite
COVARIANCE(:,:,1) = P0; 
% same for normalized coordinates: 
P0Norm = P0./(TxNum.^2);    % for comparison with non-normalized case
PMinusNorm = P0Norm;        % to overwrite
COVARIANCENorm(:,:,1) = P0Norm; 

% measurement uncertainty: 
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
% fine tuning: 
R(4,4) = 5E3*R(4,4);    % SIN
R(5,5) = 5E2*R(5,5);    % TS
RNorm = R./(TyNum.^2); 

% process uncertainty: 
Q = diag([0.016, 0.555, 0.563, 958.4, 1.263, 2.816, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
Q(4,4) = 1E-3*Q(4,4);       % unit change for h2o [g/l] --> [kg/l]
Q(10,10) = 1E-3*Q(10,10);   % unit change for ash [g/l] --> [kg/l]  
QNorm = Q./(TxNum.^2);

% obtain feeding information:
inputMat = MESS.inputMat;   % [tEvents,feed vol flow,inlet concentrations]
tEvents = inputMat(:,1);    % times of feeding events (on/off)

% if inputMat has first feeding entry only at t > t0, add zero row at t0:
if tEvents(1) > 0
    nColcInputMat = size(inputMat,2);   % number of columns in inputMat
    inputMat = [zeros(nColcInputMat,1),inputMat]; 
end

nIntervals = length(tEvents); % # constant feeding regimes

% normalize the different columns of inputMat appropriately:
inputMatNorm = inputMat;    % copy
inputMatNorm(:,2) = inputMatNorm(:,2)./TuNum; % normalized inputs
inputMatNorm(:,3:end) = inputMatNorm(:,3:end)./repmat(TxNum',nIntervals,1); % normalize inlet concentrations with state normalization

%% derive all function handles for EKF from symbolic model equations
% define symbolic variables ("S") (all vectors are column vectors):
xS = sym('x', [nStates,1]);         % states as col. vector
syms uS real                        % input
xiS = sym('xi', [nStates,1]);       % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));   % time-variant parameters (theta)
cS = sym('c', size(params.c));      % time-invariant parameters (known & constant)
aS = sym('a', size(params.a));      % petersen matrix with stoichiometric constants

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

%% derive all function handles in normalized coordinates from symbolilc expressions:
nOutputs = size(MESS.yClean,2); 
xNormS = sym('xNorm', [nStates 1]);  % normalized states as col. vector
syms uNorm real;                % normalized input
xiNormS = sym('xi', [nStates,1]);    % normalized inlet concentrations 
TxS = sym('Tx', [nStates,1]);        % normalization matrix for states
TyS = sym('Ty', [nOutputs,1]);         % normalization matrix for outputs
syms Tu real                    % normalization variable for input

dynamicsNorm = BMR4_AB_frac_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = BMR4_AB_frac_norm_mgl_sym(xNormS, cS, TxS, TyS); 

% partial derivatives for EKF (symbolic): 
dfdxNormSym = jacobian(dynamicsNorm,xNormS); 
dhdxNormSym = jacobian(outputsNorm,xNormS); 

% convert to numeric function handles:
fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 
dfdxNorm = matlabFunction(dfdxNormSym, 'Vars', {xNormS, uNorm, thS, cS, aS, TxS, Tu}); 
dhdxNorm = matlabFunction(dhdxNormSym, 'Vars', {xNormS, cS, TxS, TyS}); 

%% call EKF iteratively
tic
% integrate over all (online) measurement intervals (like in reality):
for k = 1:nSamplesMinor 
    % Mind: for now I work with parallel structure involving both absolute 
    % and normalized coordinates: 
    tSpan = [tKF(k) tKF(k+1)]; % measurement interval. In reality, this is the 
    tk = tSpan(1);      % t_k
    tkp1 = tSpan(2);    % t_k+1
    
    % get most recent measurement:
    yMeas = MEASUnite(k,:);    % simulated measurement
%     % check if you're at minor or major instance: 
%     flagMajor = all(~isnan(yMeas)); % 0: minor instance, 1: major instance
%     
%     if flagMajor == 1
%         disp('now major instance!')
%     end    

    %% get feeding information:
    % pass only relevant feedings during the measurement interval, because
    % only those would be known in reality:
    idxRelEvents = find(tEvents >= tk & tEvents <= tkp1);
    tRelEvents = tEvents(idxRelEvents); % Auswertung anhand Index
    
    % find the critcal last feeding event before current measurement interval:
    idxLastEvent = find(tEvents < tk,1,'last');

    % Case a: constant feeding during measurement interval:
    if isempty(tRelEvents) 
        feedInfo = inputMat(idxLastEvent,:);        % abs. coordinates
        feedInfoNorm = inputMatNorm(idxLastEvent,:);% rel. coordinates
    % Case b: feeding changes during measurement interval:
    else
        % use all relevant feeding events during measurement interval and 
        % the critical last one before:
        feedInfo = inputMat([idxLastEvent,idxRelEvents],:);         % abs. coordinates
        feedInfoNorm = inputMatNorm([idxLastEvent,idxRelEvents],:); % rel. coordinates
    end

    %% execute multirate EKF
    % absolute coordinates: 
    [xPlus,PPlus] = extendedKalmanFilterSlice(xMinus,PMinus, ...
        feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx,tSpan,nStates); 
    
    % normalized coordinates:  
    [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormSlice(xMinusNorm,PMinusNorm, ...
        feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm, ...
        TxNum,TyNum,TuNum,tSpan,nStates); 
    
    % save results, abs. coordinates:
    ESTIMATES(k+1,:) = xPlus(1:nStates)';   % only keep the "real" states, not the sample state in case of augmentation
    COVARIANCE(:,:,k+1) = PPlus(1:nStates,1:nStates);   % same for P-matrix
    % same for normalized coordinates:
    ESTIMATESNorm(k+1,:) = xPlusNorm(1:nStates); 
    COVARIANCENorm(:,:,k+1) = PPlusNorm(1:nStates,1:nStates);
    
    % Update for next iteration, abs. coordinates:  
    xMinus = xPlus;     % estimated state from Kalman Filter
    PMinus = PPlus;     % state error covariance matrix
    % same for normalized coordinates:
    xMinusNorm = xPlusNorm; 
    PMinusNorm = PPlusNorm; 

end
toc

%% compute system outputs (both absolute and normalized):
EKFOutput = nan(nSamplesMinor + 1,nOutputs);    % allocate memory
EKFOutputNorm = nan(nSamplesMinor + 1,nOutputs);% allocate memory
for k = 1:nSamplesMinor + 1 % also consider the initial value at t = t0
    EKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
    EKFOutputNorm(k,:) = gNorm(ESTIMATESNorm(k,:)',params.c,TxNum,TyNum);
end

%% de-normalize state estimates and outputs:
xEKFDeNorm = repmat(TxNum',nSamplesMinor+1,1).*ESTIMATESNorm;
yEKFDeNorm = repmat(TyNum',nSamplesMinor+1,1).*EKFOutputNorm;

% normalize state estimates formerly in abs. coordinates (for comparison): 
xEKFNorm = ESTIMATES./repmat(TxNum',nSamplesMinor+1,1);
% states and outputs computed through both ways (absolute/normalized 
% coordinates) are numerically identical (except for roundoff errors)

%% compute goodness of fit for all measurements
RMSSE = zeros(nOutputs,1);         % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:nOutputs
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = yEKFDeNorm(2:end,kk);  % base: normalized coordinates
%     estimatedMeasurements = EKFOutput(2:end,kk); % base: absolute coordinates
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

RMSSE_mean = mean(RMSSE); 

%% Plot results

% plot model output based on EKF estimates, compare with real measurements:
yCleanOn = MESS.yCleanOn;   % clean online model outputs
yCleanStd = MESS.yCleanStd; % clean offline model outputs without delay
feedVolFlow = inputMat(:,2);% [l/d]

colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figure

% gas volume flow: 
subplot(3,2,1)
scatter(tMinor, MEASOn(:,1)/24,'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMinor,yCleanOn(:,1)/24,'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(tKF,yEKFDeNorm(:,1)/24,'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0,15])
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
plot(tKF,yEKFDeNorm(:,2),'DisplayName','EKF-Output',...
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
plot(tKF,yEKFDeNorm(:,3),'DisplayName','EKF-Output',...
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
scatter(tStdSample,MEASStd(:,1),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample; tMinor(end)],[yCleanStd(:,1);yCleanStd(end,1)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,yEKFDeNorm(:,4),'DisplayName','EKF-Output',...
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
scatter(tStdSample,MEASStd(:,2),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample; tMinor(end)],[yCleanStd(:,2);yCleanStd(end,2)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,yEKFDeNorm(:,5),'DisplayName','EKF-Output',...
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
scatter(tStdSample,MEASStd(:,3),'DisplayName','noisy measurements',...
        'Marker','o', 'Color', colorPaletteHex(1)); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample; tMinor(end)],[yCleanStd(:,3);yCleanStd(end,3)],...
     'DisplayName','clean, un-delayed output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
plot(tKF,yEKFDeNorm(:,6),'DisplayName','EKF-Output',...
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
plot(tMinor,STATES(:,9),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
hold on
plot(tKF,xEKFDeNorm(:,9),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(tKF,xEKFDeNorm(:,9)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
plot(tKF,xEKFDeNorm(:,9)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
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
plot(tKF,xEKFDeNorm(:,3),'DisplayName','estimate',...
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
plot(tKF,xEKFDeNorm(:,4),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([600,1200])
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
plot(tKF,xEKFDeNorm(:,5),'DisplayName','estimate',...
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
plot(tKF,xEKFDeNorm(:,6),'DisplayName','estimate',...
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
plot(tKF,xEKFDeNorm(:,9),'DisplayName','estimate',...
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
plot(tKF,xEKFDeNorm(:,10),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.015,0.025])
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