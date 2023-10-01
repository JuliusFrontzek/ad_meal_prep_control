%% Version
% (R2022b) Update 5
% Erstelldatum: 19.09.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-frac-noWater in normierten Koordinaten

close all
clear all
clc

addpath('modelEquations/');

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX
counterSigmaInit = 0;
counterSigmaProp = 0;
counterSigmaX = 0; 
counterX = 0; 

% Load Measurement Data:
load Messung_ADM1_R4_frac_noWater_norm

% get numeric values of normalization matrices: 
TxNum = TNum.Tx; % states
TyNum = TNum.Ty; % outputs
TuNum = TNum.Tu; % inputs

%% Initialization and Tuning of UKF 
% initial state values:
x0Init = MESS.x0;       % intitial value from createMess-file 
x0 = x0Init;            % x0 will be replaced in every iteration later
% false initial estimate (add random number but ensure positive concentrations)
nStates = length(x0); 
nTheta = length(params.th); % number of time-variant parameters
rng('default');     % fix seed for random number generation (for replicable results)
% xHat = x0.*(1 + 0.1*abs(randn(nStates,1))); 
xHat = x0.*abs(randn(nStates,1)); 
xMinus = xHat;      % to be overwritten
% xMinus = x0Init;    % XY Rania

% initial state in normalized coordinates: 
x0InitNorm = MESS.x0Norm'; 
x0Norm = x0InitNorm; 
xHatNorm = xHat./TxNum; 
xMinusNorm = xHatNorm;      % to be overwritten

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

% same for normalized coordinates: 
STATESNorm = MESS.xNorm; 
ESTIMATESNorm = zeros(nSamples + 1,nStates); 
COVARIANCENorm = zeros(nStates,nStates,nSamples + 1); 

% Initialize Kalman Filter:
ESTIMATES(1,:) = xHat;
P0 = eye(nStates); % XY: sicher besseres Tuning möglich
COVARIANCE(:,:,1) = P0; 
PMinus = P0;      % to overwrite
% same for normalized coordinates: 
ESTIMATESNorm(1,:) = xHatNorm;
P0Norm = eye(nStates);
COVARIANCENorm(:,:,1) = P0Norm; 
PMinusNorm = P0Norm;  % to overwrite

buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.R; 
% fine-tuning of Kalman Filter - measurement uncertainty: 
% R(4,4) = 5E1*R(4,4);    % SIN
% R(5,5) = 5E2*R(5,5);    % TS
RNorm = R./(TyNum.^2); 
% fine-tuning of Kalman Filter - process uncertainty: 
Q = diag([0.016, 0.555, 0.563, 1.263, 2.816, 2.654, 0.972, 2.894, 0.374, 0.948]);
% Q = eye(nStates); 
% Q(4,4) = 1E-3*Q(4,4);       % unit change for h2o [g/l] --> [kg/l]
% Q(10,10) = 1E-3*Q(10,10);   % unit change for ash [g/l] --> [kg/l]  
% QNorm = eye(nStates);   % XY: sicher besser auszulegen
QNorm = Q./(TxNum.^2); % for comparison with non-normalized case

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
xS = sym('x', [nStates,1]);     % states as col. vector
syms uS real                    % input
xiS = sym('xi', [nStates,1]);   % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));   % time-variant parameters (theta)
cS = sym('c', size(params.c));      % time-invariant parameters (known & constant)
aS = sym('a', size(params.a));      % petersen matrix with stoichiometric constants

% obtain symbolic objects of model equations:
dynamics = BMR4_AB_frac_noWater_ode_sym(xS, uS, xiS, thS, cS, aS); 
outputs = BMR4_AB_frac_noWater_mgl_sym(xS,cS);

% derive function handles of ODEs and output equation: 
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% % partial derivatives for EKF (symbolic): 
% dfdxSym = jacobian(dynamics,xS); 
% dhdxSym = jacobian(outputs,xS); 
% 
% % convert to numeric function handles: 
% dfdxNew = matlabFunction(dfdxSym, 'Vars', {xS, uS, thS, cS, aS}); 
% dhdxNew = matlabFunction(dhdxSym, 'Vars', {xS, cS}); 

%% derive all function handles in normalized coordinates from symbolilc expressions:
xNormS = sym('xNorm', [10 1]);  % normalized states as col. vector
syms uNorm real;                % normalized input
xiNormS = sym('xi', [10,1]);    % normalized inlet concentrations 
TxS = sym('Tx', [10,1]);        % normalization matrix for states
TyS = sym('Ty', [6,1]);         % normalization matrix for outputs
syms Tu real                    % normalization variable for input

dynamicsNorm = BMR4_AB_frac_noWater_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = BMR4_AB_frac_noWater_norm_mgl_sym(xNormS, cS, TxS, TyS); 

% convert to numeric function handles:
fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 

%% call UKF iteratively
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
        feedInfoNorm = inputMatNorm(idxLastEvent,:);
    % Case b: feeding changes during measurement interval:
    else
        % use all relevant feeding events during measurement interval and 
        % the critical last one before:
        feedInfo = inputMat([idxLastEvent,idxRelEvents],:);
        feedInfoNorm = inputMatNorm([idxLastEvent,idxRelEvents],:);
    end
    %% execute UKF

%     [xPlus,PPlus,Kv] = extendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdxNew,dhdxNew); % standard EKF
    [xPlus,PPlus] = unscKalmanFilterKolasAdditive(xMinus,PMinus,...
                            tSpan,feedInfo,yMeas,params,Q,R,f,g);
    [xPlusNorm,PPlusNorm] = unscKalmanFilterKolasAdditiveNorm(xMinusNorm,PMinusNorm,...
                        tSpan,feedInfoNorm,yMeas,params,QNorm,RNorm,...
                        fNorm,gNorm,TxNum,TyNum,TuNum); 
    
    % save results:
    ESTIMATES(k+1,:) = xPlus';
    COVARIANCE(:,:,k+1) = PPlus; 
%     GAIN(k,:) = Kv;     % Kalman Gain * Innovation
    % same for normalized coordinates:
    ESTIMATESNorm(k+1,:) = xPlusNorm'; 
    COVARIANCENorm(:,:,k+1) = PPlusNorm; 
    
    % Update for next iteration:  
    xMinus = xPlus;     % estimated state from Kalman Filter
    PMinus = PPlus;     % state error covariance matrix
    % same for normalized coordinates:
    xMinusNorm = xPlusNorm; 
    PMinusNorm = PPlusNorm; 

end
toc

%% compute system outputs from estimated states:
q = size(MESS.yMeas,2);     % number of measurement signals
UKFOutput = nan(nSamples + 1,q);   % allocate memory
UKFOutputNorm = nan(nSamples + 1,q);    % allocate memory
for k = 1:nSamples + 1
    UKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
    UKFOutputNorm(k,:) = gNorm(ESTIMATESNorm(k,:)',params.c,TxNum,TyNum);
end

yClean = MESS.yClean; 
feedVolFlow = inputMat(:,2);    % [l/d]

%% de-normalize estimates of states and outputs:

xUKFDeNorm = repmat(TxNum',nSamples+1,1).*ESTIMATESNorm;
yUKFDeNorm = repmat(TyNum',nSamples+1,1).*UKFOutputNorm;

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
scatter(tMeas, MESS.yMeas(:,1)/24,'DisplayName','noisy measurements',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tMeas,yClean(:,1)/24,'DisplayName','clean model output',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
plot(t,yUKFDeNorm(:,1)/24,'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);  
% plot(t,UKFOutput(:,1)/24,'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylim([0,15])
ylabel('gas vol flow [l/h]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
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
plot(t,yUKFDeNorm(:,2),'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);  
% plot(t,UKFOutput(:,2),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
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
plot(t,yUKFDeNorm(:,3),'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
% plot(t,UKFOutput(:,3),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
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
plot(t,yUKFDeNorm(:,4),'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
% plot(t,UKFOutput(:,4),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('inorg. nitrogen [g/L]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% % TS:  
% subplot(3,2,5)
% scatter(tMeas, MESS.yMeas(:,5),'DisplayName','noisy measurements',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% hold on; 
% plot(tMeas,yClean(:,5),'DisplayName','clean model output',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% plot(t,UKFOutput(:,5),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
% ylim([0,0.1])
% ylabel('total solids [-]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% ylabel('feed vol flow [l/h]')
% set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% % legend('Location','NorthEast'); 
% 
% % VS:  
% subplot(3,2,6)
% scatter(tMeas, MESS.yMeas(:,6),'DisplayName','noisy measurements',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% hold on; 
% plot(tMeas,yClean(:,6),'DisplayName','clean model output',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% plot(t,UKFOutput(:,6),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
% ylim([0,1])
% ylabel('volatile solids [-]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% ylabel('feed vol flow [l/h]')
% set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% % legend('Location','NorthEast'); 

sgtitle('Comparison of UKF and clean model output')

%% Biomass X_bac:
% sigmaBacArray = sqrt(COVARIANCE(8,8,:)); 
% sigmaBac = reshape(sigmaBacArray,nSamples + 1,1);
% 
% figure
% plot(tMeas,STATES(:,8)','DisplayName','true',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
% hold on
% plot(t,ESTIMATES(:,8)','DisplayName','estimate',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
% plot(t,ESTIMATES(:,8)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
% plot(t,ESTIMATES(:,8)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
% hold off
% ylabel('biomass X_{bac} [g/L]')
% title('Estimated and true biomass concentration')
% legend()

%% Plot trajectories of relevant states: 
trueStates = MESS.x; 

figure()

% S_IN:
subplot(2,2,1)
plot(tMeas,trueStates(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,3),'DisplayName','estimate',...
     'LineStyle',':', 'Color', colorPaletteHex(4), 'LineWidth',1.2); 
plot(t,xUKFDeNorm(:,3),'DisplayName','estimate deNorm',...
     'LineStyle','-.', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('S_{IN} [g/L]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 

% % S_h2o:
% subplot(3,2,2)
% plot(tMeas,trueStates(:,4),'DisplayName','true',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
% hold on; 
% plot(t,ESTIMATES(:,4),'DisplayName','estimate',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.8,1.1])
% ylabel('S_{h2o} [g/L]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% ylabel('feed vol flow [l/h]')
% set(gca, "YColor", 'k')     % make right y-axis black 
% % xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_ch_fast:
subplot(2,2,2)
plot(tMeas,trueStates(:,4),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,4),'DisplayName','estimate',...
     'LineStyle',':', 'Color', colorPaletteHex(4), 'LineWidth',1.2); 
plot(t,xUKFDeNorm(:,4),'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,fast} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% X_ch_slow:
subplot(2,2,3)
plot(tMeas,trueStates(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,5),'DisplayName','estimate',...
     'LineStyle','-.', 'Color', colorPaletteHex(4), 'LineWidth',1.2); 
plot(t,xUKFDeNorm(:,5),'DisplayName','estimate deNorm',...
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

% X_bac:
subplot(2,2,4)
plot(tMeas,trueStates(:,8),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
hold on; 
plot(t,ESTIMATES(:,8),'DisplayName','estimate',...
     'LineStyle','-.', 'Color', colorPaletteHex(4), 'LineWidth',1.2); 
plot(t,xUKFDeNorm(:,8),'DisplayName','estimate deNorm',...
     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); % ylim([0.4,0.7])
ylabel('X_{bac} [g/L]')
yyaxis right
stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
ylabel('feed vol flow [l/h]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
% 
% % S_ch4_gas:
% subplot(3,2,6)
% plot(tMeas,trueStates(:,5),'DisplayName','true',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
% hold on; 
% plot(t,ESTIMATES(:,5),'DisplayName','estimate',...
%      'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
% % ylim([0.4,0.7])
% ylabel('S_{ch4}^{gas} [g/L]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24,'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% ylabel('feed vol flow [l/h]')
% set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% % legend('Location','NorthEast'); 

sgtitle('Comparison of true and UKF-estimates of states')

%% Kalman Gains of EKF
% figure
% subplot(321)
% stairs(t(2:end),GAIN(1,:),'b')
% ylabel('S_{ch4} [g/L]')
% grid on
% %
% subplot(323)
% stairs(t(2:end),GAIN(3,:),'r')
% ylabel('S_{IN} [g/L]')
% grid on
% %
% subplot(325)
% stairs(t(2:end),GAIN(9,:),'g')
% ylabel('X_{ash} [g/L]')
% xlabel('Zeit [h]')
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%
% subplot(322)
% stairs(t(2:end),GAIN(5,:),'b')
% title('Korrektur der Zustände: K*v')
% ylabel('X_{ch} [g/L]')
% grid on
% %
% subplot(324)
% stairs(t(2:end),GAIN(6,:),'r')
% ylabel('X_{pr} [g/L]')
% grid on
% %
% subplot(326)
% stairs(t(2:end),GAIN(7,:),'g')
% ylabel('X_{li} [g/L]')
% xlabel('Zeit [h]')
% grid on
% sgtitle('Korrektur der Zustände: K*v')