%% Version
% (R2022b) Update 6
% Erstelldatum: 31.07.2023
% last modified: 21.11.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-frac-norm mit Online/Offline Messwerten
% but only one sample during delay period each.

close all
clear all
clc

global counter
counter = 0; 

% Load Measurement Data:
load generatedOutput/Messung_ADM1_R4_frac_norm_MultiRate.mat

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

% obtain true time instances:
tOfflineSample = MESS.tOfflineSample;   
tOfflineArrival = MESS.tOfflineArrival; 
% shift/pad sampling and arrival times of offline measurements into fine 
% time grid of tMinor: 
tOfflineSampleShift = interp1(tMinor,tMinor,tOfflineSample,'next'); 
tOfflineArrivalShiftPre = interp1(tMinor,tMinor,tOfflineArrival,'next'); % includes NaN where extrapolation would be necessary
idxKeepOfflineMeasurements = ~isnan(tOfflineArrivalShiftPre); % remove NaNs
tOfflineArrivalShift = tOfflineArrivalShiftPre(idxKeepOfflineMeasurements); % remove NaNs

%% create united measurement array with NaNs where no (offline) measurement available
MEASOn = MESS.yMeasOn;
MEASOffPre = MESS.yMeasOff; % contains measurements arriving AFTER the end of simulation...
MEASOff = MEASOffPre(idxKeepOfflineMeasurements,:); % ... so drop those
qOn = size(MEASOn,2);       % # secondary (online) measurements
qOff = size(MEASOffPre,2);  % # primary (online) measurements
MEASUnite = nan(nSamplesMinor,qOn+qOff);     % allocate memory
% insert minor instances in first qOn columns: 
MEASUnite(:,1:qOn) = MEASOn; 
% insert minor instances (at right timing) in last qOff cols of MEASUnite:
[Lia,idxTOfflineArrPre] = ismember(tOfflineArrivalShift,tMinor); % determine correct row positioning of major instances in MEASUnite (right timing)
MEASUnite(idxTOfflineArrPre,qOn+1:qOn+qOff) = MEASOff; 

%% Initialize EKF:
% initial state values:
x0Init = MESS.x0;       % intitial value from createMess-file 
x0 = x0Init;            % x0 will be replaced in every iteration later
% false initial estimate (add random number but ensure positive concentrations)
nStates = length(x0); 
nTheta = length(params.th); % number of time-variant parameters
rng('default');         % fix seed for random number generation (for replicable results)
xHat = x0Init.*abs(randn(nStates,1)); % enforce false estimate in beginning
xMinus = xHat;          % to be overwritten
% xMinus = x0Init;      % XY Rania

% initial estimate in normalized coordinates: 
xHatNorm = xHat./TxNum; 
xMinusNorm = xHatNorm;      % to be overwritten

% include parametric plant-model mismatch: 
paramsOfSyntheticMeasData = params.th; 
perturbedParams = paramsOfSyntheticMeasData.*(1 + abs(randn(size(params.th))));
perturbedParams(end) = 0.7;     % make sure fracChFast is not too far off
params.th = perturbedParams; 

% set up raw matrices. all matrices which require initialization have 1 
% column more than the number of measurement instances:
ESTIMATESNorm = zeros(nSamplesMinor + 1,nStates);
COVARIANCENorm = zeros(nStates,nStates,nSamplesMinor + 1); 

%% Tune EKF: 
P0 = diag((xHat-x0).^2);    % Schneider und Georgakis
P0Norm = P0./(TxNum.^2);    % for comparison with non-normalized case
PMinusNorm = P0Norm;        % to overwrite
COVARIANCENorm(:,:,1) = P0Norm; 

% measurement uncertainty: 
buffer = 1.5;   % conservative safety margin of 50% for measurement noise covariance
R = buffer * MESS.C; 
RNorm = R./(TyNum.^2);
% fine-tuning:  
RNorm(1:3,1:3) = RNorm(1:3,1:3)*1E3; % reduce shakyness of gas measurements
RNorm(5,5) = RNorm(5,5)*2;   % TS a bit noisy
RNorm(6,6) = RNorm(6,6)*1E-2;% VS not so noisy

% process uncertainty: 
Q = diag([0.016, 0.555, 0.563, 958.4, 1.263, 2.816, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
QNorm = eye(nStates); 

% obtain feeding information:
inputMat = MESS.inputMat;   % [tEvents,feed vol flow,inlet concentrations]
tEvents = inputMat(:,1);    % times of feeding events (on/off)

% if inputMat has first feeding entry only at t > t0, add zero row at t0:
if tEvents(1) > 0
    nColcInputMat = size(inputMat,2);   % number of columns in inputMat
    inputMat = [zeros(nColcInputMat,1);inputMat]; 
end

% normalize the different columns of inputMat appropriately:
inputMatNorm = inputMat;    % copy
inputMatNorm(:,2) = inputMatNorm(:,2)./TuNum; % normalized inputs
nIntervals = length(tEvents); % # constant feeding regimes
inputMatNorm(:,3:end) = inputMatNorm(:,3:end)./repmat(TxNum',nIntervals,1); % normalize inlet concentrations with state normalization

%% derive all function handles for EKF from symbolic model equations (abs. coordinates)
% define symbolic variables ("S") (all vectors are column vectors):
xS = sym('x', [nStates,1]);         % states as col. vector
% syms uS real                        % input
% xiS = sym('xi', [nStates,1]);       % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));   % time-variant parameters (theta)
cS = sym('c', size(params.c));      % time-invariant parameters (known & constant)
aS = sym('a', size(params.a));      % petersen matrix with stoichiometric constants

% get extended output function with CH4/CO2 volFlows as numeric function: 
outputsExt = BMR4_AB_frac_mgl_gasVolFlows_sym(xS,cS);
% transform into numeric function handle: 
gExt = matlabFunction(outputsExt, 'Vars', {xS, cS}); 

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
flagAugmented = 0;      % 0: non-augmented; 1: augmented
flagDelayPeriod = 0;    % 0: not in delay period; 1: in delay period

% integrate over all (online) measurement intervals (like in reality):
for k = 1:nSamplesMinor 

    tSpan = [tKF(k);tKF(k+1)]; % measurement interval. In reality, this is the 
    tk = tSpan(1);      % t_k
    tkp1 = tSpan(2);    % t_k+1
        
    % check if you're AT primary sampling (k==s): 
    if ismember(tk,tOfflineSampleShift)
        flagAugmented = 1; 
        flagDelayPeriod = 1; 
        disp('sampling'); 
        % augment and initialize state x & state err. cov. matrix P
        xAugMinusNorm = [xMinusNorm;xMinusNorm]; 
        PAugMinusNorm = repmat(PMinusNorm,2,2); 
        xMinusNorm = xAugMinusNorm; 
        PMinusNorm = PAugMinusNorm; 
    end
    
    if flagDelayPeriod == 1
        disp(['in delay period... at time ', num2str(tk)]);
    end

    % get most recent measurement:
    yMeas = MEASUnite(k,:);    % simulated measurement
    % check if you're at minor or major instance: 
    flagMajor = all(~isnan(yMeas)); % 0: minor instance, 1: major instance
    
    if flagMajor == 1
        disp(['now major instance!', newline,...
        'Delay period over.'])
        flagDelayPeriod = 0;
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
        feedInfoNorm = inputMatNorm(idxLastEvent,:);% rel. coordinates
    % Case b: feeding changes during measurement interval:
    else
        % use all relevant feeding events during measurement interval and 
        % the critical last one before:
        feedInfoNorm = inputMatNorm([idxLastEvent;idxRelEvents],:); % rel. coordinates
    end

    %% execute multirate EKF
    
    % normalized coordinates:  
    [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormMultiRateSingleDelay(xMinusNorm,PMinusNorm, ...
        feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm, ...
        TxNum,TyNum,TuNum,tSpan,nStates,qOn,qOff,flagAugmented,flagDelayPeriod,flagMajor); 
    
    % save results in normalized coordinates:
    ESTIMATESNorm(k+1,:) = xPlusNorm(1:nStates); 
    COVARIANCENorm(:,:,k+1) = PPlusNorm(1:nStates,1:nStates);
    
    % Update for next iteration:  
    xMinusNorm = xPlusNorm; 
    PMinusNorm = PPlusNorm; 

    % remove state augmentation once primary measurement (& delay period) is over:
    if ismember(tk,tOfflineArrivalShift)
        flagAugmented = 0; 
%         flagDelayPeriod = 0;
        disp(['primary measurement over', newline, ...
              'remove augmentation']); 
        % remove augmentation of state and state err. cov. matrix P:
        xMinusUnAugNorm = xMinusNorm(1:nStates); 
        PMinusUnAugNorm = PMinusNorm(1:nStates,1:nStates); 
        % overwrite old (augmented) enteties: 
        xMinusNorm = xMinusUnAugNorm; 
        PMinusNorm = PMinusUnAugNorm;  

    end

end
toc

%% compute system outputs (both absolute and normalized):
EKFOutputNorm = nan(nSamplesMinor + 1,nOutputs);% allocate memory
for k = 1:nSamplesMinor + 1 % also consider the initial value at t = t0
%     EKFOutput(k,:) = g(ESTIMATES(k,:)',params.c);
    EKFOutputNorm(k,:) = gNorm(ESTIMATESNorm(k,:)',params.c,TxNum,TyNum);
end

%% de-normalize state estimates and outputs:
xEKFDeNorm = repmat(TxNum',nSamplesMinor+1,1).*ESTIMATESNorm;
yEKFDeNorm = repmat(TyNum',nSamplesMinor+1,1).*EKFOutputNorm;

% additionally compute estimated gasVolFlows of CH4/CO2:
yEKFDeNormExt = nan(nSamplesMinor+1,qOn + 2 + qOff);  % allocate memory for two additional outputs
for k = 1:nSamplesMinor + 1
    yEKFDeNormExt(k,:) = gExt(xEKFDeNorm(k,:)',params.c);
end 

%% compute goodness of fit for all measurements: 
RMSSE = zeros(nOutputs,1);         % allocate memory

% compute RMSSE for each measurement signal:
for kk = 1:nOutputs
    measurements = MESS.yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = yEKFDeNorm(2:end,kk);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

RMSSE_mean = mean(RMSSE); 

%% compute filter efficiency for selected variables: 
% volFlowCH4:
[NSE_VCH4,eNSE_VCH4] = compute_NSE(MESS.yMeasExt(:,4), yEKFDeNormExt(2:end,4)); 

% TS:
% compute interpolated measurements first (without EKF, you would only have
% access to them): 
TSMeasInterp = interp1(tOfflineArrivalShift,MESS.yMeasOff(:,2),tMinor,'previous');
% reduce all vectors to time period for which there are measurements available: 
idxMeasAvailable = ~isnan(TSMeasInterp); 
TSMeasInterpAvail = TSMeasInterp(idxMeasAvailable); 
TSEKF = yEKFDeNormExt(2:end,7);  % leave out time t0
TSEKFAvail = TSEKF(idxMeasAvailable);
[NSE_TS,eNSE_TS] = compute_NSE(TSMeasInterpAvail,TSEKFAvail); 

%% Plot results

% plot model output based on EKF estimates, compare with real measurements:
yClean = MESS.yClean; 
yCleanOn = MESS.yCleanOn; 
yCleanExt = MESS.yCleanExt; 
yCleanOff = MESS.yCleanOff; % clean model outputs without delay (use with tOfflineSampleShift) 
yMeasExt = MESS.yMeasExt; 
feedVolFlow = inputMat(:,2);% [m^3/d]

colorPaletteHex = ["#fe9f6d","#de4968","#8c2981","#3b0f70","#000004"]; 
colorPaletteHexMagma = ["#fcfdbf","#fc8961","#b73779","#51127c","#000004"];

figOutputs = figure; 
% set text appearance to latex rendering: 
set(figOutputs,'defaultAxesTickLabelInterpreter','latex');  
set(figOutputs,'defaulttextinterpreter','latex');
set(figOutputs,'defaultLegendInterpreter','latex');

% gas volume flow: 
subplot(3,2,1)
scatter(tMinor, MEASOn(:,1),'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on; 
plot(tMinor,yCleanOn(:,1),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
plot(tKF,yEKFDeNorm(:,1),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0,15])
ylabel('gas vol flow [m^3/d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight  % draw axis only as long as time vectors
% xlabel('time [d]')
legend('Location','NorthEast'); 

% V_ch4: 
subplot(3,2,2)
scatter(tMinor, yMeasExt(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on
plot(tMinor,yCleanExt(:,4),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
plot(tKF,yEKFDeNormExt(:,4),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0,50])
ylabel("$\dot V_{ch4} \, \mathrm{[m^3\,d^{-1}]}$")
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis on

% % pch4: 
% subplot(3,2,2)
% scatter(tMinor, MEASOn(:,2),'DisplayName','noisy measurements',...
%         'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
% hold on
% plot(tMinor,yCleanOn(:,2),'DisplayName','clean output',...
%      'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
% plot(tKF,yEKFDeNorm(:,2),'DisplayName','EKF-Output',...
%      'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.4,0.85])
% ylabel('p_{ch4} [bar]')
% yyaxis right
% stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
% ylabel('feed vol flow [m^3/d]')
% set(gca, "YColor", 'k')     % make right y-axis black 
% % xlabel('time [d]')
% % legend('Location','NorthEast'); 
% axis tight % draw axis only as long as time vectors

% pco2:
subplot(3,2,3)
scatter(tMinor, MEASOn(:,3),'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on; 
plot(tMinor,yCleanOn(:,3),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5)
plot(tKF,yEKFDeNorm(:,3),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8)
ylim([0.2,0.6])
ylabel('p_{co2} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight % draw axis only as long as time vectors
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% SIN:  
subplot(3,2,4)
plot(tMinor,yClean(:,4),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5)
hold on; 
scatter(tOfflineSampleShift,MEASOffPre(:,1),'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
scatter(tOfflineArrivalShift,MEASOff(:,1),'DisplayName','noisy arrival',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
plot(tKF,yEKFDeNorm(:,4),'DisplayName','EKF-Output deNorm',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8)
ylabel('inorg. nitrogen [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
axis tight % draw axis only as long as time vectors
% xlabel('time [d]')
legend('Location','NorthEast'); 

% TS:  
subplot(3,2,5)
plot(tMinor,yClean(:,5)*100,'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5)
hold on; 
scatter(tOfflineSampleShift,MEASOffPre(:,2)*100,'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
scatter(tOfflineArrivalShift,MEASOff(:,2)*100,'DisplayName','noisy arrival',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
plot(tKF,yEKFDeNorm(:,5)*100,'DisplayName','EKF-Output deNorm',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8)
ylim([0,0.1]*100)
ylabel('total solids [%-FM]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
axis tight % draw axis only as long as time vectors
% legend('Location','NorthEast'); 

% VS:  
subplot(3,2,6)
plot(tMinor,yClean(:,6)*100,'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5)
hold on; 
scatter(tOfflineSampleShift,MEASOffPre(:,3)*100,'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
scatter(tOfflineArrivalShift,MEASOff(:,3)*100,'DisplayName','noisy arrival',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
plot(tKF,yEKFDeNorm(:,6)*100,'DisplayName','EKF-Output deNorm',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8)
ylim([0,1]*100)
ylabel('volatile solids [%-TS]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
axis tight % draw axis only as long as time vectors
% legend('Location','NorthEast'); 

sgtitle('Comparison of EKF and clean model output')
fontsize(figOutputs, 14, 'points')

%% Biomass X_bac:
STATES = MESS.xSolOn; 

covarianceBac = reshape(COVARIANCENorm(9,9,:),nSamplesMinor+1,1); 
covarianceBacDeNorm = covarianceBac*TxNum(9)^2;  
sigmaBac = sqrt(abs(covarianceBacDeNorm)); % exclude negative values because of numerical issues

figure
plot(tMinor,STATES(:,9),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5)
hold on
plot(tKF,xEKFDeNorm(:,9),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8)
plot(tKF,xEKFDeNorm(:,9)+sigmaBac,'--b','DisplayName','+1 \sigma boundary')
plot(tKF,xEKFDeNorm(:,9)-sigmaBac,':b','DisplayName','-1 \sigma boundary')
hold off
ylabel('biomass X_{bac} [kg/m^3]')
title('Estimated and true biomass concentration')
legend()

%% Plot trajectories of relevant states: 

figStates = figure; 

% S_IN:
subplot(3,2,1)
plot(tMinor,STATES(:,3),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,3),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('S_{IN} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% S_h2o:
subplot(3,2,2)
plot(tMinor,STATES(:,4),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,4),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
ylim([600,1200])
ylabel('S_{h2o} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
% xlabel('time [d]')
legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_ch_fast:
subplot(3,2,3)
plot(tMinor,STATES(:,5),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,5),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,fast} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_ch_slow:
subplot(3,2,4)
plot(tMinor,STATES(:,6),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,6),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{ch,slow} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_bac:
subplot(3,2,5)
plot(tMinor,STATES(:,9),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,9),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.4,0.7])
ylabel('X_{bac} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

% X_ash:
subplot(3,2,6)
plot(tMinor,STATES(:,10),'DisplayName','true',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
hold on; 
plot(tKF,xEKFDeNorm(:,10),'DisplayName','estimate',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0.015,0.025])
ylabel('X_{ash} [kg/m^3]')
yyaxis right
stairs(tEvents, feedVolFlow,'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
ylabel('feed vol flow [m^3/d]')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
% legend('Location','NorthEast'); 
axis tight % draw axis only as long as time vectors

sgtitle('Comparison of true and EKF-estimates of states')
fontsize(figStates, 14, 'points')

%% export data for plotting in matplotlib: 
targetPath = '\\dbfz-user.leipzig.dbfz.de\user$\shellmann\GIT\testFiles\plotting'; 

% tMinor und cleane/Messwerte yCleanExt/yMeasExt: 
trueValuesTabOn = array2table([tMinor,yCleanExt(:,4),yMeasExt(:,4),yClean(:,5)], 'VariableNames',{'tMinor','volFlowCH4Clean','volFlowCH4Meas','TSClean'});
trueValuesTabOff = array2table([tOfflineSampleShift, tOfflineArrivalShift, MEASOff(:,2)], 'VariableNames',{'tOffSampleShift','tOffArrShift','TSMeas'});
fileNameTrueValuesOn = 'trueOnlineValuesR4'; 
fileNameTrueValuesOff = 'trueOfflineValuesR4'; 
pathAndFileNameTrueOn = fullfile(targetPath,fileNameTrueValuesOn); 
pathAndFileNameTrueOff = fullfile(targetPath,fileNameTrueValuesOff); 
writetable(trueValuesTabOn,pathAndFileNameTrueOn)
writetable(trueValuesTabOff,pathAndFileNameTrueOff)

% tKF und Schätzwerte yEKFDeNormExt: 
estValuesTab = array2table([tKF,yEKFDeNormExt(:,4),yEKFDeNormExt(:,7)], 'VariableNames',{'tKF','volFlowCH4Est','TSEst'});
fileNameEstValues = 'estValuesR4'; 
pathAndFileNameEst = fullfile(targetPath,fileNameEstValues); 
writetable(estValuesTab,pathAndFileNameEst) 
