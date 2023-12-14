%% Version
% (R2022b) Update 6
% Erstelldatum: 23.11.2023
% last modified: 14.12.2023
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-frac-norm mit Online/Offline Messwerten
% with multiple samples during delay period

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
tMinor = MESS.tOnline(2:end);   % time grid of minor measurements. no measurement in initial state! 
nSamplesMinor = length(tMinor); % number of online measurements taken
% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMinor); 
dt = diff_t(1);               % sampling time
tKF = [tMinor(1)-dt;tMinor];  % time vector for Kalman Filter, starting at initial state

% obtain true time instances (ignore the index shift possibly caused by 
% excluding the first online measurement at t0):
tOfflineSample = MESS.tOfflineSample;   
tOfflineArrival = MESS.tOfflineArrival; 
% shift/pad sampling and arrival times of offline measurements into fine 
% time grid of tMinor: 
tOfflineSampleShift = interp1(tMinor,tMinor,tOfflineSample,'next'); 
tOfflineArrivalShiftPre = interp1(tMinor,tMinor,tOfflineArrival,'next'); % includes NaN where extrapolation would be necessary
idxKeepOfflineMeasurements = ~isnan(tOfflineArrivalShiftPre); % remove NaNs
tOfflineArrivalShift = tOfflineArrivalShiftPre(idxKeepOfflineMeasurements); % remove NaNs

%% create united measurement array with NaNs where no (offline) measurements available
MEASOn = MESS.yMeasOn(2:end,:); % also ignore measurement at initial value t0
MEASOffPre = MESS.yMeasOff; % contains measurements arriving AFTER the end of simulation...
MEASOff = MEASOffPre(idxKeepOfflineMeasurements,:); % ... so drop those
qOn = size(MEASOn,2);       % # secondary (online) measurements
qOff = size(MEASOffPre,2);  % # primary (online) measurements
MEASUnite = nan(nSamplesMinor,4+qOn+qOff);     % allocate memory [tMinor, tSamples, tArrivals, idxSample, yOn, yOff]
% insert minor sampling times and corresponding measurements: 
MEASUnite(:,1) = tKF(2:end); 
MEASUnite(:,4+1:4+qOn) = MEASOn; 
% insert major instances (at right timing) in last qOff cols of MEASUnite:
[~,idxTOfflineArrShift] = ismember(tOfflineArrivalShift,tMinor); % determine correct row positioning of major instances in MEASUnite (right timing)
MEASUnite(idxTOfflineArrShift,4+qOn+1:4+qOn+qOff) = MEASOff; 
% insert corresponding offline sample times: 
MEASUnite(idxTOfflineArrShift,2) = tOfflineSampleShift; 
% insert corresponding offline arrival times: 
MEASUnite(idxTOfflineArrShift,3) = tOfflineArrivalShift; 
% insert a 1 in 4th col whenever there are samples taken: 
[~,idxTOfflineSampleShift] = ismember(tOfflineSampleShift,tMinor); % determine correct row positioning of samples in MEASUnite (right timing)
MEASUnite(idxTOfflineSampleShift,4) = ones(numel(tOfflineSampleShift),1); 

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

%% call MR-EKF iteratively
tic
flagDelayPeriod = 0;    % 0: not in delay period; 1: in delay period
flagArrival = 0;        % 0: no offline measurement arrival, 1: arrival
nAug = 0;               % number of augmentations
tActiveSamples = [];    % allocate memory
tActiveArrivals = [];   % allocate memory

% integrate over all (online) measurement intervals (like in reality):
for k = 1:nSamplesMinor 
    
    % make index shift; note that tKF starts at t0 = 0, while the first 
    % measurement comes in at t1. The latter is what counts for the
    % measurement update within the Kalman Filter updates
    tkm1 = tKF(k); 
    tk = tKF(k+1); 
    tSpan = [tkm1;tk]; 

    knownMeasUnite = MEASUnite(1:k,:); 

    % check if you're AT a primary sampling (k==s). If so, remember this
    % moment as sample time:
    % XY: note that in theory, multiple samples could be drawn in one step
    if ismember(tk,tOfflineSampleShift)
        % assume that only one sample is drawn in every minor step!
        nAug = nAug + 1; % increase level of augmentation
        % add new knwon sampling time: 
        if isempty(tActiveSamples)
            tActiveSamples = tk;
        else
            tActiveSamples = unique([tActiveSamples;tk]);   
        end
        disp(['took new sample at time ',num2str(tk),'. Level of augmentation: ',num2str(nAug)]); 
        % augment and initialize state x & state err. cov. matrix P:
        xAugMinusNorm = [xMinusNorm;xMinusNorm]; 
        PAugMinusNorm = blkdiag(PMinusNorm,PMinusNorm); 
        xMinusNorm = xAugMinusNorm; 
        PMinusNorm = PAugMinusNorm; 
    end
    nActiveSamples = nnz(~isnan(tActiveSamples)); % number of non-zero elements
    
    if nAug > 0
        disp(['in delay period... at time ', num2str(tk)]); 
    end

    % from this knowledge, create matrix of active samples (only times matter): 
    if nActiveSamples > 0
    % out of knownSamplesMat, extract rows that belong to active samples:
        [~,idxActSamples] = ismember(tActiveSamples,knownMeasUnite(:,1)); 
        actSamplesMat = knownMeasUnite(idxActSamples,[1,3]); % [tSample, tArrival] 
    else 
        actSamplesMat = []; 
    end

    % find active arrival times: 
    if ismember(tk,tOfflineArrivalShift)
        disp('received offline measurement!'); 
        if isempty(tActiveArrivals)
            tActiveArrivals = tk;   
        else
            tActiveArrivals = unique([tActiveArrivals;tk]);   
        end
        flagArrival = 1; 
    end
    nActiveArrivals = nnz(~isnan(tActiveArrivals)); % number of non-zero elements
    
    % insert those active arrivals in measurement array: 
    if nActiveArrivals > 0 % otherwise leave actSamplesMat unchanged!
        [~,idxActArrivals] = ismember(tActiveArrivals,knownMeasUnite(:,3));
        % save corresponding entries for now ...: 
        activeArrivals = knownMeasUnite(idxActArrivals,[2,3]); % [tCorrespondingSample, tArrival] 
        % ... and add values at the corresponding rows of actSamplesMat: 
        [~,idxCorrectRowIn_actSampleMat] = ismember(activeArrivals(:,1),actSamplesMat(:,1)); 
        actSamplesMat(idxCorrectRowIn_actSampleMat,:) = activeArrivals; 
    end
    
    % retrieve most recent measurement (get anything that comes in):
    yMeas = MEASUnite(k,4+1:end);    % obtained full measurement
    
    %% get feeding information:
    % pass on only relevant feedings during the measurement interval, 
    % because only those would be known in reality:
    idxRelEvents = find(tEvents >= tkm1 & tEvents <= tk);
    tRelEvents = tEvents(idxRelEvents); % Auswertung anhand Index
    
    % find the critcal last feeding event before current measurement interval:
    idxLastEvent = find(tEvents < tkm1,1,'last');

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
    [xPlusNorm,PPlusNorm] = extendedKalmanFilterNormMultiRateMultiDelay(xMinusNorm,PMinusNorm, ...
        feedInfoNorm,yMeas,params,QNorm,RNorm,fNorm,gNorm,dfdxNorm,dhdxNorm, ...
        TxNum,TyNum,TuNum,tSpan,nStates,nAug,actSamplesMat); 
    disp('executed real EKF step.')

    % save results in normalized coordinates:
    ESTIMATESNorm(k+1,:) = xPlusNorm(1:nStates); 
    COVARIANCENorm(:,:,k+1) = PPlusNorm(1:nStates,1:nStates);
    
    % Update for next iteration:  
    xMinusNorm = xPlusNorm; 
    PMinusNorm = PPlusNorm; 

    % after primary measurement (& delay period), remove redeemed samples, 
    % corresponding arrivals and augmentation:
    % XY: note that in theory, measurements of multiple measurements could 
    % be returned in one step (resulting in multiple augmentations dropped)
    if flagArrival
        % assume that only one arrival is used every minor step!
        nAug = nAug - 1; 
        disp(['primary measurement over', newline, ...
              'remove 1 augmentation. Level of augmentation: ', num2str(nAug)]); 
        % remove entries from active samples ... : 
        [~,idxRemove] = ismember(tk,actSamplesMat(:,2)); 
        tSamplesToBeRemoved = actSamplesMat(idxRemove,1);
        [~,idxSampleRemove] = ismember(tSamplesToBeRemoved,tActiveSamples); 
        tActiveSamples = tActiveSamples(~idxSampleRemove); 
        
        % ... and remove entries from active arrivals: 
        tArrivalsToBeRemoved = actSamplesMat(idxRemove,2);
        [~,idxArrivalRemove] = ismember(tArrivalsToBeRemoved,tActiveArrivals);
        tActiveArrivals = tActiveArrivals(~idxArrivalRemove); 
        
        % remove augmentation of state and state err. cov. matrix P:
        xMinusUnAugNorm = xMinusNorm(1:nStates*(1 + nAug)); 
        PMinusUnAugNorm = PMinusNorm(1:nStates*(1 + nAug),1:nStates*(1 + nAug)); 
        % overwrite old (augmented) enteties: 
        xMinusNorm = xMinusUnAugNorm; 
        PMinusNorm = PMinusUnAugNorm;
        
        flagArrival = 0; 
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
    estimatedMeasurements = yEKFDeNorm(:,kk);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

RMSSE_mean = mean(RMSSE); 

% %% compute filter efficiency for selected variables: 
% % volFlowCH4:
% [NSE_VCH4,eNSE_VCH4] = compute_NSE(MESS.yMeasExt(:,4), yEKFDeNormExt(2:end,4)); 
% 
% % TS:
% % compute interpolated measurements first (without EKF, you would only have
% % access to them): 
% TSMeasInterp = interp1(tOfflineArrivalShift,MESS.yMeasOff(:,2),tMinor,'previous');
% % reduce all vectors to time period for which there are measurements available: 
% idxMeasAvailable = ~isnan(TSMeasInterp); 
% TSMeasInterpAvail = TSMeasInterp(idxMeasAvailable); 
% TSEKF = yEKFDeNormExt(2:end,7);  % leave out time t0
% TSEKFAvail = TSEKF(idxMeasAvailable);
% [NSE_TS,eNSE_TS] = compute_NSE(TSMeasInterpAvail,TSEKFAvail); 

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
% % set text appearance to latex rendering: 
% set(figOutputs,'defaultAxesTickLabelInterpreter','latex');  
% set(figOutputs,'defaulttextinterpreter','latex');
% set(figOutputs,'defaultLegendInterpreter','latex');

% gas volume flow: 
subplot(3,2,1)
scatter(tMinor, MEASOn(:,1),'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on; 
plot(tKF,yCleanOn(:,1),'DisplayName','clean output',...
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
scatter(tKF, yMeasExt(:,4),'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on
plot(tKF,yCleanExt(:,4),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
plot(tKF,yEKFDeNormExt(:,4),'DisplayName','EKF-Output',...
     'LineStyle','-', 'Color', colorPaletteHexMagma(3), 'LineWidth',0.8); 
% ylim([0,50])
ylabel("vol flow CH4 [m^3/d]")
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
plot(tKF,yCleanOn(:,3),'DisplayName','clean output',...
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
plot(tKF,yClean(:,4),'DisplayName','clean output',...
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
plot(tKF,yClean(:,5)*100,'DisplayName','clean output',...
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
plot(tKF,yClean(:,6)*100,'DisplayName','clean output',...
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
plot(tKF,STATES(:,9),'DisplayName','true',...
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
plot(tKF,STATES(:,3),'DisplayName','true',...
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
plot(tKF,STATES(:,4),'DisplayName','true',...
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
plot(tKF,STATES(:,5),'DisplayName','true',...
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
plot(tKF,STATES(:,6),'DisplayName','true',...
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
plot(tKF,STATES(:,9),'DisplayName','true',...
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
plot(tKF,STATES(:,10),'DisplayName','true',...
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

% %% export data for plotting in matplotlib: 
% targetPath = '\\dbfz-user.leipzig.dbfz.de\user$\shellmann\GIT\testFiles\plotting_/&blaDontOverwrite!'; 
% 
% % tMinor und cleane/Messwerte yCleanExt/yMeasExt: 
% trueValuesTabOn = array2table([tMinor,yCleanExt(:,4),yMeasExt(:,4),yClean(:,5)], 'VariableNames',{'tMinor','volFlowCH4Clean','volFlowCH4Meas','TSClean'});
% trueValuesTabOff = array2table([tOfflineSampleShift, tOfflineArrivalShift, MEASOff(:,2)], 'VariableNames',{'tOffSampleShift','tOffArrShift','TSMeas'});
% fileNameTrueValuesOn = 'trueOnlineValuesR4'; 
% fileNameTrueValuesOff = 'trueOfflineValuesR4'; 
% pathAndFileNameTrueOn = fullfile(targetPath,fileNameTrueValuesOn); 
% pathAndFileNameTrueOff = fullfile(targetPath,fileNameTrueValuesOff); 
% writetable(trueValuesTabOn,pathAndFileNameTrueOn)
% writetable(trueValuesTabOff,pathAndFileNameTrueOff)
% 
% % tKF und Schätzwerte yEKFDeNormExt: 
% estValuesTab = array2table([tKF,yEKFDeNormExt(:,4),yEKFDeNormExt(:,7)], 'VariableNames',{'tKF','volFlowCH4Est','TSEst'});
% fileNameEstValues = 'estValuesR4'; 
% pathAndFileNameEst = fullfile(targetPath,fileNameEstValues); 
% writetable(estValuesTab,pathAndFileNameEst) 
