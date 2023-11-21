%% Version
% (R2022b) Update 6
% Erstelldatum: 25.07.2023
% last modified: 21.11.2023
% Autor: Simon Hellmann

% create the multirate synthetic measurement data for Kalman Filtering. 
% Modell: ADM1-R4-frac-norm (including ash and 2 carbohydrate fractions)

clc
clear
close all

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
k_ch_F = 0.25;        % war vorher der Wert für kch
k_ch_S = 1E-1*k_ch_F;   % selbst gewählt
k_pr = 0.2; 
k_li = 0.1; 
k_dec = 0.02; 
fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: [mol/l/bar] = [kmol/m³/bar] (Tab. B.7 in Sörens Diss)
K_H_ch4 = 0.0011;      
K_H_co2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
p_h2o = 0;      % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p_atm = 1.0130; % atmospheric pressure [bar]
k_La = 200;     % mass transfer coefficient [1/d]
k_p = 5e4;      % friction parameter [l/bar/d] --> upscaling assumption: remains same value in [m³/bar/d]
T = 273.15;     % operating temperature [K]
V_liq = 163;    % liquid volume FBGA [m³]
V_gas = 26;     % gas volume FBGA [m³]
rho = 1000;     % mass density of digestate [kg/m³]
Mch4 = 16;      % molar mass CH4 [kg/kmol]
Mco2 = 44;      % molar mass CO2 [kg/kmol]

%% order model parameters in the rights structures (prepare simulation)
c1 = 1/V_liq; 
c2 = k_La; 
c3 = k_La*K_H_ch4*R*T; 
c4 = k_La*K_H_co2*R*T; 
c5 = k_La*V_liq/V_gas; 
c6 = k_p/p_atm*(R*T/Mch4)^2;
c7 = 2*k_p/p_atm*(R*T)^2/Mch4/Mco2;
c8 = k_p/p_atm*(R*T/Mco2)^2;
c9 = k_p/p_atm*R*T/Mch4*(2*p_h2o - p_atm); 
c10 = k_p/p_atm*R*T/Mco2*(2*p_h2o - p_atm); 
c11 = k_p/p_atm*(p_h2o - p_atm)*p_h2o; 
c12 = R*T/Mch4;
c13 = R*T/Mco2;
c14 = rho; 
c15 = -k_p/p_atm/V_gas*(R*T/Mch4)^2;
c16 = -2*k_p/p_atm/V_gas*(R*T)^2/Mch4/Mco2;
c17 = -k_p/p_atm/V_gas*(R*T/Mco2)^2;
c18 = -k_p/p_atm/V_gas*(R*T/Mch4)*(2*p_h2o - p_atm);
c19 = -k_p/p_atm/V_gas*(R*T/Mco2)*(2*p_h2o - p_atm);
c20 = -k_La*V_liq/V_gas*K_H_ch4*R*T - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
c21 = -k_La*V_liq/V_gas*K_H_co2*R*T - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
c22 = V_liq/V_gas;
c23 = -k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
% combine all in column vector: 
cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23]';

% petersen matrix acc. to Weinrich 2021
% Note:all aij with |aij ~= 1|, only take the absolute value (see arXiv paper). 
% Index "Num" for Numeric
aNum = [0.2482,  0.6809,   0.0207,     0.0456,    -1,      0,      0,      0,       0.1372,    0,      0,          0; 
        0.2482,  0.6809,   0.0207,     0.0456,     0,     -1,      0,      0,       0.1372,    0,      0,          0;
        0.3221,  0.7954,   0.1689,     0.4588,     0,      0,     -1,      0,       0.1723,    0,      0,          0; 
        0.6393,  0.5817,   0.0344,     0.4152,     0,      0,      0,     -1,       0.2286,    0,      0,          0; 
        0,       0,        0,          0,          0.18,   0,      0.77,   0.05,   -1,         0,      0,          0;
       -1,       0,        0,          0,          0,      0,      0,      0,       0,         0,      c22,        0; 
        0,      -1,        0,          0,          0,      0,      0,      0,       0,         0,      0,          c22]';

% inlet concentrations [GitHub Sören], vmtl. Rindergülle
%      S_ch4, S_IC,S_IN,  S_h2o,   X_chF,  X_chS, X_pr,  X_li,  X_bac,X_ash,  S_ch4,g, S_co2,g
xIn = [0,     0,   0.592, 960.512, 23.398, 0,     4.75,  1.381, 0,    17,     0,    0]'; % [kg/m³], 
% xAshIn = 17 selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
thNum = [k_ch_F, k_ch_S, k_pr, k_li, k_dec, fracChFast]'; % [kchF, kchS, kpr, kli, kdec, fracChFast] 
params.th = thNum; 

%% create feeding pattern: 
intShort = 0.25/2;    % [d]
intMed = 0.5/2;       % [d]
intLong = 1.25/2;      % [d]
intPattern = [intLong;intMed;intLong;intLong;intShort;intShort;intLong];
% start with feeding right after steady state and some pause: 
tFeedOnWeek = cumsum([intPattern;intPattern]);  
feedingDuration = 20/60/24;     % [min] --> [d]
tFeedOn = tFeedOnWeek;% [tFeedOnWeek;tFeedOnWeek + 7]; 
tFeedOff = tFeedOn + feedingDuration; % end times of feedings

% set other times: 
tEnd = 7;   % [d] End of Simulation
tSS = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 
dt = 20/60/24;          % sample time [min], converted to [d]. 
tOnline = (0:dt:tEnd)'; % time grid for online measurements. highest frequency in this script
tEvents = sort([0;tFeedOn;tFeedOff]); 
tMinor = unique([tOnline; tEvents]);% Join and sort timestamps

% construct vector of feed volume flows at times tFeedOn (we start in
% steady state but wait with first feeding for a little)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 1212;  % max. feed volume flow [m³/d]; results in about Raumbelastung 4 on average
% feed factors of max. feed volume flow [%]:
facHigh = 100;  
facMed = 50; 
facLow = 25;
facPattern = [facHigh;facLow;facHigh;facHigh;facMed;facLow;facMed]; 
feedFactorsWeek = [facPattern;facPattern]./100; % 1 single week
feedFactors = repmat(feedFactorsWeek,1,1);  % mulitiples of 1 week
portions = feedFactors*feedMax; % [m³/d]         	
% steady state feed volume flow [m³/d] should be the average of what is fed
% during dynamic operation:
totalFeed = sum(feedingDuration.*portions);  
avFeedVolFlow = totalFeed/tEnd; 
feedVolFlow = zeros(nIntervals,1);  % allocate memory
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nStates = length(xIn);
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % allocate memory 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn',nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; 

%% derive odeFun from symbolic model definition and determine steady state
% define symbolic ("S") variables (all vector are defined as column vectors)
xS = sym('x', [12 1]);   % states as col. vector
syms uS real             % input
xiS = sym('xi', [12,1]); % inlet concentrations (assumed known) 
thS = sym('th', [6,1]);  % 6 time-variant parameters (theta)
cS = sym('c', [21,1]);   % 21 known & constant time-invariant parameters 
aS = sym('a', [12,7]);% petersen matrix with stoichiometric constants

dynamics = BMR4_AB_frac_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = BMR4_AB_frac_mgl_sym(xS,cS);
outputsExt = BMR4_AB_frac_mgl_gasVolFlows_sym(xS,cS); % extended output with volume flows of CH4 and CO2

% transform into numeric function handle. Note that the independet
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 
gExt = matlabFunction(outputsExt, 'Vars', {xS, cS}); 

% define ODE function handle to determine steady state: 
tSpanSS = [0,tSS]; 
feedVolFlowSS = avFeedVolFlow;
% feedVolFlowSS = 21.65; % hiermit bekommt man Raumbelastung 4
odeFunSS = @(t,x) f(x,avFeedVolFlow,xIn,thNum,cNum,aNum); 

% initial value: allocation of carbohydrates to fast/slow fraction acc. to
% fracChFast; assume 1 kg/m^3 ash concentration:
x0ch = 3.26;    % total carbohydrates initial value
x0SS = [0.091, 0.508, 0.944, 956.97, fracChFast*x0ch, (1-fracChFast)*x0ch, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 

[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0Init = xSS(end,:)';   % start actual simulation in steady state
x0 = x0Init;

%% derive normalized system equations
xNormS      = sym('xNorm', [12 1]);  % normalized states as col. vector
syms uNorm real;                     % normalized input
xiNormS     = sym('xi', [12,1]);     % normalized inlet concentrations 
TxS         = sym('Tx', [12,1]);     % normalization matrix for states
TyS         = sym('Ty', [ 6,1]);     % normalization matrix for outputs
syms Tu real                         % normalization variable for input

dynamicsNorm = BMR4_AB_frac_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = BMR4_AB_frac_norm_mgl_sym(xNormS, cS, TxS, TyS); 

% turn into numeric function handles: 
fNorm       = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm       = matlabFunction(outputsNorm,  'Vars', {xNormS, cS, TxS, TyS}); 

% define numeric values for normalization with steady state: 
TxNum       = x0;    % ensure no numeric rounoff errors lead to neg. concentrations!
y0          = g(x0,cNum);    % steady state output
TyNum       = y0; 
TuNum       = feedVolFlowSS; 
% summarize in stuct to save them: 
TNum        = struct; 
TNum.Tx     = TxNum; 
TNum.Ty     = TyNum; 
TNum.Tu     = TuNum;

% normalization of simulation inputs:
uNorm       = feedVolFlowSS./TuNum; 
xInNorm     = xIn./TxNum; 
x0SSNorm    = x0SS./TxNum; 

% simulate transition into steady state in normalized coordinates: 
odeFunNormSS = @(t,xNorm) fNorm(xNorm,uNorm,xInNorm,thNum,cNum,aNum,TxNum,TuNum); 
[tVecSSNorm,xSSNorm] = ode15s(odeFunNormSS,tSpanSS,x0SSNorm);
x0InitNorm  = xSSNorm(end,:);
x0Norm      = x0InitNorm;

% compute normalized system output at steady state: 
ySSNorm     = gNorm(xSSNorm(end,:)', cNum, TxNum, TyNum); 

% perform de-normalization to check if normalization works properly: 
xSSDeNorm   = repmat(TxNum',size(xSSNorm,1),1).*xSSNorm;
ySSDeNorm   = TyNum.*ySSNorm;

%% Dynamic simulation
% via iterative solution of constant feeding regimes (on or off)

xSimNorm = zeros(length(tMinor), nStates);% allocate memory
tSim = zeros(length(tMinor),1);           % allocate memory

% integrate ODEs for each interval (=time when feeding constantly =on or
% =off):
tic
for cI = 1:nIntervals
    if cI == nIntervals   % letztes Intervall geht von letztem Fütterungsimpuls (aus) bis Simulationsende ...
        tCurrent   = tEvents(end);
        tNext      = tEnd;
    else    % ... alle anderen Fütterungsintervalle:
        tCurrent   = tEvents(cI);
        tNext      = tEvents(cI + 1);
    end
    
    % Get current feeding volume flow and inlet concentrations:
    inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
    % split into current values of feedVolFlow and xIn: 
    feedVolFlowCurr = inputVector(1); 
    xInCurr = inputVector(2:end)'; 
    
    % apply normalization:
    uCurrNorm   = feedVolFlowCurr./TuNum;   % inputs
    xInCurrNorm = xInCurr./TxNum;           % normalization of inlet concentrations just like state normalization

    % normalized function handle: 
    odeFunNorm = @(t,xNorm) fNorm(xNorm,uCurrNorm,xInCurrNorm,thNum,cNum,aNum,TxNum,TuNum); 
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval = (tMinor >= tCurrent & tMinor <= tNext);
    t_ode           = tMinor(idxTimeInterval); 
    if length(t_ode) == 2   % in this case, the solver would interpret 
        % t_ode as a time span and choose integration time points on his own
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % request to evaluate at exactly 3 time points
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm);
        % of the 3 time points evaluated, only save 2 (first and last):
        xSimNorm(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode has more than two time instants, which are the times 
        % when the integration is evaluated:
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm); % Returns >= 3 values
        xSimNorm(idxTimeInterval,:) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0Norm = solVec(end, :);    % update initial value for next interval
end
toc

% Evaluate xSol only in tOnline, discard the rest
idxGridOn = ismember(tMinor, tOnline); 
xSolOnNorm = xSimNorm(idxGridOn,:);

% de-normalize online evaluation of states: 
NOn = size(xSolOnNorm,1);    % number of survived online sampling points
xSolOn = repmat(TxNum',NOn,1).*xSolOnNorm;

%% compute all-online output variables
% Note: all-online measurements (unrealistic)
q = 6;   
yCleanNorm = zeros(NOn,q);  % allocate memory
yCleanExt = nan(NOn,q + 2); % allocate memory for two additional outputs
for k = 1:NOn
    yCleanNorm(k,:) = gNorm(xSolOnNorm(k,:)',cNum,TxNum,TyNum)'; % Simons Implementierung (arXiv)
    yCleanExt(k,:) = gExt(xSolOn(k,:)',cNum);   % additionally with volume flows of CH4 and CO2
end
% de-normalize normalized outputs:
yClean = repmat(TyNum',NOn,1).*yCleanNorm;

%% separate measurements into online and delayed offline
dtOff = 1;      % sample time for offline measurements [d]
% times when samples were taken:
tOfflineSample = (0.45:dtOff:tEnd)';  % offset offline measurements from online ones
NOff = numel(tOfflineSample); % # offline sample points

% interpolate xSimNorm at offline sample times: 
xSolOffNorm = interp1(tMinor,xSimNorm,tOfflineSample);
% de-normalize offline evaluation of states: 
xSolOff = repmat(TxNum',NOff,1).*xSolOffNorm;
 
qOn = 3;  % # online measurement signals
qOff = 3; % # offline measurement signals
yCleanOnNorm = zeros(NOn,qOn);      % allocate memory
yCleanOffNorm = zeros(NOff,qOff);
% evaluate model outputs at online sampling times:
for k = 1:NOn
    fullOutputNorm = gNorm(xSolOnNorm(k,:)',cNum,TxNum,TyNum)';
    yCleanOnNorm(k,:) = fullOutputNorm(1:qOn); 
end
% de-normalize:
yCleanOn = repmat(TyNum(1:qOn)',NOn,1).*yCleanOnNorm;

% ... and at offline sampling times:
for kk = 1:NOff
    fullOutputNorm = gNorm(xSolOffNorm(kk,:)',cNum,TxNum,TyNum)';
    yCleanOffNorm(kk,:) = fullOutputNorm(qOn+1:qOn+qOff); 
end
% de-normalize:
yCleanOff = repmat(TyNum(qOn+1:qOn+qOff)',NOff,1).*yCleanOffNorm;

% add time delays for arrival times of samplings
% construct times when offline and atline measurements return from lab: 
delayMin = 0.5;     % [d]
delaySpread = 2/24; % [min], converted to [d] (manually chosen) 
rng('default');     % fix seed for random number generation (for replicable results)
delayOff = delayMin + rand(NOff,1)*delaySpread;
tOfflineArrival = tOfflineSample + delayOff;

% optional clipping of return times by tEnd: 
% idxCliptOff = tOfflineRet <= tEnd; 
% tOfflineRetClip = tOfflineRet(idxCliptOff); 
% % correct remaining number of elements
% NOffClip = numel(tOfflineRetClip); 
% % incorporate knowledge of sampling sequence and incoming sequence: 
% [tOfflineRetSorted,indexOffline] = sort(tOfflineArrival); 
% trueIndexOffline = 1:NOff; 
% indexDevOffline = indexOffline' - trueIndexOffline; 

%% add noise to measurements according to sensor data sheets
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
sigmaV = 0.08*24; % Trommelgaszähler FBGA [m³/h] -> [m³/d]
% sigmaV = 0.2*24;    % Trommelgaszähler Labor [L/h] --> [L/d]
sigmapCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmapCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaVCh4 = 0.98;   % [m³/d], derived from Trommelgaszähler FBGA, mean volume flows of 100m³/d and 50/50% CH4/CO2 content
sigmaVCo2 = 0.98;   % [m³/d], see above
sigmaSIN = 0.12;    % NH4-N [kg/m³]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmapCh4, sigmapCo2, sigmaSIN, sigmaTS, sigmaVS]; 
sigmasExt = [sigmaV, sigmapCh4, sigmapCo2, sigmaVCh4, sigmaVCo2, sigmaSIN, sigmaTS, sigmaVS]; 
sigmasOn = [sigmaV, sigmapCh4, sigmapCo2]; 
sigmasOff = [sigmaSIN, sigmaTS, sigmaVS];

sigmaMat = repmat(sigmas,NOn,1);
sigmaMatExt = repmat(sigmasExt,NOn,1);
sigmaMatOn = repmat(sigmasOn,NOn,1);
sigmaMatOff = repmat(sigmasOff,NOff,1);

% create normally distributed measurement noise matrices:
% zero mean for all measurements 
yMean = zeros(NOn,q); 
yMeanExt = zeros(NOn,q+2);
yMeanOn = zeros(NOn,qOn); 
yMeanOff = zeros(NOff,qOff); 

% rng('default');     % fix seed for random number generation (for replicable results)
normalMeasNoise = normrnd(yMean,sigmaMat);
normalMeasNoiseExt = normrnd(yMeanExt,sigmaMatExt);
normalMeasNoiseOn = normrnd(yMeanOn,sigmaMatOn);
normalMeasNoiseOff = normrnd(yMeanOff,sigmaMatOff);

% add noise to clean model outputs:
yMeas = yClean + normalMeasNoise; 
yMeasExt = yCleanExt + normalMeasNoiseExt; 
yMeasOn = yCleanOn + normalMeasNoiseOn; 
yMeasOff = yCleanOff + normalMeasNoiseOff; 

% construct measurement noise covariance matrices:
noiseCovMat = diag(sigmas.^2);
noiseCovMatOn = diag(sigmasOn.^2);  
noiseCovMatOff = diag(sigmasOff.^2); 

%% Plot results (separated into Online/Offline Measurements)
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figOutputs = figure; 

% online measurements:
% gas volume flow: 
subplot(3,2,1)
scatter(tOnline,yMeasOn(:,1),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yClean(:,1),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [m³/d]')
ylim([0,400])
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
legend('Location','SouthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tOnline,yMeasOn(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tOnline,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(3,2,3)
scatter(tOnline,yMeasOn(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yClean(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{co2} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Offline measurements:
% SIN:  
subplot(3,2,4)
scatter(tOfflineSample,yMeasOff(:,1),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSample;tEnd],[yCleanOff(:,1);yCleanOff(end,1)],...
       'DisplayName','clean','LineStyle','-.',...
       'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('inorg. nitrogen in kg/m^3')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% TS:  
subplot(3,2,5)
scatter(tOfflineSample,yMeasOff(:,2),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSample;tEnd],[yCleanOff(:,2);yCleanOff(end,2)],...
    'DisplayName','clean','LineStyle','-.',...
    'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('total solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color',colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(3,2,6)
scatter(tOfflineSample,yMeasOff(:,3),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
% hold last measurement value till end of simulation:
stairs([tOfflineSample;tEnd],[yCleanOff(:,3);yCleanOff(end,3)],...
       'DisplayName','clean', 'LineStyle','-.', ...
       'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('volatile solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m^3/d]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4-frac')
fontsize(figOutputs, 14,'points'); 
%% save results in struct MESS: 

MESS.tOnline = tOnline;
MESS.tOfflineSample = tOfflineSample; 
MESS.tOfflineArrival = tOfflineArrival; 
MESS.x0 = x0Init; 
MESS.x0Norm = x0Norm; 
% state trajectories evaluated at diff. sampling times
MESS.xSolOn = xSolOn;   
MESS.xSolOff = xSolOff;
MESS.xSolOnNorm = xSolOnNorm;   % normalized state trajectories
MESS.xSolOffNorm = xSolOffNorm;
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [m^3/d]
% clean outputs: 
MESS.yClean = yClean; 
MESS.yCleanExt = yCleanExt; % extended by volFlows of CH4/CO2
MESS.yCleanOn = yCleanOn;  
MESS.yCleanOff = yCleanOff;
MESS.yCleanOnNorm = yCleanOnNorm;   % normalized outputs
MESS.yCleanOffNorm = yCleanOffNorm; 
% noisy outputs: 
MESS.yMeas = yMeas; 
MESS.yMeasExt = yMeasExt;   % extended by volFlows of CH4/CO2
MESS.yMeasOn = yMeasOn; 
MESS.yMeasOff = yMeasOff; 
MESS.C = noiseCovMat;       % accurate values from sensor data sheets
MESS.COn = noiseCovMatOn;
MESS.COff = noiseCovMatOff;

% order fiels alphabetically: 
MESS = orderfields(MESS); 

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
fileName = 'Messung_ADM1_R4_frac_norm_MultiRate.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params', 'TNum')

% save feeding information as CSVs in output-folder:
targetPath = '\\dbfz-user.leipzig.dbfz.de\user$\shellmann\GIT\testFiles\plotting'; 
feedInfo = array2table([tEvents,feedVolFlow], 'VariableNames',{'tEvents','feedVolumeFlow'});
fileNameTab = 'feedInfoTab'; 
writetable(feedInfo,fullfile(targetPath,fileNameTab))
