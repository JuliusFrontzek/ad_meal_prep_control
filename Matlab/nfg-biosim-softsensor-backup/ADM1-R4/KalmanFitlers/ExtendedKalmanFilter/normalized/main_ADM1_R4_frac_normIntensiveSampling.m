%% Version
% (R2022b) Update 5
% Erstelldatum: 17.08.2023
% Autor: Simon Hellmann

% create synthetic measurement data for Kalman Filtering
% scenario: no delays, online + offline rates (Std.-Analyses)
% model: ADM1-R4-frac-norm (including ash and 2 carbohydrate fractions) 

clc
clear
% close all

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
k_ch_F = 0.25;        % war vorher der Wert für kch
k_ch_S = 1E-1*k_ch_F;   % selbst gewählt
k_pr = 0.2; 
k_li = 0.1; 
k_dec = 0.02; 
fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
K_H_ch4 = 0.0011;      
K_H_co2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
p_h2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p_atm = 1.0130;    % atmospheric pressure [bar]
k_La = 200;      % mass transfer coefficient [1/d]
k_p = 5e4;       % friction parameter [l/bar/d]
T = 273.15;     % operating temperature [K]
V_liq = 100;       % liquid volume, aus Sörens GitHub [l]
V_gas = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % mass density of digestate [g/l]
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
xIn = [0,     0,   0.592, 960.512, 23.398, 0,     4.75,  1.381, 0,    17,     0,    0]'; % [g/l], 
% xAshIn = 17 selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

% initial value: allocation of carbohydrates to fast/slow fraction acc. to
% fracChFast; assume 1 g/l ash concentration:
x0ch = 3.26;    % total carbohydrates initial value
x0SS = [0.091, 0.508, 0.944, 956.97, fracChFast*x0ch, (1-fracChFast)*x0ch, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 

% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
thNum = [k_ch_F, k_ch_S, k_pr, k_li, k_dec, fracChFast]'; % [kchF, kchS, kpr, kli, kdec, fracChFast] 
params.th = thNum; 

% times: 
tEnd = 14;  % [d] End of Simulation
tSS = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.25;    % [d]
intLong = 1.5;      % [d]
intMed = 0.5;       % [d]
% start with feeding right after steady state and some pause: 
% ints = [intLong,intShort,intLong,intMed]';
% ... and transform intervals to absolute times:
% cumInts = cumsum(ints);     % cumulated intervals
tFeedOnWeek = [1;2.5;4.5;6];  % beginning times of feedings (1 week)
tFeedOn = [tFeedOnWeek;tFeedOnWeek + 7]; 
feedingDurationsWeek = [intMed; intShort; intMed; intShort];
feedingDurations = repmat(feedingDurationsWeek,2,1); 
tFeedOff = tFeedOn + feedingDurations; % end times of feedings
tEvents = sort([0;tFeedOn;tFeedOff]); 
dt = 5/60/24;           % sample time [min], converted to [d]
tOnline = (0:dt:tEnd)'; % time grid. The model outputs will be evaluated here later on
tOverall = unique([tOnline; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactorsWeek = [70,30,50,40]'/100; 
feedFactors = repmat(feedFactorsWeek,2,1); 
portions = feedFactors*feedMax; % [l/d]         	
% steady state feed volume flow [l/d] should be the average of what is fed
% during dynamic operation:
totalFeed = sum(feedingDurations.*portions);  
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

% transform into numeric function handle. Note that the independet
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% define ODE function handle to determine steady state: 
feedVolFlowSS = totalFeed/tEnd; 
tSpanSS = [0,tSS]; 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
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
[tVecSSNormSym,xSSNormSym] = ode15s(odeFunNormSS,tSpanSS,x0SSNorm);
x0InitNorm  = xSSNormSym(end,:);
x0Norm      = x0InitNorm;

% compute normalized system output at steady state: 
ySSNorm     = gNorm(xSSNormSym(end,:)', cNum, TxNum, TyNum); 

% perform de-normalization to check if normalization works properly: 
xSSDeNorm   = repmat(TxNum',size(xSSNormSym,1),1).*xSSNormSym;
ySSDeNorm   = TyNum.*ySSNorm;

%% Solve ODE via iterative solution of constant feeding regimes (on or off)
xSimNorm = zeros(length(tOverall), nStates);% allocate memory
tSim = zeros(length(tOverall),1);           % allocate memory

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
    idxTimeInterval = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode           = tOverall(idxTimeInterval); 
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
idxGridOn = ismember(tOverall, tOnline); 
xSolOnNorm = xSimNorm(idxGridOn,:);

% de-normalize online evaluation of states: 
NOn = size(xSolOnNorm,1);    % number of survived online sampling points
xSolOn = repmat(TxNum',NOn,1).*xSolOnNorm;

%% compute output variables
% case 1: only online measurements
q = 6;   
yCleanNorm = zeros(NOn,q); % allocate memory
for k = 1:NOn
    yCleanNorm(k,:) = gNorm(xSolOnNorm(k,:)',cNum,TxNum,TyNum)'; % Simons Implementierung (arXiv)
end
% de-normalize outputs:
yClean = repmat(TyNum',NOn,1).*yCleanNorm;

% case 2: online and offline measurements (no delay)
dtStd = 3;          % sample time for offline measurements [d]
tStdOffset = 0.5;   % dont start first std. measurement right at beginning 
% times when standard samples were taken (TS, VS, NH4-N):
tStdSample = (tStdOffset:dtStd:tEnd)';  % no delay 
NStd = numel(tStdSample);   % # std. measurement samples

% interpolate xSimNorm at standard offline sample times: 
xSolStdNorm = interp1(tOverall,xSimNorm,tStdSample);
% de-normalize standard offline evaluation of states: 
xSolStd = repmat(TxNum',NStd,1).*xSolStdNorm;
 
qOn = 3;  % # online measurement signals
qOff = 3; % # offline measurement signals
yCleanOnNorm = zeros(NOn,qOn);      % allocate memory
yCleanStdNorm = zeros(NStd,qOff);
% evaluate model outputs at online sampling times:
for k = 1:NOn
    fullOutputNorm = gNorm(xSolOnNorm(k,:)',cNum,TxNum,TyNum)';
    yCleanOnNorm(k,:) = fullOutputNorm(1:qOn); 
end
% de-normalize:
yCleanOn = repmat(TyNum(1:qOn)',NOn,1).*yCleanOnNorm;

% ... then at standard offline sampling times:
for kk = 1:NStd
    fullOutputNorm = gNorm(xSolStdNorm(kk,:)',cNum,TxNum,TyNum)';
    yCleanStdNorm(kk,:) = fullOutputNorm(qOn+1:qOn+qOff); 
end
% de-normalize:
yCleanStd = repmat(TyNum(qOn+1:qOn+qOff)',NStd,1).*yCleanStdNorm;

%% add noise to measurements acc to sensor data sheets
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2*24;    % Trommelgaszähler Labor [L/h] --> [L/d]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.12;    % NH4-N [g/L]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
sigmasOn = [sigmaV, sigmaCh4, sigmaCo2]; 
sigmasOff = [sigmaSIN, sigmaTS, sigmaVS];

sigmaMat = repmat(sigmas,NOn,1);
sigmaMatOn = repmat(sigmasOn,NOn,1);
sigmaMatOff = repmat(sigmasOff,NStd,1);

% create normally distributed measurement noise matrices:
% zero mean for all measurements 
yMean = zeros(NOn,q); 
yMeanOn = zeros(NOn,qOn); 
yMeanStd = zeros(NStd,qOff); 

rng('default');     % fix seed for random number generation (for replicable results)
normalMeasNoise = normrnd(yMean,sigmaMat);
normalMeasNoiseOn = normrnd(yMeanOn,sigmaMatOn);
normalMeasNoiseOff = normrnd(yMeanStd,sigmaMatOff);

% add noise to clean model outputs:
yMeas = yClean + normalMeasNoise; 
yMeasOn = yCleanOn + normalMeasNoiseOn; 
yMeasStd = yCleanStd + normalMeasNoiseOff; 

% construct measurement noise covariance matrices:
noiseCovMat = diag(sigmas.^2);
noiseCovMatOn = diag(sigmasOn.^2);  
noiseCovMatOff = diag(sigmasOff.^2); 

%% Plot results (separated into Online/Offline Measurements)
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figure()

% online measurements:
% gas volume flow: 
subplot(3,2,1)
scatter(tOnline,yMeasOn(:,1)/24,'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yClean(:,1)/24,'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [l/h]')
ylim([0,18])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tOnline,yMeasOn(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tOnline,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} in bar')
ylim([0,0.85])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
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
ylim([0.2,0.6])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Offline measurements:
% SIN:  
subplot(3,2,4)
scatter(tStdSample,yMeasStd(:,1),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample;tEnd],[yCleanStd(:,1);yCleanStd(end,1)],...
       'DisplayName','clean','LineStyle','-.',...
       'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('inorg. nitrogen in g/L')
ylim([0.35,0.85])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% TS:  
subplot(3,2,5)
scatter(tStdSample,yMeasStd(:,2),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample;tEnd],[yCleanStd(:,2);yCleanStd(end,2)],...
    'DisplayName','clean','LineStyle','-.',...
    'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('total solids [-]')
xlabel('time [d]')
ylim([0.02,0.065])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color',colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(3,2,6)
scatter(tStdSample,yMeasStd(:,3),'DisplayName','noisy',...
        'Marker','o', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
% hold last measurement value till end of simulation:
stairs([tStdSample;tEnd],[yCleanStd(:,3);yCleanStd(end,3)],...
       'DisplayName','clean', 'LineStyle','-.', ...
       'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('volatile solids [-]')
xlabel('time [d]')
ylim([0.56,0.58])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
% legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4-frac-norm')

%% save results in struct MESS: 

MESS.tOnline = tOnline;
MESS.tStdSample = tStdSample; 
% MESS.tOfflineArrival = tOfflineArrival; 
MESS.x0 = x0Init; 
MESS.x0Norm = x0Norm; 
% state trajectories evaluated at diff. sampling times
MESS.xSolOn = xSolOn;   
MESS.xSolStd = xSolStd;
MESS.xSolOnNorm = xSolOnNorm;   % normalized state trajectories
MESS.xSolOffNorm = xSolStdNorm;
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
% clean outputs: 
MESS.yClean = yClean; 
MESS.yCleanOn = yCleanOn;  
MESS.yCleanOff = yCleanStd;
MESS.yCleanOnNorm = yCleanOnNorm;   % normalized outputs
MESS.yCleanOffNorm = yCleanStdNorm; 
% noisy outputs: 
MESS.yMeas = yMeas; 
MESS.yMeasOn = yMeasOn; 
MESS.yMeasOff = yMeasStd; 
MESS.C = noiseCovMat; % accurate values from sensor data sheets
MESS.COn = noiseCovMatOn;
MESS.COff = noiseCovMatOff;

save('Messung_ADM1_R4_frac_norm_MultiRate.mat', 'MESS', 'params','TNum')
