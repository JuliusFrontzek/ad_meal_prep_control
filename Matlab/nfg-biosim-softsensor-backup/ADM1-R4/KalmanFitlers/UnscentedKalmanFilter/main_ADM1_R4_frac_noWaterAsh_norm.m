%% Version
% (R2022b) Update 5
% Erstelldatum: 19.09.2023
% Autor: Simon Hellmann

% create the MESS-vector (synthetic measurement data) for Kalman Filtering. 
% Modell: ADM1-R4-frac_noWater (no water, no ash, but 2 carbohydrate fractions)
% use model in normalized coordinates.

addpath('modelEquations/')

clc
clear
close all

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kchF = 0.25;        % war vorher der Wert für kch
kchS = 1E-1*kchF;   % selbst gewählt
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 
fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011;      
Kco2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0130;    % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [l/bar/d]
T = 273.15;     % operating temperature [K]
Vl = 100;       % liquid volume, aus Sörens GitHub [l]
Vg = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % mass density of digestate [kg/l]
Mch4 = 16;      % molar mass CH4 [kg/kmol]
Mco2 = 44;      % molar mass CO2 [kg/kmol]

%% order model parameters in the rights structures (prepare simulation)
c1 = 1/Vl; 
c2 = kla; 
c3 = kla*Kch4*R*T; 
c4 = kla*Kco2*R*T; 
c5 = kla*Vl/Vg; 
c6 = kp/p0*(R*T/Mch4)^2;
c7 = 2*kp/p0*(R*T)^2/Mch4/Mco2;
c8 = kp/p0*(R*T/Mco2)^2;
c9 = kp/p0*R*T/Mch4*(2*ph2o - p0); 
c10 = kp/p0*R*T/Mco2*(2*ph2o - p0); 
c11 = kp/p0*(ph2o - p0)*ph2o; 
c12 = R*T/Mch4;
c13 = R*T/Mco2;
c14 = rho; 
c15 = -kp/p0/Vg*(R*T/Mch4)^2;
c16 = -2*kp/p0/Vg*(R*T)^2/Mch4/Mco2;
c17 = -kp/p0/Vg*(R*T/Mco2)^2;
c18 = -kp/p0/Vg*(R*T/Mch4)*(2*ph2o - p0);
c19 = -kp/p0/Vg*(R*T/Mco2)*(2*ph2o - p0);
c20 = -kla*Vl/Vg*Kch4*R*T - kp/p0/Vg*(ph2o - p0)*ph2o;
c21 = -kla*Vl/Vg*Kco2*R*T - kp/p0/Vg*(ph2o - p0)*ph2o;
c22 = Vl/Vg;
c23 = -kp/p0/Vg*(ph2o - p0)*ph2o;
% combine all in column vector: 
cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23]';

% petersen matrix acc. to Weinrich 2021
% Note:all aij with |aij ~= 1|, only take the absolute value (see arXiv paper). 
% Index "Num" for Numeric
aNum = [0.2482,  0.6809,   0.0207,     -1,     0,      0,      0,       0.1372,    0,          0; 
        0.2482,  0.6809,   0.0207,     0,     -1,      0,      0,       0.1372,    0,          0;
        0.3221,  0.7954,   0.1689,     0,      0,     -1,      0,       0.1723,    0,          0; 
        0.6393,  0.5817,   0.0344,     0,      0,      0,     -1,       0.2286,    0,          0; 
        0,       0,        0,          0.18,   0,      0.77,   0.05,   -1,         0,          0;
       -1,       0,        0,          0,      0,      0,      0,       0,         c22,        0; 
        0,      -1,        0,          0,      0,      0,      0,       0,         0,          c22]';

% inlet concentrations [GitHub Sören], vmtl. Rindergülle
%      S_ch4, S_IC,S_IN,  X_chF,  X_chS, X_pr,  X_li,  X_bac,S_ch4,g, S_co2,g
xIn = [0,     0,   0.592, 23.398, 0,     4.75,  1.381, 0,    0,    0]'; % [g/l], 
% xAshIn = 17 selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

% initial value: 50/50 allocation of carbohydrates to fast/slow fraction; 
% assume 1 g/l ash concentration:
x0ch = 3.26;    % total carbohydrates initial value
x0SS = [0.091, 0.508, 0.944, fracChFast*x0ch, (1-fracChFast)*x0ch, 0.956, 0.413, 2.569,0.315, 0.78]'; 

% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
thNum = [kchF, kchS, kpr, kli, kdec, fracChFast]'; % [kchF, kchS, kpr, kli, kdec, fracChFast] 
params.th = thNum; 

% times: 
tEnd = 7;   % [d] End of Simulation
tSS = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.5;     % [d]
intLong = 2.5;      % [d]
intMed = 1;        % [d]
% start with feeding right after steady state and some pause: 
ints = [intLong,intShort,intLong,intMed]';
% ... and transform intervals to absolute times:
cumInts = cumsum(ints);     % cumulated intervals
tFeedOn = [cumInts(1);cumInts(3)];  % beginning times of feedings
tFeedOff = [cumInts(2);cumInts(4)]; % end times of feedings
feedingDurations = [intShort; intMed];  
tEvents = sort([0;tFeedOn;tFeedOff]); 
dt = 0.5/24;              % sample time [h], converte to [d]
tGrid = (0:dt:tEnd)';     % time grid. The model outputs will be evaluated here later on
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactors = [70,30]'/100; 
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
xS = sym('x', [10 1]);   % states as col. vector
syms uS real             % input
xiS = sym('xi', [10,1]); % inlet concentrations (assumed known) 
thS = sym('th', [6,1]);  % 6 time-variant parameters (theta)
cS = sym('c', [21,1]);   % 21 known & constant time-invariant parameters 
aS = sym('a', [10,7]);% petersen matrix with stoichiometric constants

dynamics = BMR4_AB_frac_noWaterAsh_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = BMR4_AB_frac_noWaterAsh_mgl_sym(xS,cS); 

% transform into numeric function handle. Note that the independet
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

% define ODE function handle and integrate system dynamics: 
feedVolFlowSS = totalFeed/tEnd; 
tSpanSS = [0,tSS]; 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0Init = xSS(end,:)';   % start actual simulation in steady state
x0 = x0Init; 

%% derive normalized system equations
xNormS      = sym('xNorm', [10 1]);  % normalized states as col. vector
syms uNorm real;                     % normalized input
xiNormS     = sym('xi', [10,1]);     % normalized inlet concentrations 
TxS         = sym('Tx', [10,1]);     % normalization matrix for states
TyS         = sym('Ty', [ 6,1]);     % normalization matrix for outputs
syms Tu real                         % normalization variable for input

dynamicsNorm = BMR4_AB_frac_noWaterAsh_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = BMR4_AB_frac_noWaterAsh_norm_mgl_sym(xNormS, cS, TxS, TyS); 

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
x0SSNorm    = x0SS./TxNum; 
xInNorm     = xIn./TxNum; 

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
tSim = zeros(length(tOverall),1);       % allocate memory

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
    uCurrNorm   = feedVolFlowCurr./TuNum; % inputs
    xInCurrNorm = xInCurr./TxNum;       % normalization of inlet concentrations just like state normalization

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

% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSolNorm = xSimNorm(idxGrid,:);

%% compute output variables
N = size(xSolNorm,1);    % number of survived sampling points
q = 4;   
yCleanNorm = zeros(N,q); % allocate memory
for k = 1:N
    yCleanNorm(k,:) = gNorm(xSolNorm(k,:)',cNum, TxNum, TyNum)'; % Simons Implementierung (arXiv)
end

%% de-normalize states and outputs
xSol = repmat(TxNum',N,1).*xSolNorm;
yClean = repmat(TyNum',N,1).*yCleanNorm;

%% add noise to measurements acc to sensor data sheets
 
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2*24;    % Trommelgaszähler Labor [L/h] --> [L/d]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.12;    % NH4-N [g/L]
% sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
% sigmaVS = 0.31/100; % [%] -> [-]

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN];%, sigmaTS, sigmaVS]; 
sigmaMat = repmat(sigmas,N,1);
noiseCovMat = diag(sigmas.^2);  % measurement noise covariance matrix

% create normally distributed measurement noise matrix:
yMean = zeros(N,q); % zero mean for all online measurements 
rng('default');     % fix seed for random number generation (for replicable results)
normalMeasNoise = normrnd(yMean,sigmaMat);
yMeas = yClean + normalMeasNoise; 

%% Plot results (partial pressures, gas volume flow, SIN and feedings)
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figure()

% gas volume flow: 
subplot(2,2,1)
scatter(tGrid,yMeas(:,1)/24,'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,1)/24,'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [l/h]')
ylim([0,16])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 

% pch4: 
subplot(2,2,2)
scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tGrid,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(2,2,3)
scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{co2} in bar')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% SIN:  
subplot(2,2,4)
scatter(tGrid,yMeas(:,4),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,4),'r',...
     'LineWidth',1.2,'DisplayName','clean')
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% % TS:  
% subplot(3,2,5)
% scatter(tGrid,yMeas(:,5),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
% hold on; 
% plot(tGrid,yClean(:,5),'DisplayName','clean',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
% ylabel('total solids [-]')
% xlabel('time [d]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% set(gca, "YColor", 'k')     % make right y-axis black 
% ylabel('feed vol flow [l/h]')
% legend('Location','NorthEast'); 
% set(gca, "YColor", 'k')
% 
% % VS:  
% subplot(3,2,6)
% scatter(tGrid,yMeas(:,6),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
% hold on; 
% plot(tGrid,yClean(:,6),'DisplayName','clean',...
%      'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
% ylabel('volatile solids [-]')
% xlabel('time [d]')
% yyaxis right
% stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
%        'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
% set(gca, "YColor", 'k')     % make right y-axis black 
% ylabel('feed vol flow [l/h]')
% legend('Location','NorthEast'); 
% set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4-frac-noWater')

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0Init; 
MESS.x = xSol;
MESS.xNorm = xSolNorm;      % normalized state trajectories
MESS.x0Norm = x0InitNorm;   % normalized initial state
% MESS.xSim = [tSim,xSimNorm]; 
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = yClean;
MESS.yCleanNorm = yCleanNorm;  
MESS.yMeas = yMeas; 
MESS.R = noiseCovMat; % accurate values from sensor data sheets

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
fileName = 'Messung_ADM1_R4_frac_noWaterAsh_norm.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params', 'TNum')