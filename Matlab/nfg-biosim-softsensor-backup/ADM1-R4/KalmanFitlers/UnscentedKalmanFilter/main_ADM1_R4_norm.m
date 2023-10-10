%% Version
% (R2022b) Update 5
% Erstelldatum: 1.10.2023
% last modified: 07.10.2023
% Autor: Simon Hellmann

% create the MESS-vector (synthetic measurement data) for Kalman Filtering. 
% Modell: ADM1-R4 (mit Asche, nur 1 CH-Fraktion)

addpath('modelEquations/');

close all
clear
clc

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch    = 0.25;
kpr     = 0.2; 
kli     = 0.1; 
kdec    = 0.02; 
% fracChFast  = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: [mol/L/bar] (Tab. B.7 in Sörens Diss)
K_H_ch4 = 0.0011;      
K_H_co2 = 0.0250; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R       = 0.08315; % id. gas constant [bar L/mol/K]
p_h2o   = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p_atm   = 1.0130;  % atmospheric pressure [bar]
k_La    = 200;     % mass transfer coefficient [1/d]
k_p     = 5e4;     % friction parameter [L/bar/d]

T_op    = 273.15;  % operating temperature [K]
V_liq   = 100E-3;  % liquid volume, aus Sörens GitHub [m3]
V_gas   = 10E-3;   % gas volume, aus Sörens GitHub [m3]
rho_liq = 1000;    % mass density of digestate [kg/m3]
M_ch4   = 16;      % molar mass CH4 [kg/kmol]
M_co2   = 44;      % molar mass CO2 [kg/kmol]

%% combine model parameters for compact control notation: 
c       = NaN(23,1);
c( 1)   = 1/V_liq; 
c( 2)   = k_La; 
c( 3)   = k_La*K_H_ch4*R*T_op; 
c( 4)   = k_La*K_H_co2*R*T_op; 
c( 5)   = k_La*V_liq/V_gas; 
c( 6)   = k_p/p_atm*(R*T_op/M_ch4)^2;
c( 7)   = 2*k_p/p_atm*(R*T_op)^2/M_ch4/M_co2;
c( 8)   = k_p/p_atm*(R*T_op/M_co2)^2;
c( 9)   = k_p/p_atm*R*T_op/M_ch4*(2*p_h2o - p_atm); 
c(10)   = k_p/p_atm*R*T_op/M_co2*(2*p_h2o - p_atm); 
c(11)   = k_p/p_atm*(p_h2o - p_atm)*p_h2o; 
c(12)   = R*T_op/M_ch4;
c(13)   = R*T_op/M_co2;
c(14)   = rho_liq; 
c(15)   = -k_p/p_atm/V_gas*(R*T_op/M_ch4)^2;
c(16)   = -2*k_p/p_atm/V_gas*(R*T_op)^2/M_ch4/M_co2;
c(17)   = -k_p/p_atm/V_gas*(R*T_op/M_co2)^2;
c(18)   = -k_p/p_atm/V_gas*(R*T_op/M_ch4)*(2*p_h2o - p_atm);
c(19)   = -k_p/p_atm/V_gas*(R*T_op/M_co2)*(2*p_h2o - p_atm);
c(20)   = -k_La*V_liq/V_gas*K_H_ch4*R*T_op - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
c(21)   = -k_La*V_liq/V_gas*K_H_co2*R*T_op - k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
c(22)   = V_liq/V_gas;
c(23)   = -k_p/p_atm/V_gas*(p_h2o - p_atm)*p_h2o;
cNum    = c;

% petersen matrix acc. to Weinrich 2021
% Note: for the arXiv version, all aij with |aij ~= 1|, only take the
% absolute value. Index "Num" for Numeric
aNum = [0.2482, 0.6809,   0.0207,     0.0456,    -1,      0,      0,       0.1372,    0,      0,          0; 
        0.3221, 0.7954,   0.1689,     0.4588,     0,     -1,      0,       0.1723,    0,      0,          0; 
        0.6393, 0.5817,   0.0344,     0.4152,     0,      0,     -1,       0.2286,    0,      0,          0; 
        0,      0,        0,          0,          0.18,   0.77,   0.05,   -1,         0,      0,          0;
       -1,      0,        0,          0,          0,      0,      0,       0,         0,      c(22),      0; 
        0,     -1,        0,          0,          0,      0,      0,       0,         0,      0,          c(22)]';

% inlet concentrations [GitHub Sören], vmtl. Rindergülle
%          S_ch4, S_IC,S_IN,  S_h2o,   X_ch,   X_pr,  X_li,  X_bac, X_ash,  S_ch4,g, S_co2,g
xIn     = [0,     0,   0.592, 960.512, 23.398, 4.75,  1.381, 0,     17,     0,       0]'; % [kg/m3]
% xAshIn = 17 selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

% combine constant parameters in struct: 
params      = struct;    % allocate memory
params.a    = aNum; 
params.c    = cNum; 
% fester Parametersatz (Hydrolysekonstanten) (index "Num" for numeric values)
thNum       = [kch, 0, kpr, kli, kdec, 0]'; % [kch, (kchS), kpr, kli, kdec, (fracChFast)] 
% füge den zum struct hinzu: 
params.th   = thNum; 

% times: 
tEnd        = 7;    % [d] End of Simulation
tSS         = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% design feeding pattern by defining interval lengths: 
intShort    = 0.5;    % [d]
intLong     = 2.5;    % [d]
intMed      = 1;      % [d]
% start with feeding right after steady state and some pause: 
ints        = [intLong,intShort,intLong,intMed]';  % feeding intervals
% ... and transform intervals to absolute times:
cumInts     = cumsum(ints);             % cumulated intervals
tFeedOn     = [cumInts(1);cumInts(3)];  % beginning times of feedings
tFeedOff    = [cumInts(2);cumInts(4)];  % end times of feedings
feedingDurations = [intShort; intMed];  
tEvents     = sort([0;tFeedOn;tFeedOff]); 
dt          = 0.5/24;                   % sample time [h], converte to [d]
tGrid       = (0:dt:tEnd)';             % time grid. The model outputs will be evaluated here later on
tOverall    = unique([tGrid; tEvents]); % Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals  = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax     = 10*24*1E-3;  % max. feed volume flow [L/h] converted to [m3/d]
feedFactors = [70,30]'/100; 
portions    = feedFactors*feedMax; % [L/d]         	
% steady state feed volume flow [L/d] should be the average of what is fed
% during dynamic operation:
totalFeed   = sum(feedingDurations.*portions);  
feedVolFlow = zeros(nIntervals,1);  % placeholder
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at time instances tFeedOn: 
nStates     = length(xIn);
nFeedings   = length(tFeedOn); 
xInMat      = zeros(nIntervals,nStates);         % placeholder 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn',nFeedings,1); 
inputMat    = [feedVolFlow,xInMat]; % structure required by right-hand side of ODE: u + xIn

%% derive system equations from symbolic model files (non-normalized): 
% define symbolic ("S") variables (all vector are defined as column vectors)
% syms t real              % dummy variable for time
xS          = sym('x', size(xIn));      % states as col. vector
syms uS real                            % input
xiS         = sym('xi', size(xIn));     % inlet concentrations (assumed known) 
thS         = sym('th', size(thNum));   % 6 time-variant parameters (theta), use only 4 of them here
cS          = sym('c',  size(cNum));    % 21 known & constant time-invariant parameters 
aS          = sym('a',  size(aNum));    % petersen matrix with stoichiometric constants

dynamics    = BMR4_AB_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs     = BMR4_AB_mgl_sym(xS,cS); % therein: vol flow in L/d

% transform into numeric function handles. Note that the independentxNorm
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f           = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g           = matlabFunction(outputs,  'Vars', {xS, cS}); 

%% determine steady state as initial value for simulation: 
tSpanSS     = [0,tSS]; 
feedVolFlowSS = totalFeed/tEnd;    
% inputVectorSS = [feedVolFlowSS,xIn']; 

% define function handle to determine the steady state:
odeFunSS    = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 

% Quelle: Sörens GitHub. Selbst angepasst: xAsh0, XCH je zu 50% auf XCHFast
% x0SS = [0.091, 0.508, 0.944, 956.97, 3.26, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 
x0SS        = [0.091, 0.508, 0.944, 956.97, 3.26, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0SS); 

x0Init      = xSS(end,:)';  % start dynamic simulation in steady state
x0          = abs(x0Init);  % ensure no numeric rounoff errors lead to neg. concentrations!

%% derive normalized system equations
xNormS      = sym('xNorm', size(xIn));  % normalized states as col. vector
syms uNorm real;                        % normalized input
xiNormS     = sym('xi', size(xIn));     % normalized inlet concentrations 
TxS         = sym('Tx', size(xIn));     % normalization matrix for states
TyS         = sym('Ty', size(thNum));   % normalization matrix for outputs
syms Tu real                            % normalization variable for input

dynamicsNorm = BMR4_AB_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = BMR4_AB_norm_mgl_sym(xNormS, cS, TxS, TyS); 

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

%% Solve ODE in a loop: iterative solution between all feeding points (both on & off)
xSimNorm    = zeros(length(tOverall), nStates);    % allocate memory
tSim        = zeros(length(tOverall),1);           % allocate memory

% integriere die System-DGLs abschnittsweise (jeder Bereich mit 
% Fütterung =on oder =off ist ein Abschnitt):
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
    xInCurr     = inputVector(2:end)';
    
    % apply normalization:
    uCurrNorm   = feedVolFlowCurr./TuNum; % inputs
    xInCurrNorm = xInCurr./TxNum;       % normalization of inlet concentrations just like state normalization

    % normalized function handle: 
    odeFunNorm = @(t,xNorm) fNorm(xNorm,uCurrNorm,xInCurrNorm,thNum,cNum,aNum,TxNum,TuNum); 

    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode           = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % der Solver würde t_ode hier als Zeitspanne 
        % interpretieren und seine eigenen Integrationszeitpunkte wählen.
        % request to evaluate at exactly 3 time points:
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % lege 3 äquidistante Zeitpunkte zur Auswertung der Integration fest
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm);
        % of the 3 time points evaluated, only save 2 (first and last):
        xSimNorm(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode hat mehr als zwei Zeitpunkte. Das sind die Integrations-Zeitpunkte für ode15s
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm); % Returns >= 3 values
        xSimNorm(idxTimeInterval,:) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0Norm  = solVec(end, :);    % update initial value for next interval
end
toc

% Evaluate xSim only in tGrid, discard the rest
idxGrid     = ismember(tOverall, tGrid); 
xSolNorm    = xSimNorm(idxGrid,:);

%% compute output variables
N           = size(xSolNorm,1);    % number of survived sampling points
q           = 6;   
yCleanNorm  = zeros(N,q); % allocate memory
for k = 1:N
    yCleanNorm(k,:) = gNorm(xSolNorm(k,:)', cNum, TxNum, TyNum); 
end

%% de-normalize states and outputs
xSol        = repmat(TxNum',N,1).*xSolNorm;
yClean      = repmat(TyNum',N,1).*yCleanNorm;

%% add noise to measurements acc to noise covariances: 
volFlowClean = yClean(:,1);
pCh4Clean   = yClean(:,2); 
pCo2Clean   = yClean(:,3);
SINClean    = yClean(:,4); 
TSClean     = yClean(:,5); 
VSClean     = yClean(:,6); 
 
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000;    % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV      = 0.2*24*1E-3;   % Trommelgaszähler Labor [L/h] -> [m3/d]
sigmaCh4    = 0.2/100;  % [Vol-%] -> [-]; für ideale Gase und p_atm ungefähr 1 bar: -> [bar]
sigmaCo2    = 0.2/100;  % [Vol-%] -> [-]; für ideale Gase und p_atm ungefähr 1 bar: -> [bar]
sigmaSIN    = 0.12;     % NH4-N [kg/m3]
sigmaTS     = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS     = 0.31/100; % [%] -> [-]
% combine all in sigma matrix and covariance matrix:
sigmas      = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
sigmaMat    = repmat(sigmas,N,1);
noiseCovMat = diag(sigmas.^2); 

nMeas       = length(pCh4Clean);  % number of measurements
rng('default');     % fix seed for random number generation (for replicable results)

% create normally distributed measurement noise matrix:
yMean       = zeros(N,q);      % zero mean for all online measurements 
normalMeasNoise = normrnd(yMean,sigmaMat);
yMeas       = yClean + normalMeasNoise; 

volFlowMeas = volFlowClean  + sigmaV*randn(nMeas,1);
pCh4Meas    = pCh4Clean     + sigmaCh4*randn(nMeas,1); 
pCo2Meas    = pCo2Clean     + sigmaCo2*randn(nMeas,1); 
SINMeas     = SINClean      + sigmaSIN*randn(nMeas,1);
TSMeas      = TSClean       + sigmaTS*randn(nMeas,1);
VSMeas      = VSClean       + sigmaVS*randn(nMeas,1);

%% Plot results (partial pressures, gas volume flow, SIN and feedings)
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figure()

% gas volume flow: 
subplot(3,2,1)
scatter(tGrid,yMeas(:,1),'.','DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,1),'DisplayName','noisy',...
%         'Marker','.', 'MarkerEdgeColor', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,1),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [m3/d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tGrid,yMeas(:,2),'.', 'DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tGrid,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(3,2,3)
scatter(tGrid,yMeas(:,3),'.', 'DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{co2} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% SIN:  
subplot(3,2,4)
scatter(tGrid,yMeas(:,4),'.', 'DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,4),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,4),'r',...
     'LineWidth',1.2,'DisplayName','clean')
ylabel('inorg. nitrogen in kg/m3')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% TS:  
subplot(3,2,5)
scatter(tGrid,yMeas(:,5),'.', 'DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,5),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,5),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('total solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(3,2,6)
scatter(tGrid,yMeas(:,6),'.', 'DisplayName','noisy', 'MarkerEdgeColor',colorPaletteHex(1), 'LineWidth',1);
% scatter(tGrid,yMeas(:,6),'DisplayName','noisy',...
%         'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
plot(tGrid,yClean(:,6),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('volatile solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4-frac-norm')

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0Init; 
MESS.x0Norm = x0InitNorm;   % normalized initial state
MESS.x = xSol; 
MESS.xNorm = xSolNorm;      % normalized state trajectories
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = yClean;
MESS.yCleanNorm = yCleanNorm; % normalized outputs
MESS.yMeas = yMeas; 
MESS.R = noiseCovMat; 

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
fileName = 'Messung_ADM1_R4_norm.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params', 'TNum')
