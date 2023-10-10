%% Version
% (R2022b) Update 5
% Erstelldatum: 06.10.2023
% Autor: Simon Hellmann

% create the MESS-vector (synthetic measurement data) for Kalman Filtering. 
% Modell: ADM1-R4_Core (no water, no ash, no nitrogen, no gas phase)

addpath('modelEquations/');

clc
clear
close all

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch = 0.25;        % war vorher der Wert für kch
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 
% fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011;      
Kco2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0130;    % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [m3/bar/d] (actually [l/bar/d])
T = 273.15;     % operating temperature [K]
Vl = 100E-3;    % liquid volume, aus Sörens GitHub [m^3]
Vg = 10E-3;     % gas volume, aus Sörens GitHub [l]
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
aNum = [0.2482,  0.6809,   -1,     0,      0,       0.1372; 
        0.3221,  0.7954,   0,     -1,      0,       0.1723; 
        0.6393,  0.5817,   0,      0,     -1,       0.2286; 
        0,       0,        0.18,   0.77,   0.05,   -1]';

% inlet concentrations [GitHub Sören], vmtl. Rindergülle
%      S_ch4, S_IC,X_ch, X_pr,    X_li,  X_bac
xIn = [0,     0,   0.592, 23.398, 4.75,  1.381]'; % [kg/m3], 

% initial value: 
%       S_ch4, S_IC,  X_ch, X_pr,  X_li,  X_bac
x0SS = [0.091, 0.508, 3.26, 0.956, 0.413, 2.569]'; 

% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
thNum = [kch, 0, kpr, kli, kdec, 0]'; % [kch, (kchS), kpr, kli, kdec, (fracChFast)] 
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
feedMax = 10*24*1E-3;  % max. feed volume flow [l/h] converted to [m3/d]
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
xS = sym('x', size(xIn));       % states as col. vector
syms uS real                    % input
xiS = sym('xi', size(xIn));     % inlet concentrations (assumed known) 
thS = sym('th', size(thNum));   % 6 time-variant parameters (theta)
cS = sym('c', size(cNum));      % 21 known & constant time-invariant parameters 
aS = sym('a', size(aNum));      % petersen matrix with stoichiometric constants

dynamics = ADM1_R4_Core_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = ADM1_R4_Core_mgl_sym(xS); 

% transform into numeric function handle. Note that the independet
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS}); 

% define ODE function handle and integrate system dynamics: 
feedVolFlowSS = totalFeed/tEnd; 
tSpanSS = [0,tSS]; 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0Init = xSS(end,:)';   % start actual simulation in steady state
x0 = x0Init; 

%% Solve ODE via iterative solution of constant feeding regimes (on or off)
xSim = zeros(length(tOverall), nStates);% allocate memory
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
    odeFun = @(t,x) f(x,feedVolFlowCurr,xInCurr,thNum,cNum,aNum); 
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode           = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % in this case, the solver would interpret 
        % t_ode as a time span and choose integration time points on his own
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % request to evaluate at exactly 3 time points
        [tVec, solVec] = ode15s(odeFun, t_ode, x0);
        % of the 3 time points evaluated, only save 2 (first and last):
        xSim(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode has more than two time instants, which are the times 
        % when the integration is evaluated:
        [tVec, solVec] = ode15s(odeFun, t_ode, x0); % Returns >= 3 values
        xSim(idxTimeInterval,:) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0 = solVec(end, :);    % update initial value for next interval
end
toc

% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSol = xSim(idxGrid,:);

%% compute output variables
N = size(xSol,1);    % number of survived sampling points
q = 3;   
yClean = zeros(N,q); % allocate memory
for k = 1:N
    yClean(k,:) = g(xSol(k,:)')'; % Simons Implementierung (arXiv)
end

%% add noise to measurements acc to sensor data sheets
 
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaCh4 = 0.8; % assumed noise on liquid CH4 measurements [kg/m3] 
sigmaCo2 = 1; % assumed noise on liquid CO2 measurements [kg/m3]
sigmaBac = 0.4; % NH4-N [kg/m3]

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaCh4, sigmaCo2, sigmaBac];
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

% Sch4: 
subplot(3,1,1)
scatter(tGrid,yMeas(:,1),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tGrid,yClean(:,1),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('S_{ch4} in kg/m3')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Sco2:
subplot(3,1,2)
scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('S_{co2} in kg/m3')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Xbac:  
subplot(3,1,3)
scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,3),'r-.',...
     'LineWidth',1.2,'DisplayName','clean')
ylabel('biomass in kg/m3')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m3/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4-frac-noWater')

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0Init; 
MESS.x = xSol; 
MESS.xSim = [tSim,xSim]; 
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = yClean;  
MESS.yMeas = yMeas; 
MESS.C = noiseCovMat; % accurate values from sensor data sheets

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
fileName = 'Messung_ADM1_R4_Core.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params')
