%% Version
% (R2022b) Update 5
% Erstelldatum: 05.10.2023
% Autor: Simon Hellmann

% create the MESS-vector (synthetic measurement data) for Kalman Filtering. 
% Modell: ADM1-R4 (including ash) using m³ and kg as units

addpath('modelEquations/');

close all
clear
clc

%% define numeric values of parameters
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch = 0.25;        % war vorher der Wert für kch
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 
% fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

% Henry coefficients: (Tab. B.7 in Sörens Diss) 
Kch4 = 0.0011; % [mol/l/bar]  
Kco2 = 0.025;  % [mol/l/bar]  

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar L/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0130;    % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e1;       % friction parameter [m³/bar/d]
T = 273.15;     % operating temperature [K]
Vl = 100e-3;    % liquid volume, aus Sörens GitHub. unit change [l] -> [m³]
Vg = 10e-3;     % gas volume, aus Sörens GitHub.  unit change [l] -> [m³]
rho = 1e3;      % mass density of digestate. unit change [kg/l] -> [kg/m³] 
Mch4 = 16;      % molar mass CH4 [g/mol] 
Mco2 = 44;      % molar mass CO2 [g/mol]

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
aNum = [0.2482, 0.6809,   0.0207,     0.0456,    -1,      0,      0,       0.1372,    0,      0,          0; 
        0.3221, 0.7954,   0.1689,     0.4588,     0,     -1,      0,       0.1723,    0,      0,          0; 
        0.6393, 0.5817,   0.0344,     0.4152,     0,      0,     -1,       0.2286,    0,      0,          0; 
        0,      0,        0,          0,          0.18,   0.77,   0.05,   -1,         0,      0,          0;
       -1,      0,        0,          0,          0,      0,      0,       0,         0,      c22,        0; 
        0,     -1,        0,          0,          0,      0,      0,       0,         0,      0,          c22]';

% inlet concentrations [GitHub Sören], vmtl. Rindergülle
%      S_ch4, S_IC,S_IN,  S_h2o,   X_chF,  X_pr,  X_li,  X_bac,X_ash,  S_ch4,g, S_co2,g
xIn = [0,     0,   0.592, 960.512, 23.398, 4.75,  1.381, 0,    17,     0,    0]'; % [g/l], 
% xAshIn = 17 selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)

% initial value:
% assume 1 g/l ash concentration:
x0SS = [0.091, 0.508, 0.944, 956.97, 3.26, 0.956, 0.413, 2.569, 1, 0.315, 0.78]'; 

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
feedingDurations = [intShort; intMed];  % [d]  
tEvents = sort([0;tFeedOn;tFeedOff]); 
dt = 0.5/24;              % sample time [h], converte to [d]
tGrid = (0:dt:tEnd)';     % time grid. The model outputs will be evaluated here later on
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 5*24*1e-3;  % max. feed volume flow [l/h] -> [m³/d]
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
xS = sym('x', size(xIn));   % states as col. vector
syms uS real                    % input
xiS = sym('xi', size(xIn));     % inlet concentrations (assumed known) 
thS = sym('th', size(params.th));% 6 time-variant parameters (theta)
cS = sym('c', size(params.c));  % 21 known & constant time-invariant parameters 
aS = sym('a', size(params.a));  % petersen matrix with stoichiometric constants

dynamics = BMR4_AB_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = BMR4_AB_mgl_sym(xS,cS); 

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
q = 6;   
yClean = zeros(N,q); % allocate memory
for k = 1:N
    yClean(k,:) = g(xSol(k,:)',cNum)'; % Simons Implementierung (arXiv)
end

%% add noise to measurements acc to sensor data sheets
 
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2*24*1e-3; % Trommelgaszähler Labor [L/h] --> [m³/d]
% sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCh4 = 20*0.2/100; % [Vol-%] -> [-]; erhöhe künstlich das Messrauschen; ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.12;    % NH4-N [kg/m³] = [g/L]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
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
subplot(3,2,1)
scatter(tGrid,yMeas(:,1),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,1), 'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [m³/d]')
% ylim([0,16])
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tGrid,yClean(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(3,2,3)
scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{co2} in bar')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% SIN:  
subplot(3,2,4)
scatter(tGrid,yMeas(:,4),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,4),'r',...
     'LineWidth',1.2,'DisplayName','clean')
ylabel('inorg. nitrogen in kg/m³')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% TS:  
subplot(3,2,5)
scatter(tGrid,yMeas(:,5),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tGrid,yClean(:,5),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('total solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(3,2,6)
scatter(tGrid,yMeas(:,6),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
plot(tGrid,yClean(:,6),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('volatile solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [m³/d]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

sgtitle('clean and noisy measurements from ADM1-R4')

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
fileName = 'Messung_ADM1_R4_m3.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params')
