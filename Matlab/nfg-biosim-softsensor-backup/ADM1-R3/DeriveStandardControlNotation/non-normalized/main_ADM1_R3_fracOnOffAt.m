%% Version
% (R2022b) Update 5
% Erstelldatum: 16.6.2023
% Autor: Simon Hellmann

% run the ADM1-R3-frac in a realistic scenario with online/offline/atline
% measurements

clc
clear
close all

%% define parameters like Sören: 
% Load standard model parameters
load('Model_data\ADM1_parameters.mat')
% Load experimental data
load('Model_data\ADM1_input_data.mat')
% load steady-state values from Sören's implementation: 
load('SteadyState_ADM1-R3_Soeren.mat')

% renaming for easier understanding: 
s = system.Variables;   % s = systemParameters
systemInput = input.ADM1_R3.Variables;
parameters = parameters.ADM1_R3.Variables; % modelParameters

% extract system parameter values: 
Vl = s(1);  % gas volume
Vg = s(2);  % liquid volume
% p0 = s(3);  % atmospheric pressure [bar] % XY: Achtung: Typo im GitHub
p0 = 1.0133; 

% extract model parameter values:  
Kch4 =  parameters(1);  % Henry parameter ch4 [mol/l/bar]
Kco2 =  parameters(2);  % Henry parameter co2 [mol/l/bar]
KSIN =  parameters(3);  % half-saturation constant nitrogen limitation [g/l]
KInh3 =  parameters(4); % ammonia inbibition constant [g/l]
% KaIN =  parameters(5);  % dissociation constant ammonia [mol/l]
% Kaac =  parameters(6);  % dissociation constant acetic acid [mol/l]
% Kaco2 =  parameters(7); % dissociation constant carbonate [mol/l]
% corrected values acc. to Simons computation (Sören made slight mistake):
KaIN = 1.3809E-9; 
Kaac = 1.7378E-5; 
Kaco2 = 5.0981E-7;
KSac =  parameters(8);  % half-saturation constant acetoclastic methanogenesis [g/l]
Kw =  parameters(9);    % ion product of water [mol/l]
R =  parameters(10);    % universal gas constant [bar l/mol/K]
T =  parameters(11);    % operating temperature [K]
kABIN =  parameters(12);    % kin. dissociation constant ammonia [l/mol/d] 
kABac =  parameters(13);    % kin. dissociation constant acetic acid [l/mol/d]
kABco2 =  parameters(14);   % kin. dissociation constant carbonate [l/mol/d]
kLa =  parameters(15);  % mass transfer coefficient [1/d]
kch =  parameters(16);  % hydrolysis constant carbohydrates [1/d]
kdec =  parameters(17); % hydrolysis constant biomass decay [1/d]
kli =  parameters(18);  % hydrolysis constant lipids [1/d]
muM =  parameters(19);  % max. growth rate [1/d]
kp =  parameters(20);   % friction parameter [l/bar/d]
kpr =  parameters(21);  % hydrolysis constant proteins [1/d]
pHLL =  parameters(22); % lower pH boundary  
pHUL =  parameters(23); % upper pH boundary
ph2o =  parameters(24); % partial pressure of water in gas phase (saturated) [bar]

rho = 1000;     % mass density of digestate [g/l]
Mch4 = 16;      % molar mass CH4 [kg/kmol]
Mco2 = 44;      % molar mass CO2 [kg/kmol]
nac = 3/(pHUL - pHLL); 

fracChFast = 1; % fraction of fast inlet cabohydrates

%% order model parameters in the rights structures (prepare simulation)
c1 = 1/Vl; 
c2 = nac; 
c3 = 10^(-(3/2)*(pHUL + pHLL)/(pHUL - pHLL)); 
c4 = 4*Kw; 
c5 = kLa; 
c6 = kLa*Kch4*R*T; 
c7 = kLa*Kco2*R*T; 
c8 = KSIN; 
c9 = kABac;
c10 = kABco2; 
c11 = kABIN; 
c12 = kLa*Vl/Vg; 
c13 = kp/p0*(R*T/Mch4)^2;
c14 = 2*kp/p0*(R*T)^2/Mch4/Mco2;
c15 = kp/p0*(R*T/Mco2)^2;
c16 = kp/p0*R*T/Mch4*(2*ph2o - p0); 
c17 = kp/p0*R*T/Mco2*(2*ph2o - p0); 
c18 = kp/p0*(ph2o - p0)*ph2o; 
c19 = R*T/Mch4;
c20 = R*T/Mco2;
c21 = rho; 
c22 = -kp/Vg/p0*(R*T/Mch4)^2;
c23 = -2*kp/Vg/p0*(R*T)^2/Mch4/Mco2;
c24 = -kp/Vg/p0*(R*T/Mco2)^2;
c25 = -kp/Vg/p0*(R*T/Mch4)*(2*ph2o - p0);
c26 = -kp/Vg/p0*(R*T/Mco2)*(2*ph2o - p0);
c27 = -kLa*Vl/Vg*Kch4*R*T - kp/Vg/p0*(ph2o - p0)*ph2o;
c28 = -kLa*Vl/Vg*Kco2*R*T - kp/Vg/p0*(ph2o - p0)*ph2o;
c29 = kABac*Kaac; 
c30 = kABco2*Kaco2;
c31 = kABIN*KaIN; 
c32 = Vl/Vg;
c33 = -kp/Vg/p0*(ph2o - p0)*ph2o;
% combine all in column vector: 
cNum = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,...
        c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33]';

% petersen matrix acc. to Weinrich 2021
% Note:all aij with |aij ~= 1|, only take the absolute value (see arXiv paper). 
% Index "Num" for Numeric
aNum = [0.6555, 0.081837, 0.2245,  0.016932, 0.057375, -1,      0,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
        0.6555, 0.081837, 0.2245,  0.016932, 0.057375,  0,     -1,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
        0.9947, 0.069636, 0.10291, 0.17456,  0.47666,   0,      0,     -1,      0,      0.13486,0, 0, 0, 0, 0, 0, 0,    0;
        1.7651, 0.19133,  0.64716, 0.024406, 0.44695,   0,      0,      0,     -1,      0.1621, 0, 0, 0, 0, 0, 0, 0,    0;
        26.5447,6.7367,  18.4808,  0.15056,  0.4778,    0,      0,      0,      0,      0,      1, 0, 0, 0, 0, 0, 0,    0; 
        0,      0,        0,       0,        0,         0.18,   0,      0.77,   0.05,  -1,      0, 0, 0, 0, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0.18,   0,      0.77,   0.05,   0,     -1, 0, 0, 0, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0,-1, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0,-1, 0, 0,    0; 
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0,-1, 0,    0;
        0,     -1,        0,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0, 0, c32,  0; 
        0,      0,       -1,       0,        0,         0,      0,      0,      0,      0,      0, 0, 0, 0, 0, 0, 0,    c32;]';

%% inlet concentrations:
feedVolFlowSS = resultsSoeren.input(2); 
xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
% beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
% fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
% der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
% ist es folgendermaßen:
ScatIn = xInPre(11); 
SanIn = xInPre(12);
SionIN = ScatIn - SanIn; 
% adapt inlet concentrations for slightly different state indexing in std. 
% control notation of ADM1-R3-frac (X_chS/F, X_ash, S_ion = S_cat - S_an): 
nStates = length(xInPre) + 1;    % add second CH fraction
XInCh = xInPre(6);      % all carbohydrates in
xIn = zeros(nStates,1); % allocate memory
xIn(1:5) = xInPre(1:5)';    % all states up to h2o
xIn(8:end) = xInPre(7:end)';% all other states except ch
% assign all carbohydrates_in to the fast fraction, distribute them later 
% in the ode model file:
xIn(6) = XInCh; 
xIn(7) = 0; 
% overwrite values in positions for S_an and S_cat with those for S_ion and X_ash:
xAshIn = 14; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)
xIn(12) = xAshIn;  
xIn(13) = SionIN; 

%% initial condition: 
x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
Scat0 = x0Pre(11); 
San0 = x0Pre(12); 
Sion0 = Scat0 - San0; 
% adapt initial condition for slightly different state indexing in
% Simon's model (X_ash, S_ion = S_cat - S_an): 
x0SS = zeros(nStates,1); % x0Pre'
Xch0 = x0Pre(6);   
x0SS(1:5) = x0Pre(1:5); 
x0SS(6) = fracChFast*Xch0; 
x0SS(7) = (1-fracChFast)*Xch0;
x0SS(8:end) = x0Pre(7:end); 
% overwrite values in positions for S_an and S_cat with those for S_ion and X_ash:
xAsh0 = 1; % assume 1 g/l initial ash concentration
x0SS(12) = xAsh0;
x0SS(13) = Sion0;   

%% miscellaneous parameters for simulation
% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
kchS = 1e-1*kch; 
kchF = kch; 
thNum = [kchF, kchS, kpr, kli, kdec, muM, KSac, KInh3, fracChFast]';
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
dt = 1/60/24;             % sample time [min], converted to [d]
tOnline = (0:dt:tEnd)';   % time grid for online measurements. highest frequency in this script
tOverall = unique([tOnline; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 10*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactors = [70,30]'/100; 
portions = feedFactors*feedMax; % [l/d]         	
% steady state feed volume flow [l/d] should be the average of what is fed
% during dynamic operation:
totalFeed = sum(feedingDurations.*portions);  
feedVolFlow = zeros(nIntervals,1);  % allocate memory
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % allocate memory 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn',nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; 

%% derive odeFun from symbolic model definition and determine steady state
% define symbolic ("S") variables (all vector are defined as column vectors)
xS = sym('x', [18 1]);   % states as col. vector
syms uS real             % input
xiS = sym('xi', [18,1]); % inlet concentrations (assumed known) 
thS = sym('th', [9,1]);  % time-variant parameters (theta)
cS = sym('c', [33,1]);   % known & constant time-invariant parameters 
aS = sym('a', [18,12]);% petersen matrix with stoichiometric constants

dynamics = ADM1_R3_frac_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = ADM1_R3_frac_mgl_sym(xS,cS); 

% transform into numeric function handle. Note that the independet
% variables are explicitely defined. Their order must be followed when 
% using the function handle!
f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
g = matlabFunction(outputs, 'Vars', {xS, cS}); 

%% determine steady state as initial value for simulation:  
tSpanSS = [0,tSS]; 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
% odeFunSSNum = @(t,x) ADM1_R3_ode(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
[tVecSS,xVecSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0DynR3frac = xVecSS(end,:)';   % start actual simulation in steady state
x0 = x0DynR3frac;         % to be overwritten

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
    t_ode = tOverall(idxTimeInterval); 
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

%% compute output variables
% case 1: only online measurements
% Evaluate xSol only in tOnline, discard the rest
idxGridOn = ismember(tOverall, tOnline); 
xSolOn = xSim(idxGridOn,:);
NOn = size(xSolOn,1);       % number of survived sampling points
q = 8;   
yClean = zeros(NOn,q);    % allocate memory
for k = 1:NOn
    yClean(k,:) = g(xSolOn(k,:)',cNum)'; % Simons Implementierung (arXiv)
end

% case 2: online, and delayed offline and atline measurements
dtOff = 1;      % sample time for offline measurements [d]
dtAt = 1/24;    % sample time for atline measurements [d]
% times when samples were taken:
tOfflineSample = (0.45:dtOff:tEnd)';  % offset offline from online measurements 
tAtlineSample = (0.25:dtAt:tEnd)';   % offset atline from online measurements 

% interpolate xSim at online/offline/atline sample times: 
xSolOff = interp1(tOverall,xSim,tOfflineSample);
xSolAt = interp1(tOverall,xSim,tAtlineSample);

NOff = numel(tOfflineSample); 
NAt = numel(tAtlineSample); 

% number of measurement signals: 
qOn = 4; 
qOff = 3; 
qAt = 1; 
% allocate memory:
yCleanOn = zeros(NOn,qOn);
yCleanOff = zeros(NOff,qOff);
yCleanAt = zeros(NAt,qAt);
% evaluate model outputs at different sampling times:
for k = 1:NOn
    fullOutput = g(xSolOn(k,:)',cNum)';
    yCleanOn(k,:) = fullOutput(1:qOn); 
end
for kk = 1:NOff
    fullOutput = g(xSolOff(kk,:)',cNum)';
    yCleanOff(kk,:) = fullOutput(qOn+1:qOn+qOff); 
end
for kkk = 1:NAt
    fullOutput = g(xSolAt(kkk,:)',cNum)';
    yCleanAt(kkk,:) = fullOutput(qOn+qOff+1:end);
end

%% add time delays for return times of sample times
% construct times when offline and atline measurements return from lab: 
delayMin = 0.8;     % [d]
delaySpread = 1/7;  % [d] (manually chosen) 
rng('default');     % fix seed for random number generation (for replicable results)
delayOff = delayMin + rand(NOff,1)*delaySpread;
delayAt = delayMin + rand(NAt,1)*delaySpread;
tOfflineRet = tOfflineSample + delayOff;  % Ret = return
tAtlineRet = tAtlineSample + delayAt;     % Ret = return

% optional clipping of return times by tEnd: 
% idxCliptOff = tOfflineRet <= tEnd; 
% idxCliptAt = tAtlineRet <= tEnd; 
% tOfflineRetClip = tOfflineRet(idxCliptOff); 
% tAtlineRetClip = tAtlineRet(idxCliptAt); 
% % correct remaining number of elements
% NOffClip = numel(tOfflineRetClip); 
% NAtClip = numel(tAtlineRetClip); 

% incorporate knowledge of sampling sequence and incoming sequence: 
[tOfflineRetSorted,indexOffline] = sort(tOfflineRet); 
[tAtlineRetSorted,indexAtline] = sort(tAtlineRet); 
trueIndexOffline = 1:NOff; 
trueIndexAtline = 1:NAt;
indexDevOffline = indexOffline' - trueIndexOffline; 
indexDevAtline = indexAtline' - trueIndexAtline; 

%% add noise to measurements acc to sensor data sheets
% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2*24;    % Trommelgaszähler Labor [L/h] --> [L/d]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaPh = 0.1;      % pH [-]
sigmaSIN = 0.12;    % NH4-N [g/L]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]
sigmaAc = 0.04;     % FOS [g/L]. Aber Achtung: brauchbarere Messgröße für 
% S_ac ist eher die Einzelsäure aus GC, diese hat aber sigma=0,01 g/L 

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh, sigmaSIN, sigmaTS, sigmaVS, sigmaAc]; 
sigmasOn = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh]; 
sigmasOff = [sigmaSIN, sigmaTS, sigmaVS];
sigmasAt = [sigmaAc]; 

sigmaMat = repmat(sigmas,NOn,1);
sigmaMatOn = repmat(sigmasOn,NOn,1);
sigmaMatOff = repmat(sigmasOff,NOff,1);
sigmaMatAt = repmat(sigmasAt,NAt,1);

% create normally distributed measurement noise matrices:
% zero mean for all measurements 
yMean = zeros(NOn,q); 
yMeanOn = zeros(NOn,qOn); 
yMeanOff = zeros(NOff,qOff); 
yMeanAt = zeros(NAt,qAt); 
normalMeasNoise = normrnd(yMean,sigmaMat);
normalMeasNoiseOn = normrnd(yMeanOn,sigmaMatOn);
normalMeasNoiseOff = normrnd(yMeanOff,sigmaMatOff);
normalMeasNoiseAt = normrnd(yMeanAt,sigmaMatAt);

% add noise to clean model outputs:
yMeas = yClean + normalMeasNoise; 
yMeasOn = yCleanOn + normalMeasNoiseOn; 
yMeasOff = yCleanOff + normalMeasNoiseOff; 
yMeasAt = yCleanAt + normalMeasNoiseAt; 

% construct measurement noise covariance matrices:
noiseCovMat = diag(sigmas.^2);  
noiseCovMatOn = diag(sigmasOn.^2);  
noiseCovMatOff = diag(sigmasOff.^2); 
noiseCovMatAt = diag(sigmasAt.^2); 

%% Plot results (separated into Online/Offline/Atline Measurements)
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figure()

% online measurements:
% gas volume flow: 
subplot(4,2,1)
scatter(tOnline,yMeasOn(:,1)/24,'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yCleanOn(:,1)/24,'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
ylabel('gas vol flow [l/h]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
title('Online')

% pch4: 
subplot(4,2,2)
scatter(tOnline,yMeasOn(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(tOnline,yCleanOn(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(4,2,3)
scatter(tOnline,yMeasOn(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yCleanOn(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('p_{co2} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pH:  
subplot(4,2,4)
scatter(tOnline,yMeasOn(:,4),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
plot(tOnline,yCleanOn(:,4),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('pH value [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Offline measurements:
% SIN:  
subplot(4,2,5)
scatter(tOfflineSample,yMeasOff(:,1),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOfflineSample,yCleanOff(:,1),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('SIN [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')
title('Offline')

% TS:  
subplot(4,2,6)
scatter(tOfflineSample,yMeasOff(:,2),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
hold on; 
plot(tOfflineSample,yCleanOff(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('total solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(4,2,7)
scatter(tOfflineSample,yMeasOff(:,3),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
plot(tOfflineSample,yCleanOff(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('volatile solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Atline measurements:
% Sac:  
subplot(4,2,8)
scatter(tAtlineSample,yMeasAt(:,1),'DisplayName','noisy',...
        'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on; 
plot(tAtlineSample,yCleanAt(:,1),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylabel('acetic acid [g/l]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')
title('Atline')

sgtitle('Clean and noisy simulation outputs (4 Online/3 Offline/1 Atline)')

%% save results in struct MESS: 
MESS.tOnline = tOnline;
MESS.tOfflineSend = tOfflineSample; % sampling times 
MESS.tAtlineSend = tAtlineSample;     
MESS.tOfflineRet = tOfflineRet;     % return times  
MESS.tAtlineRet = tAtlineRet; 
MESS.x0 = x0DynR3frac; 
MESS.xSolOn = xSolOn;   % state trajectories evaluated at diff. sampling times
MESS.xSolOff = xSolOff; 
MESS.xSolAt = xSolAt; 
% MESS.xSim = [tSim,xSim]; 
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yCleanOn = yCleanOn;  
MESS.yCleanOff = yCleanOff;  
MESS.yCleanAt = yCleanAt;  
MESS.yMeasOn = yMeasOn; 
MESS.yMeasOff = yMeasOff; 
MESS.yMeasAt = yMeasAt; 
MESS.C = noiseCovMat; % accurate values from sensor data sheets
MESS.COn = noiseCovMatOn;
MESS.COff = noiseCovMatOff;
MESS.CAt = noiseCovMatAt;

save('Messung_ADM1_R3_frac_OnOffAt.mat', 'MESS', 'params')