%% Version
% (R2022b) Update 5
% Erstelldatum: 11.07.2023
% Autor: Simon Hellmann

% create synthetic measurement data for Kalman Filtering
% user can select the model 

clc
clear
close all

addpath('ODEs/')
addpath('measurementEquations/')

%% choose model type, normalization, and optionally -frac
flagModel = 4;  % 3: ADM1-R3; 4: ADM1-R4
flagNorm = 0;   % 0: absolute coordinates; 1: normalized coordinates
flagFrac = 0;   % 0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)

switch flagModel
    case 3
        model = 'ADM1_R3';
    case 4
        model = 'ADM1_R4';
end

%% define parameters like Sören: 
% Load standard model parameters
load('SoerensFiles/Model_data/ADM1_parameters.mat')
% Load experimental data
load('SoerensFiles/Model_data/ADM1_input_data.mat')
% load steady-state values from Sören's implementation: 
load(['generatedOutput/SteadyState_',num2str(model),'_Soeren.mat'])

%% inlet concentrations:
feedVolFlowSS = resultsSoeren.input(2); 

% obtain inlet concentrations and initial value for steady-state transition:
[xIn,x0SS] = getXInX0SS(flagModel,flagFrac,resultsSoeren); 
nStates = numel(xIn); 

% obtain parameters:
params = getParameters(flagModel,flagFrac,system,parameters); 
cNum = params.c; 
thNum = params.th; 
aNum = params.a; 

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
dt = 0.5/24;            % sample time [h], converte to [d]
tGrid = (0:dt:tEnd)';   % time grid. The model outputs will be evaluated here later on
N = numel(tGrid);       % # sampling points surviving at the end
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 10*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactors = [70,30]'/100; 
portions = feedFactors*feedMax; % [l/d]         	
% steady state feed volume flow [l/d] should be the average of what is fed
% during dynamic operation:
% totalFeed = sum(feedingDurations.*portions);  
feedVolFlow = zeros(nIntervals,1);  % allocate memory
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % allocate memory 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn',nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; 

% obtain system equations from symbolic formulation:
% define sizes of a, c and th: 
sza = size(aNum); 
szc = size(cNum); 
szth = size(thNum); 
[f,g] = getSystemEquations(flagModel,flagFrac,nStates,sza,szc,szth); 

%% determine steady state as initial value for simulation: 
tSpanSS = [0,tSS]; 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
% odeFunSSNum = @(t,x) ADM1_R3_ode(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
[tVecSS,xVecSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0DynR3 = xVecSS(end,:)';   % start actual simulation in steady state 
x0 = x0DynR3;         % to be overwritten

%% derive normalized system equations
q = 8;    % # of measurement signals; XY: nur fürs R3!
[fNorm,gNorm] = getSystemEquationsNorm(flagModel,flagFrac,nStates,sza,szc,szth);

% define numeric values for normalization with steady state: 
TxNum = x0; 
y0 = g(x0,cNum);    % steady state output
TyNum = y0; 
TuNum = feedVolFlowSS; % u0

% summarize in stuct to save them: 
TNum = struct; 
TNum.Tx = TxNum; 
TNum.Ty = TyNum; 
TNum.Tu = TuNum;

% normalize simulation inputs:
uNorm = feedVolFlowSS./TuNum; 
x0SSNorm = x0SS./TxNum; 
xInNorm = xIn./TxNum; 

%% simulate transition into steady state in normalized coordinates:
% XY: das ist hier nicht unbedingt nötig, denn du könntest auch direkt x0
% normieren. War eher zur Validierung gedacht!
odeFunNormSS = @(t,xNorm) fNorm(xNorm,uNorm,xInNorm,thNum,cNum,aNum,TxNum,TuNum); 
[tVecSSNorm,xSSNorm] = ode15s(odeFunNormSS,tSpanSS,x0SSNorm);
nSamplesTillSS = size(xSSNorm,1); 
x0DynNorm = xSSNorm(end,:)'; % extract only last value as actual steady state
x0Norm = x0DynNorm;
 
% compute state trajectories during dynamic feeding scenario:
if flagNorm == 1
    % case a: in normalized coordinates:
    xSolNorm = computeXTrajectoriesNorm(x0Norm,fNorm,inputMat,tEvents,tOverall,tGrid,tEnd,thNum,cNum,aNum,nStates,TxNum,TuNum); 
    xSol = repmat(TxNum',N,1).*xSolNorm;    % de-normalize states
else
    % case b: in absolute coordinates:
    xSol = computeXTrajectories(x0,f,inputMat,tEvents,tOverall,tGrid,tEnd,thNum,cNum,aNum,nStates); 
end

%% compute output variables
if flagNorm == 1
    yCleanNorm = zeros(N,q); % allocate memory
    for k = 1:N
        yCleanNorm(k,:) = gNorm(xSolNorm(k,:)', cNum, TxNum, TyNum);
    end
    
    % de-normalize outputs:
    yClean = repmat(TyNum',N,1).*yCleanNorm;
else 
    yClean = zeros(N,q); % allocate memory
    for k = 1:N
        yClean(k,:) = g(xSol(k,:)',cNum)'; % Simons Implementierung (arXiv)
    end
end

% add measurement noise to clean outputs (adapt sigmas therein):
[yMeas,noiseCovMat] = addNoiseToMeasurements(flagModel,yClean); 

% plot results:
plotOutputs(flagModel,flagFrac,yClean,yMeas,feedVolFlow,tGrid,tEvents)

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0DynR3;  % steady state as initial value for dynamic simulation
% MESS.xSim = [tSim,xSimNorm]; 
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = yClean; 
MESS.yMeas = yMeas; 
MESS.R = noiseCovMat; % accurate values from sensor data sheets

if flagNorm == 1
    MESS.x0Norm = x0DynNorm;   % normalized initial state
    MESS.xNorm = xSolNorm; 
    MESS.yCleanNorm = yCleanNorm; % normalized outputs
end

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
mkdir(pathToResults);   % create subfolder (gives warning if exists yet)
fileName = 'Messung_ADM1_R3_norm.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params', 'TNum')
