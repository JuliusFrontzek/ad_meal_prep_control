%% Version
% (R2022b) Update 5
% Erstelldatum: 25.09.2023
% Autor: Simon Hellmann

% create synthetic measurement data for Kalman Filtering
% scenario: online/offline rates, but offline only with single augmentation
% model: ADM1-R3-frac-norm (including ash and 2 carbohydrate fractions) 

clc
clear
close all

%% define parameters like Sören: 
% Load standard model parameters
load('SoerensFiles\Model_data\ADM1_parameters.mat')
% Load experimental data
load('SoerensFiles\Model_data\ADM1_input_data.mat')
% load steady-state values from Sören's implementation: 
load('generatedOutput\SteadyState_ADM1-R3_Soeren.mat')

% renaming for easier understanding: 
s = system.Variables;   % s = systemParameters
systemInput = input.ADM1_R3.Variables;
parameters = parameters.ADM1_R3.Variables; % modelParameters

% extract system parameter values: 
V_liq = s(1);  % gas volume
V_gas = s(2);  % liquid volume
% p0 = s(3);  % atmospheric pressure [bar] XY: Achtung: Typo im GitHub
p_atm = 1.0133; 

% extract model parameter values:  
K_H_ch4 =  parameters(1);  % Henry parameter ch4 [mol/l/bar]
K_H_co2 =  parameters(2);  % Henry parameter co2 [mol/l/bar]
K_S_IN =  parameters(3);  % half-saturation constant nitrogen limitation [g/l]
K_I_nh3 =  parameters(4); % ammonia inbibition constant [g/l]
% KaIN =  parameters(5);  % dissociation constant ammonia [mol/l]
% Kaac =  parameters(6);  % dissociation constant acetic acid [mol/l]
% Kaco2 =  parameters(7); % dissociation constant carbonate [mol/l]
% corrected values acc. to Simons computation (Sören made slight mistake):
K_a_IN = 1.3809E-9; 
K_a_ac = 1.7378E-5; 
K_a_co2 = 5.0981E-7;
K_S_ac =  parameters(8);  % half-saturation constant acetoclastic methanogenesis [g/l]
K_w =  parameters(9);    % ion product of water [mol/l]
R =  parameters(10);    % universal gas constant [bar l/mol/K]
T =  parameters(11);    % operating temperature [K]
k_AB_IN =  parameters(12);    % kin. dissociation constant ammonia [l/mol/d] 
k_AB_ac =  parameters(13);    % kin. dissociation constant acetic acid [l/mol/d]
k_AB_co2 =  parameters(14);   % kin. dissociation constant carbonate [l/mol/d]
k_La =  parameters(15);  % mass transfer coefficient [1/d]
k_ch =  parameters(16);  % hydrolysis constant carbohydrates [1/d]
k_dec =  parameters(17); % hydrolysis constant biomass decay [1/d]
k_li =  parameters(18);  % hydrolysis constant lipids [1/d]
k_m_ac =  parameters(19);  % max. growth rate [1/d]
k_p =  parameters(20);   % friction parameter [l/bar/d]
k_pr =  parameters(21);  % hydrolysis constant proteins [1/d]
pK_l_ac =  parameters(22); % lower pH boundary  
pK_u_ac =  parameters(23); % upper pH boundary
p_h2o =  parameters(24); % partial pressure of water in gas phase (saturated) [bar]

rho = 1000;        % mass density of digestate [g/l]
Mch4 = 16;      % molar mass CH4 [kg/kmol]
Mco2 = 44;      % molar mass CO2 [kg/kmol]
n_ac = 3/(pK_u_ac - pK_l_ac); 

fracChFast = 0.5; % fraction of fast cabohydrates Rindergülle (rel. hoher Faseranteil am Eingang)

%% order model parameters in the rights structures (prepare simulation)
c1 = 1/V_liq; 
c2 = n_ac; 
c3 = 10^(-(3/2)*(pK_u_ac + pK_l_ac)/(pK_u_ac - pK_l_ac)); 
c4 = 4*K_w; 
c5 = k_La; 
c6 = k_La*K_H_ch4*R*T; 
c7 = k_La*K_H_co2*R*T; 
c8 = K_S_IN; 
c9 = k_AB_ac;
c10 = k_AB_co2; 
c11 = k_AB_IN; 
c12 = k_La*V_liq/V_gas; 
c13 = k_p/p_atm*(R*T/Mch4)^2;
c14 = 2*k_p/p_atm*(R*T)^2/Mch4/Mco2;
c15 = k_p/p_atm*(R*T/Mco2)^2;
c16 = k_p/p_atm*R*T/Mch4*(2*p_h2o - p_atm); 
c17 = k_p/p_atm*R*T/Mco2*(2*p_h2o - p_atm); 
c18 = k_p/p_atm*(p_h2o - p_atm)*p_h2o; 
c19 = R*T/Mch4;
c20 = R*T/Mco2;
c21 = rho; 
c22 = -k_p/V_gas/p_atm*(R*T/Mch4)^2;
c23 = -2*k_p/V_gas/p_atm*(R*T)^2/Mch4/Mco2;
c24 = -k_p/V_gas/p_atm*(R*T/Mco2)^2;
c25 = -k_p/V_gas/p_atm*(R*T/Mch4)*(2*p_h2o - p_atm);
c26 = -k_p/V_gas/p_atm*(R*T/Mco2)*(2*p_h2o - p_atm);
c27 = -k_La*V_liq/V_gas*K_H_ch4*R*T - k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;
c28 = -k_La*V_liq/V_gas*K_H_co2*R*T - k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;
c29 = k_AB_ac*K_a_ac; 
c30 = k_AB_co2*K_a_co2;
c31 = k_AB_IN*K_a_IN; 
c32 = V_liq/V_gas;
c33 = -k_p/V_gas/p_atm*(p_h2o - p_atm)*p_h2o;
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
xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
% beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
% fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
% der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
% ist es folgendermaßen:
ScatIn = xInPre(11);
SanIn = xInPre(12); 
SionIN = ScatIn - SanIn; 
% adapt inlet concentrations for slightly different state indexing in
% Simon's ADM1-R3-frac model (X_chS/F, X_ash, S_ion = S_cat - S_an): 
nStates = length(xInPre) + 1;               % add second CH fraction
XInCh = xInPre(6);                          % all carbohydrates in
xIn = nan(nStates,1);                       % allocate memory
xIn([1:5,8:end]) = xInPre([1:5,7:end])';    % all states except ch
% assign all carbohydrates_in to the fast fraction, distribute them later 
% in the ode model file:
xIn(6) = XInCh; 
xIn(7) = 0; 
% overwrite values in positions for S_an and S_cat with those for S_ion and X_ash:
xAshIn = 17; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)
xIn(12) = xAshIn;  
xIn(13) = SionIN; 

%% initial condition: 
x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
Scat0 = x0Pre(11); 
San0 = x0Pre(12); 
Sion0 = Scat0 - San0; 
% adapt initial condition for slightly different state indexing in standard 
% control notation (X_ash, S_ion = S_cat - S_an): 
x0SS = nan(nStates,1); % x0Pre'
Xch0 = x0Pre(6);   
x0SS([1:5,8:end]) = x0Pre([1:5,7:end]);     % all states except ch
x0SS(6) = fracChFast*Xch0; 
x0SS(7) = (1-fracChFast)*Xch0;
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
k_ch_S = 1e-1*k_ch; 
k_ch_F = k_ch; 
thNum = [k_ch_F, k_ch_S, k_pr, k_li, k_dec, k_m_ac, K_S_ac, K_I_nh3, fracChFast]';
params.th = thNum; 

% times: 
tEnd = 7;   % [d] End of Simulation
tSS = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.25;    % [d]
intLong = 1.5;      % [d]
intMed = 0.5;       % [d]
% start with feeding right after steady state and some pause: 
% ints = [intLong,intShort,intLong,intMed]';
% ... and transform intervals to absolute times:
% cumInts = cumsum(ints);     % cumulated intervals
tFeedOnWeek = [1;2.5;4.5;6];  % beginning times of feedings (1 single week)
tFeedOn = [tFeedOnWeek];% [tFeedOnWeek;tFeedOnWeek + 7]; 
feedingDurationsWeek = [intMed; intShort; intMed; intShort];% 1 single week
feedingDurations = repmat(feedingDurationsWeek,1,1);        % keep 1 week
tFeedOff = tFeedOn + feedingDurations; % end times of feedings
tEvents = sort([0;tFeedOn;tFeedOff]); 
dt = 20/60/24;             % sample time [min], converted to [d]. Vorher: 5
tOnline = (0:dt:tEnd)';   % time grid for online measurements. highest frequency in this script
tMinor = unique([tOnline; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
nIntervals = length(tEvents); 
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 10*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactorsWeek = [70,30,50,40]'/100;       % 1 single week
feedFactors = repmat(feedFactorsWeek,1,1);  % mulitiples of 1 week
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
feedVolFlowSS = resultsSoeren.input(2); 
odeFunSS = @(t,x) f(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
% odeFunSSNum = @(t,x) ADM1_R3_ode(x,feedVolFlowSS,xIn,thNum,cNum,aNum); 
[tVecSS,xVecSS] = ode15s(odeFunSS,tSpanSS,x0SS); 
x0Init = xVecSS(end,:)';   % start actual simulation in steady state
x0 = x0Init;         % to be overwritten

%% derive normalized system equations
xNormS = sym('xNorm', [18 1]);  % normalized states as col. vector
syms uNorm real;                % normalized input
xiNormS = sym('xi', [18,1]);    % normalized inlet concentrations 
TxS = sym('Tx', [18,1]);        % normalization matrix for states
TyS = sym('Ty', [8,1]);         % normalization matrix for outputs
syms Tu real                    % normalization variable for input

dynamicsNorm = ADM1_R3_frac_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = ADM1_R3_frac_norm_mgl_sym(xNormS, cS, TxS, TyS); 

% turn into numeric function handles: 
fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 

% define numeric values for normalization with steady state: 
TxNum = x0; 
y0 = g(x0,cNum);    % steady state output
TyNum = y0; 
TuNum = feedVolFlowSS; 
% summarize in stuct to save them: 
TNum        = struct; 
TNum.Tx     = TxNum; 
TNum.Ty     = TyNum; 
TNum.Tu     = TuNum;    

% normalize simulation inputs:
uNorm = feedVolFlowSS./TuNum; 
x0SSNorm = x0SS./TxNum; 
xInNorm = xIn./TxNum; 

% simulate transition into steady state in normalized coordinates: 
odeFunNormSS = @(t,xNorm) fNorm(xNorm,uNorm,xInNorm,thNum,cNum,aNum,TxNum,TuNum); 
[tVecSSNorm,xSSNorm] = ode15s(odeFunNormSS,tSpanSS,x0SSNorm);
nSamplesTillSS = size(xSSNorm,1); 
x0InitNorm = xSSNorm(end,:)'; % extract only last value as actual steady state
x0Norm = x0InitNorm;

% compute normalized system output at steady state: 
ySSNorm = gNorm(x0InitNorm, cNum, TxNum, TyNum); 

% perform de-normalization to check if normalization works properly: 
xSSDeNorm = repmat(TxNum',nSamplesTillSS,1).*xSSNorm;
x0DynDeNorm = xSSDeNorm(end,:)';
ySSDeNorm = TyNum.*ySSNorm;     % compare relevant entries with Sören's SS:
ySSSoeren = resultsSoeren.ySS';

%% Solve ODE via iterative solution of constant feeding regimes (on or off)
xSimNorm = nan(length(tMinor), nStates);% allocate memory
tSim = nan(length(tMinor),1);       % allocate memory

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
    uCurrNorm = feedVolFlowCurr./TuNum; 
    xInCurrNorm = xInCurr./TxNum;

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

%% case 1: compute all-online output variables
% Note: all-online measurements (unrealistic)
% Evaluate xSol only in tOnline, discard the rest
q = 8;   
idxGridOn = ismember(tMinor, tOnline); 
xSolOnNorm = xSimNorm(idxGridOn,:);
NOn = size(xSolOnNorm,1);   % number of online sampling points
yCleanNorm = nan(NOn,q);    % allocate memory
for k = 1:NOn
    yCleanNorm(k,:) = gNorm(xSolOnNorm(k,:)', cNum, TxNum, TyNum)';
end

% de-normalize states and outputs:
xSolOn = repmat(TxNum',NOn,1).*xSolOnNorm;  % states
yClean = repmat(TyNum',NOn,1).*yCleanNorm;  % outputs

%% case 2: online, and delayed offline and atline measurements
% standard measurements: 
dtOff = 1;      % sample time for offline measurements [d]
% times when samples were taken:
tOfflineSample = (0.45:dtOff:tEnd)';  % offset offline from online measurements 
NOff = numel(tOfflineSample); % # offline sample points

% interpolate xSimNorm at offline sample times: 
xSolOffNorm = interp1(tMinor,xSimNorm,tOfflineSample);
% de-normalize offline evaluation of states: 
xSolOff = repmat(TxNum',NOff,1).*xSolOffNorm;

qOn = 4;  % # online measurement signals
qOff = 4; % # offline measurement signals
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
    yCleanOffNorm(kk,:) = fullOutputNorm(qOn+1:end); 
end
% de-normalize:
yCleanOff = repmat(TyNum(qOn+1:end)',NOff,1).*yCleanOffNorm;

% add time delays for arrival times of samplings
% construct times when offline measurements return from lab: 
delayMin = 0.5;     % [d] (manually chosen)
delaySpread = 2/24; % [min], converted to [d] (manually chosen) 
rng('default');     % fix seed for random number generation (for replicable results)
delayOff = delayMin + rand(NOff,1)*delaySpread;
tOfflineArrival = tOfflineSample + delayOff;

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
sigmaAc = 5*0.04;   % FOS [g/L], manuell verstärkt. Achtung: brauchbarere Messgröße für 
% S_ac ist eher die Einzelsäure aus GC, diese hat aber sigma=0,01 g/L 

% combine all in sigma matrix and covariance matrix:
sigmas = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh, sigmaSIN, sigmaTS, sigmaVS, sigmaAc]; 
sigmasOn = [sigmaV, sigmaCh4, sigmaCo2, sigmaPh]; 
sigmasOff = [sigmaSIN, sigmaTS, sigmaVS, sigmaAc];

sigmaMat = repmat(sigmas,NOn,1);
sigmaMatOn = repmat(sigmasOn,NOn,1);
sigmaMatOff = repmat(sigmasOff,NOff,1);

% create normally distributed measurement noise matrices:
% zero mean for all measurements 
yMean = zeros(NOn,q); 
yMeanOn = zeros(NOn,qOn); 
yMeanOff = zeros(NOff,qOff); 
normalMeasNoise = normrnd(yMean,sigmaMat);
normalMeasNoiseOn = normrnd(yMeanOn,sigmaMatOn);
normalMeasNoiseOff = normrnd(yMeanOff,sigmaMatOff);

% add noise to clean model outputs:
yMeas = yClean + normalMeasNoise; 
yMeasOn = yCleanOn + normalMeasNoiseOn; 
yMeasOff = yCleanOff + normalMeasNoiseOff; 

% construct measurement noise covariance matrices:
noiseCovMat = diag(sigmas.^2);  
noiseCovMatOn = diag(sigmasOn.^2);  
noiseCovMatOff = diag(sigmasOff.^2); 

%% construct united measurement array for Kalman filtering with NaNs where 
% no (offline) measurement was returned
MEASUnite = nan(NOn,qOn+qOff);     % allocate memory

% shift/pad sampling and arrival times of offline measurements into fine 
% time grid of tMinor: 
tOfflineSampleShift = interp1(tMinor,tMinor,tOfflineSample,'next'); 
tOfflineArrivalShiftPre = interp1(tMinor,tMinor,tOfflineArrival,'next'); % includes NaN where extrapolation would be necessary
idxKeepOfflineMeasurements = ~isnan(tOfflineArrivalShiftPre); % remove NaNs
tOfflineArrivalShift = tOfflineArrivalShiftPre(idxKeepOfflineMeasurements); % remove NaNs

% drop those offline measurements that return after end of simulation
yMeasOffEff = yMeasOff(idxKeepOfflineMeasurements,:);   % effective offline measurements
% insert minor instances in first qOn columns: 
MEASUnite(:,1:qOn) = yMeasOn; 
% insert minor instances (at right timing) in last qOff cols of MEASUnite:
[~,idxTOfflineArr] = ismember(tOfflineArrivalShift,tMinor); % determine correct row positioning of major instances in MEASUnite (right timing)
MEASUnite(idxTOfflineArr,qOn+1:end) = yMeasOffEff; 

%% Plot results (separated into Online/Offline/Atline Measurements)
colorPaletteHexMagma = ["#fcfdbf","#fc8961","#b73779","#51127c","#000004"];
% plot and compare the clean results with noisy measurements: 
figOutputs = figure; 

% online measurements:
% gas volume flow: 
subplot(4,2,1)
scatter(tOnline,yMeasOn(:,1)/24,'DisplayName','noisy measurements',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yCleanOn(:,1)/24,'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5); 
ylabel('gas vol flow [l/h]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
title('Online')

% pch4: 
subplot(4,2,2)
scatter(tOnline,yMeasOn(:,2),'DisplayName','noisy',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5);
hold on
plot(tOnline,yCleanOn(:,2),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
ylabel('p_{ch4} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pco2:
subplot(4,2,3)
scatter(tOnline,yMeasOn(:,3),'DisplayName','noisy',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5); 
hold on; 
plot(tOnline,yCleanOn(:,3),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
ylabel('p_{co2} [bar]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% pH:  
subplot(4,2,4)
scatter(tOnline,yMeasOn(:,4),'DisplayName','noisy',...
        'Marker','.', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1.5)
hold on; 
plot(tOnline,yCleanOn(:,4),'DisplayName','clean',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
ylabel('pH value [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Standard measurements (offline) (show sample and arriva times!)
% SIN:  
subplot(4,2,5)
plot(tMinor,yClean(:,5),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
hold on; 
% measurements at their sampling times:
scatter(tOfflineSampleShift,yMeasOff(:,1),'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
% measurements at their arrival times: 
scatter(tOfflineArrivalShift,yMeasOffEff(:,1),'DisplayName','noisy arrivals',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
ylabel('SIN [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')
title('Offline')

% TS:  
subplot(4,2,6)
plot(tMinor,yClean(:,6),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
hold on; 
% measurements at their sampling times:
scatter(tOfflineSampleShift,yMeasOff(:,2),'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
% measurements at their arrival times: 
scatter(tOfflineArrivalShift,yMeasOffEff(:,2),'DisplayName','noisy arrivals',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
ylabel('total solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% VS:  
subplot(4,2,7)
plot(tMinor,yClean(:,7),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
hold on; 
% measurements at their sampling times:
scatter(tOfflineSampleShift,yMeasOff(:,3),'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
% measurements at their arrival times: 
scatter(tOfflineArrivalShift,yMeasOffEff(:,3),'DisplayName','noisy arrvals',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
ylabel('volatile solids [-]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')

% Sac:  
subplot(4,2,8)
plot(tMinor,yClean(:,8),'DisplayName','clean output',...
     'LineStyle','-.', 'Color', colorPaletteHexMagma(2), 'LineWidth',1.5);
hold on; 
% measurements at their sampling times:
scatter(tOfflineSampleShift,yMeasOff(:,4),'DisplayName','noisy samples',...
        'Marker','o', 'MarkerEdgeColor',colorPaletteHexMagma(4), 'LineWidth',1); 
% measurements at their arrival times: 
scatter(tOfflineArrivalShift,yMeasOffEff(:,4),'DisplayName','noisy arrvals',...
        'Marker','*', 'MarkerEdgeColor',colorPaletteHexMagma(5), 'LineWidth',1);
ylabel('acetic acid [g/l]')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
       'LineStyle','-', 'Color', colorPaletteHexMagma(1), 'LineWidth',1.5); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow [l/h]')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')
title('Atline')

sgtitle('Clean and noisy simulation outputs (4 Online/3 Std/1 Acid)')

%% save results in struct MESS: 

MESS.tMinor = tMinor;                               % online sample times
MESS.tOfflineSampleShift = tOfflineSampleShift;     % offline sampling times 
MESS.tOfflineArrivalShift = tOfflineArrivalShift;   % offline arrival times  
MESS.x0 = x0Init;  
MESS.x0Norm = x0Norm;
MESS.xSolOn = xSolOn;   % state trajectories evaluated at diff. sampling times
MESS.xSolOff = xSolOff; 
MESS.xSolOnNorm = xSolOnNorm; 
MESS.xSolOffNorm = xSolOffNorm; 
MESS.inputMat = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = yClean;  
MESS.yCleanOn = yCleanOn;  
MESS.yCleanOff = yCleanOff;  
MESS.yMeas = yMeas; 
MESS.yMeasOn = yMeasOn; 
MESS.yMeasOff = yMeasOff;       % might include NaNs where last return time > tEnd
MESS.yMeasOffEff = yMeasOffEff; % NaNs removed
MESS.C = noiseCovMat; % accurate values from sensor data sheets
MESS.COn = noiseCovMatOn;
MESS.CSOff = noiseCovMatOff;

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
fileName = 'Messung_ADM1_R3_frac_norm_MultiRateSingleDelay.mat'; 
save(fullfile(pathToResults,fileName), 'MESS', 'params', 'TNum', 'MEASUnite')
