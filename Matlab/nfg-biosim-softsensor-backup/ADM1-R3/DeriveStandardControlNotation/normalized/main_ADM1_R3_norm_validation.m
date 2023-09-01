%% Version
% (R2022b) Update 5
% Erstelldatum: 29.08.2023
% Autor: Simon Hellmann

% run the ADM1-R3-norm in a simple scenario to validate arXiv model
% formulation

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
% load stead-state values from R3-frac implementation: 
load('Messung_ADM1_R3_frac_norm')

% renaming for easier understanding: 
s = system.Variables;   % s = systemParameters
systemInput = input.ADM1_R3.Variables;
parameters = parameters.ADM1_R3.Variables; % modelParameters

% extract system parameter values: 
Vl = s(1);  % gas volume
scaleFactor = 1; % scaling from lab to FBGA scale (liquid volume)
Vl = scaleFactor*Vl; 
Vg = s(2);  % liquid volume
Vg = 0.07*Vl; 
% p0 = s(3);  % atmospheric pressure [bar] XY: Achtung: Typo im GitHub
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

rho = 1000;        % mass density of digestate [g/l]
Mch4 = 16;      % molar mass CH4 [kg/kmol]
Mco2 = 44;      % molar mass CO2 [kg/kmol]
nac = 3/(pHUL - pHLL); 

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
aNum = [0.6555, 0.081837, 0.2245,  0.016932, 0.057375, -1,      0,      0,      0.11246,0, 0, 0, 0, 0, 0, 0,    0; 
        0.9947, 0.069636, 0.10291, 0.17456,  0.47666,   0,     -1,      0,      0.13486,0, 0, 0, 0, 0, 0, 0,    0;
        1.7651, 0.19133,  0.64716, 0.024406, 0.44695,   0,      0,     -1,      0.1621, 0, 0, 0, 0, 0, 0, 0,    0;
        26.5447,6.7367,  18.4808,  0.15056,  0.4778,    0,      0,      0,      0,      1, 0, 0, 0, 0, 0, 0,    0; 
        0,      0,        0,       0,        0,         0.18,   0.77,   0.05,  -1,      0, 0, 0, 0, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0.18,   0.77,   0.05,   0,     -1, 0, 0, 0, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0,-1, 0, 0, 0,    0;
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0,-1, 0, 0,    0; 
        0,      0,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0,-1, 0,    0;
        0,     -1,        0,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0, 0, c32,  0; 
        0,      0,       -1,       0,        0,         0,      0,      0,      0,      0, 0, 0, 0, 0, 0, 0,    c32;]';

%% inlet concentrations:
feedVolFlowSS = resultsSoeren.input(2); 
feedVolFlowSS = scaleFactor*feedVolFlowSS; 
xInPre = resultsSoeren.input(3:end)';  % obtain old, preliminary value
% beachte: in Sörens GitHub sind die Start- und Eingangswerte der Ionen 
% fehlerhafterweise verkehrtrum definiert. Er würde gerne die Reihenfolge
% der Zustände in der Petersen-Matrix umkehren (erst cat, dann an). Korrekt 
% ist es folgendermaßen:
ScatIn = xInPre(11);
SanIn = xInPre(12); 
SionIN = ScatIn - SanIn; 
% adapt inlet concentrations for slightly different state indexing in
% Simon's ADM1-R3 model (X_ash, S_ion = S_cat - S_an): 
nStates = length(xInPre); 
xIn = nan(nStates,1); % allocate memory
xIn([1:10,13:end]) = xInPre([1:10,13:end])'; % all states except Xash and Sion
% overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
xAshIn = 14; % selbst gewählt (grob abgeschätzt aus TS/oTS von Rindergülle/Maissilage)
xIn(11) = xAshIn;  
xIn(12) = SionIN; 

%% initial condition: 
x0Pre = resultsSoeren.x0';  % obtain old, preliminary value from Sören
Scat0 = x0Pre(11); 
San0 = x0Pre(12); 
Sion0 = Scat0 - San0; 
% adapt initial condition for slightly different state indexing in standard 
% control notation (X_ash, S_ion = S_cat - S_an): 
x0SS = nan(nStates,1); 
x0SS(1:10) = x0Pre(1:10);       % all states up to X_ac
% x0SS(9:10) = 2*x0SS(9:10);      % ensure enough biomass in beginning
x0SS(13:end) = x0Pre(13:end);   % all other states 
% overwrite values in positions of S_an/S_cat with those for S_ion/X_ash:
xAsh0 = 1; % assume 1 g/l initial ash concentration
x0SS(11) = xAsh0;
x0SS(12) = Sion0;   
% % since the dissociation constants KaIN and Kaco2 were increased, increase
% % the initial values of biomass to sustain a positive steady state:
% x0SS(9) = 1.5*x0SS(9);      % microbial biomass 
% x0SS(10) = 1.5*x0SS(10);    % methanogens 

%% miscellaneous parameters for simulation
% combine constant parameters in struct (index "Num" for numeric values): 
params = struct;    % allocate memory 
params.a = aNum; 
params.c = cNum; 
% time-variant parameters (hydrolysis constants):
% separation in slow and fast only relevant for R3-frac (so ignore)
% kchS = 0;     
% kchF = kch; 
thNum = [kch, 0, kpr, kli, kdec, muM, KSac, KInh3, 0]';
params.th = thNum; 

% times: 
% tEnd = 7;   % [d] End of Simulation
tSS = 300;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% % feeding intervals: 
% intShort = 0.5;     % [d]
% intLong = 2.5;      % [d]
% intMed = 1;        % [d]
% % start with feeding right after steady state and some pause: 
% ints = [intLong,intShort,intLong,intMed]';
% % ... and transform intervals to absolute times:
% cumInts = cumsum(ints);     % cumulated intervals
% tFeedOn = [cumInts(1);cumInts(3)];  % beginning times of feedings
% tFeedOff = [cumInts(2);cumInts(4)]; % end times of feedings
% feedingDurations = [intShort; intMed];  
% tEvents = sort([0;tFeedOn;tFeedOff]); 
% dt = 0.5/24;              % sample time [h], converte to [d]
% tGrid = (0:dt:tEnd)';     % time grid. The model outputs will be evaluated here later on
% tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state but with no feeding first)
% nIntervals = length(tEvents); 
% [~,idxFeedOn] = ismember(tFeedOn,tEvents); 
% feedMax = 10*24;  % max. feed volume flow [l/h] converted to [l/d]
% feedFactors = [70,30]'/100; 
% portions = feedFactors*feedMax; % [l/d]         	
% % steady state feed volume flow [l/d] should be the average of what is fed
% % during dynamic operation:
% % totalFeed = sum(feedingDurations.*portions);  
% feedVolFlow = zeros(nIntervals,1);  % allocate memory
% feedVolFlow(idxFeedOn) = portions;       

% % construct matrix of inlet concentrations at tFeedOn: 
% nFeedings = length(tFeedOn); 
% xInMat = zeros(nIntervals,nStates);         % allocate memory 
% % at times when there is feeding, introduce the input concentrations:
% xInMat(idxFeedOn,:) = repmat(xIn',nFeedings,1); 
% inputMat = [feedVolFlow,xInMat]; 

%% derive odeFun from symbolic model definition and determine steady state
% define symbolic ("S") variables (all vector are defined as column vectors)
xS = sym('x', [nStates 1]);   % states as col. vector
syms uS real             % input
xiS = sym('xi', [nStates,1]); % inlet concentrations (assumed known) 
thS = sym('th', [size(thNum)]);  % time-variant parameters (theta)
cS = sym('c', [size(cNum)]);   % known & constant time-invariant parameters 
aS = sym('a', [size(aNum)]);% petersen matrix with stoichiometric constants

dynamics = ADM1_R3_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
outputs = ADM1_R3_mgl_sym(xS,cS); 

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
x0DynR3 = xVecSS(end,:)';   % start actual simulation in steady state
% validation: check if this steady-state and Sören's value are identical! --> Yes
x0DynR3Soeren = resultsSoeren.xSS; 
deviationXSS = sum(abs(x0DynR3([1:10,13:end]) - x0DynR3Soeren([1:10,13:end])))
% x0DynR3frac = MESS.x0; 
x0 = x0DynR3;           % to be overwritten
y0 = g(x0,cNum);        % output at steady-state
%% derive normalized system equations
q = 8;                              % number of measurement signals    
xNormS = sym('xNorm', [nStates 1]); % normalized states as col. vector
syms uNorm real;                    % normalized input
xiNormS = sym('xi', [nStates,1]);   % normalized inlet concentrations 
TxS = sym('Tx', [nStates,1]);       % normalization matrix for states
TyS = sym('Ty', [q,1]);             % normalization matrix for outputs
syms Tu real                        % normalization variable for input

dynamicsNorm = ADM1_R3_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
outputsNorm = ADM1_R3_norm_mgl_sym(xNormS, cS, TxS, TyS); 

% turn into numeric function handles: 
fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 

% define numeric values for normalization with steady state: 
TxNum = x0; 
y0 = g(x0,cNum);    % steady state output
TyNum = y0; 
TuNum = feedVolFlowSS; % u0

% normalize simulation inputs:
uNorm = feedVolFlowSS./TuNum; 
x0SSNorm = x0SS./TxNum; 
xInNorm = xIn./TxNum; 

% simulate transition into steady state in normalized coordinates: 
odeFunNormSS = @(t,xNorm) fNorm(xNorm,uNorm,xInNorm,thNum,cNum,aNum,TxNum,TuNum); 
x0SSNorm = ones(size(x0SS)); 
x0SSNorm(x0SSNorm < 0.7) = 1; % S_ac, X_ch, X_ash, S_ac- sind zu klein!
[tVecNorm,xVecSSNorm] = ode15s(odeFunNormSS,tSpanSS,x0SSNorm); % for some reason, ode23s delivers the solution better than ode15s
nSamplesTillSS = size(xVecSSNorm,1); 
x0DynNorm = xVecSSNorm(end,:)'; % extract only last value as actual steady state
x0Norm = x0DynNorm;

% compute normalized system output at steady state: 
ySSNorm = gNorm(x0DynNorm, cNum, TxNum, TyNum); 
% compute normalized system output during transition into steady state: 
yVecSSNorm = nan(nSamplesTillSS,8); 
for k = 1:nSamplesTillSS
    xVecCurr = xVecSSNorm(k,:)'; 
    yVecSSNorm(k,:) = gNorm(xVecCurr, cNum, TxNum, TyNum)'; 
end

% perform de-normalization to check if normalization works properly: 
xVecSSDeNorm = repmat(TxNum',nSamplesTillSS,1).*xVecSSNorm;
x0DynR3DeNorm = xVecSSDeNorm(end,:)';
yVecSSDeNorm = repmat(TyNum',nSamplesTillSS,1).*yVecSSNorm;
ySSDeNorm = TyNum.*ySSNorm;     

%% Validation
% validation 1: compare relevant entries of states with Sören's SS:
deviationXDeNormSS = sum(abs(x0DynR3DeNorm([1:10,13:end]) - x0DynR3Soeren([1:10,13:end])))

% validation 2: compare relevant measurements with Sören's SS:
ySSSoeren = resultsSoeren.ySS';
deviationYDeNormSS = sum(abs(ySSSoeren(1:4) - ySSDeNorm(1:4)))

% validation 3: compare the interpolated trajectories of relevant
% measurements with Sören's 
yRelSoeren = resultsSoeren.y(:,[5,6,7,9]); 
tVecIntoSS = resultsSoeren.tVecIntoSS; % sample time 1 d
yVecSSDeNormInt = interp1(tVecNorm,yVecSSDeNorm,tVecIntoSS);

% re-order columns of yRelSoeren so they are in line with my ordering of
% outputs: 
yRelSoerenOrdered = [yRelSoeren(:,4),yRelSoeren(:,2),yRelSoeren(:,3),yRelSoeren(:,1)];

% XY: finde heraus, wie du mit ode15s auch mit der normierten
% Systemdarstellung dennoch robust einen Stationärwert finden kannst! 