% create the MESS-vector for Kalman Filtering based on Sörens/Manuels model 

clc
clear
close all

%% define parameters: 
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch = 0.25; 
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 
k = [kch, kpr, kli, kdec];

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011;      
Kco2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
patm = 1.0133;  % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [l/bar/d]

T = 273.15;     % operating temperature [K]
Vl = 100;       % liquid volume, aus Sörens GitHub[l]
Vg = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % density of digestate [g/l]

modParams    = [Kch4, Kco2, R, T, kla, kch, kdec, kli, kp, kpr, ph2o]; 
physParams   = [Vl, Vg, patm];
% combine all fixed parameters in a separate vector: 
pFix = [Kch4, Kco2, R, T, kla, kp, ph2o, Vl, Vg, patm, rho];

% combine Sörens params in struct: 
params.phyParams = physParams; 
params.modParams = modParams; 

% inlet concentrations [GitHub Sören]
%      S_ch4, S_IC, S_IN,  S_h2o,  X_ch,   X_pr, X_li, X_bac,  X_ash, S_ch4,g, S_co2,g
xIn = [0,     0,   0.592, 960.512,23.398, 4.75, 1.381,  0,     5,     0,      0]; % [g/l], xAshIn = 5 selbst gewählt

%% combine model parameters for compact control notation: 
% c1 = 1/Vl; 
% c2 = Vg; 
% c3 = kla; 
% c4 = kla*Kch4*R*T; 
% c5 = kla*Kco2*R*T; 
% c6 = kla*Vl/Vg; 
% c7 = kla*Vl/Vg*Kch4*R*T; 
% c8 = kla*Vl/Vg*Kco2*R*T; 
% c9 = kp/patm*(R*T/16)^2;
% c10 = 2*kp/patm*(R*T)^2/16/44;
% c11 = kp/patm*(R*T/44)^2;
% c12 = kp/patm*R*T/16*(2*ph2o - patm); 
% c13 = kp/patm*R*T/44*(2*ph2o - patm); 
% c14 = kp/patm*(ph2o - patm)*ph2o; 
% c15 = R*T/16;
% c16 = R*T/44;
% c17 = rho; 
% c18 = -c9/c2; 
% c19 = -c10/c2; 
% c20 = -c11/c2; 
% c21 = -c12/c2; 
% c22 = -c13/c2; 
% c23 = -(c7 + c14/c2); 
% c24 = -(c8 + c14/c2); 
% c = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24];
% 
% a = [0.2482,  0.6809,  -0.0207,    -0.0456,    -1,      0,      0,       0.137,     0,      0,          0; 
%      0.3221,  0.7954,   0.1689,    -0.4588,     0,     -1,      0,       0.1723,    0,      0,          0; 
%      0.6393,  0.5817,  -0.0344,    -0.4152,     0,      0,     -1,       0.2286,    0,      0,          0; 
%      0,       0,        0,          0,          0.19,   0.77,   0.04,   -1,         0,      0,          0;
%     -1,       0,        0,          0,          0,      0,      0,       0,         0,      Vl/Vg,      0; 
%      0,      -1,        0,          0,          0,      0,      0,       0,         0,      0,          Vl/Vg]';

% % combine constant parameters in struct: 
% AC.a = a; 
% AC.c = c; 
% AC.th = [kch, kpr, kli, kdec];

% times: 
t0 = 0;     % starting time
tEnd = 7;   % [d] End of Simulation
tSS = 100;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
% intShort = 0.5;     % [d]
% intLong = 2;        % [d]
% intNorm = 1;        % [d]
% start with feeding right after steady state: 
% ints = [0,intShort,intNorm,intLong,intNorm]';
% ... and transform intervals to absolute times:
tFeedOn = 2.5;        % [d]
durationFeed = 10/60/24;             % [min], converted to [d]
tFeedOff = tFeedOn + durationFeed; 
tEvents = sort([tFeedOn;tFeedOff]); 
tEventsFromStart = unique([t0;tEvents]);     % feeding events plus start
dt = 15/60/24;    % sampling time [min], converted to [d]
tGrid = (t0:dt:tEnd)';             % time grid every 0.5h [d]. The model outputs will be evaluated here later on
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state)
nIntervals = length(tEventsFromStart);
[~,idxFeedOn] = ismember(tFeedOn,tEventsFromStart);  % XY: überarbeiten!
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [l/d]
% feedFactors = [70,30,100,50,60]'/100; 
feedFactor = 1;
portions = feedFactor*feedMax; 
% steady state feed volume flow [l/d] should be the average of what is fed
% during dynamic operation:
totalFeed = sum(durationFeed.*portions);  
feedVolFlow = zeros(nIntervals,1);  % placeholder
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nStates = length(xIn);
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % placeholder 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn,nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; % structure required by right-hand side of ODE: u + xIn
% the mean feed vol flow during dynamic feeding regime should be fed during
% steady state, so the dynamic feeding starts in steady state:
inputVectorSS = [totalFeed/tEnd,xIn]; 

% define function handle to determine the steady state:
odeFunSS = @(t,x) BMR4_AB_ode(t, x, physParams, inputVectorSS, modParams);    

% determine steady state as initial value for simulation: 
tSpanSS = [t0,tSS]; 
x0Soeren = [0.091,0.508,0.944,956.97,3.26,0.956,0.413,2.569,1,0.315,0.78];  % Sörens GitHub, xAsh0 selbst gewählt
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0Soeren); 

x0 = xSS(end,:);   % start in steady state 
% save('x0_BMR4_AB.mat','x0')

%% Solve ODE in a loop: iterative solution between all feeding points (both on & off)
xSim = zeros(length(tOverall), nStates);% placeholder
xSimAlt = xSim; 
tSim = zeros(length(tOverall),1);       % placeholder

% integriere die System-DGLs abschnittsweise (jeder Bereich mit
% Fütterung =on oder =off ist ein Abschnitt):
tic
for cI = 1:nIntervals
    if cI == nIntervals   % letztes Intervall geht von letztem Fütterungsimpuls (aus) bis Simulationsende ...
        tCurrent   = tEventsFromStart(end);
        tNext      = tEnd;
    else    % ... alle anderen Fütterungsintervalle:
        tCurrent   = tEventsFromStart(cI);
        tNext      = tEventsFromStart(cI + 1);
    end
    
    % Get current feeding volume flow and inlet concentrations:
    inputVector = interp1(tEventsFromStart, inputMat, tCurrent, 'nearest', 0); 
    modelOdeFun = @(t, x) BMR4_AB_ode(t, x, physParams, inputVector, modParams);
        
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval  = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode       = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % der Solver würde t_ode hier als Zeitspanne 
        % interpretieren und seine eigenen Integrationszeitpunkte wählen.
        % Wir fordern aber genau 3 Auswertungspunkte:
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % lege 3 äquidistante Zeitpunkte zur Auswertung der Integration fest
        [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0);
        xSim(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode hat mehr als zwei Zeitpunkte. Das sind die Integrations-Zeitpunkte für ode15s
        [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0); % Returns >= 3 values
        xSim(idxTimeInterval, :) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    % update initial values for next interval
    x0 = solVec(end, :);
end
toc

%% Call adm1_r4_out to get output variables
% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSol = xSim(idxGrid,:);
% simDigestionResults = BMR4_AB_mgl_mat(xSol',c);
yClean = biogasmodell_mgl_mat(xSol,pFix);

%% Plot results (partial pressures, gas volume flow, SIN and feedings)

% plot and compare the results from biogasmodell_mgl and BMR4_AB_mgl: 
figure()

% pch4: 
subplot(2,2,1)
plot(tGrid,yClean(:,1),'r','DisplayName','biogasmodell-mgl')
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
set(gca, "YColor", 'k')
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 

% pco2:
subplot(2,2,2)
plot(tGrid,yClean(:,2),'r','DisplayName','biogasmodell-mgl')
ylabel('p_{co2} in bar')
legend('Location','SouthEast'); 

% gas volume flow: 
subplot(2,2,3)
plot(tGrid,yClean(:,3),'r','DisplayName','biogasmodell-mgl')
ylabel('gas vol flow in L/s')
xlabel('time [d]')
legend('Location','NorthEast'); 

% SIN:  
subplot(2,2,4)
plot(tGrid,yClean(:,4),'r','DisplayName','biogasmodell-mgl')
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
legend('Location','SouthEast'); 

sgtitle('Simulation mit Sörens/Manuels Modell')

%% add noise to measurements acc to noise covariances: 
 
pCh4Clean = yClean(:,1); 
pCo2Clean = yClean(:,2);
volFlowClean = yClean(:,3);
SINClean = yClean(:,4); 
TSClean = yClean(:,5); 
VSClean = yClean(:,6); 

% define properties of normal distributions of sensors (assume zero mean), 
% see Übersicht_Messrauschen.xlsx):
% std. deviations: 
sigmaCh4 = 0.002;   % [%]
sigmaCo2 = 0.002;   % [%]
sigmaV = 0.2;       % Trommelgaszähler [L/h]
sigmaSIN = 0.12;    % [g/L]
sigmaTS = 0.01;     % Trocknungswaage
sigmaVS = 0.0035; 

nMeas = length(pCh4Clean);  % number of measurements
pCh4Meas = pCh4Clean + sigmaCh4*randn(nMeas,1); 
pCo2Meas = pCo2Clean + sigmaCo2*randn(nMeas,1); 
volFlowMeas = volFlowClean + sigmaV*randn(nMeas,1);
SINMeas = SINClean + sigmaSIN*randn(nMeas,1);
TSMeas = TSClean + sigmaTS*randn(nMeas,1);
VSMeas = VSClean + sigmaVS*randn(nMeas,1);

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0; 
MESS.x = xSol; 
MESS.xSim = [tSim,xSim]; 
MESS.u = [tEventsFromStart, (feedVolFlow), xInMat];    % u in [L/d]
MESS.yClean = [volFlowClean, pCh4Clean, pCo2Clean, SINClean, TSClean, VSClean]; 
MESS.yMeas = [volFlowMeas, pCh4Meas, pCo2Meas, SINMeas, TSMeas, VSMeas]; 

save('Messung_ADM1_R4.mat', 'MESS', 'pFix',  'params')
