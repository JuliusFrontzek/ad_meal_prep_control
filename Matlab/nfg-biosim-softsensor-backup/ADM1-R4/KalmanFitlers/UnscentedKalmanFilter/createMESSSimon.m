% create the MESS-vector for UKF by using only models acc. to arXiv paper

clc
clear
close all

%% define parameters: 
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch = 0.25; 
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011;      
Kco2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
p0 = 1.0133;    % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [l/bar/d]

T = 273.15;     % operating temperature [K]
Vl = 100;       % liquid volume, aus Sörens GitHub[l]
Vg = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % density of digestate [g/l]

% combine all fixed parameters in a separate vector: 
pFix = [Kch4, Kco2, R, T, kla, kp, ph2o, Vl, Vg, p0, rho];

% physical and model parameters: 
physParams = [Vl, Vg, p0]; 
modParams = [Kch4, Kco2, R, T, kla, kch, kdec, kli, kp, kpr, ph2o]; 

% inlet concentrations [GitHub Sören]
%       S_ch4,  S_IC,   S_IN,  S_h2o,  X_ch,   X_pr, X_li,  X_bac, X_ash, S_ch4,g, S_co2,g
xInSS = [0,     0,      0.592, 960.512,23.398, 4.75, 1.381, 0,     10,    0,       0]; % [g/l], xAshIn = 10 selbst gewählt
xInDyn =[0,     0,      0.592, 960.512,23.398, 4.75, 1.381, 0,     10,    0,       0]; % [g/l], xAshIn = 10 selbst gewählt

%% combine model parameters for compact control notation: 
c1 = 1/Vl;  
c2 = kla; 
c3 = kla*Kch4*R*T; 
c4 = kla*Kco2*R*T; 
c5 = kla*Vl/Vg; 
c6 = kp/p0*(R*T/16)^2;
c7 = 2*kp/p0*(R*T)^2/16/44;
c8 = kp/p0*(R*T/44)^2;
c9 = kp/p0*(R*T)/16*(2*ph2o - p0);
c10 = kp/p0*(R*T)/44*(2*ph2o - p0);
c11 = kp/p0*(ph2o - p0)*ph2o; 
c12 = R*T/16;
c13 = R*T/44;
c14 = rho; 
c15 = -kp/p0/Vg*(R*T/16)^2; 
c16 = -2*kp/p0/Vg*(R*T)^2/16/44; 
c17 = -kp/p0/Vg*(R*T/44)^2; 
c18 = -kp/p0/Vg*(R*T)/16*(2*ph2o - p0);
c19 = -kp/p0/Vg*(R*T)/44*(2*ph2o - p0);
c20 = -kla*Vl/Vg*Kch4*R*T - kp/p0/Vg*(ph2o - p0)*ph2o; 
c21 = -kla*Vl/Vg*Kco2*R*T - kp/p0/Vg*(ph2o - p0)*ph2o; 
c22 = Vl/Vg; 
c23 = -kp/p0/Vg*(ph2o - p0)*ph2o; 
c = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23];

a = [0.2482,  0.6809,   0.0207,     0.0456,    -1,      0,      0,       0.1372,    0,      0,      0; 
     0.3221,  0.7954,   0.1689,     0.4588,     0,     -1,      0,       0.1723,    0,      0,      0; 
     0.6393,  0.5817,   0.0344,     0.4152,     0,      0,     -1,       0.2286,    0,      0,      0; 
     0,       0,        0,          0,          0.18,   0.77,   0.05,   -1,         0,      0,      0;
    -1,       0,        0,          0,          0,      0,      0,       0,         0,      c22,    0; 
     0,      -1,        0,          0,          0,      0,      0,       0,         0,      0,      c22]';

% combine constant parameters in struct: 
AC.a = a; 
% change unit of mass density: g/l --> kg/l
c(14) = c(14)/1000; 
AC.c = c;
AC.th = [kch, kpr, kli, kdec];

% times: 
tEnd = 7;   % [d] End of Simulation
tSS = 100;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.5;     % [d]
intLong = 2.5;      % [d]
intNorm = 1;        % [d]
% start with feeding right after steady state: 
ints = [intNorm,intShort,intNorm,intLong,intNorm]';
% ... and transform intervals to absolute times:
tFeedOn = cumsum(ints);    
durationFeed = 15/60/24;             % [min], converted to [d]
tFeedOff = tFeedOn + durationFeed; 
tEvents = sort([0;tFeedOn;tFeedOff]); 
dt = 0.5/24;                % sample time [h], converted to [d]
tGrid = (0:dt:tEnd)';       % time grid. The model outputs will be evaluated here later on
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we'll start in
% steady state)
nIntervals = length(tEvents);
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactors = [70,30,100,50,60]'/100; 
portions = feedFactors*feedMax; 
% steady state feed volume flow [l/d] should be the average of what is fed
% during dynamic operation:
totalFeed = sum(durationFeed.*portions);  
feedVolFlow = zeros(nIntervals,1);  % placeholder
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nStates = length(xInSS);
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % placeholder 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xInDyn,nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; % structure required by right-hand side of ODE: u + xIn
inputVectorSS = [totalFeed/tEnd,xInSS]; 

% define function handle to determine the steady state:
odeFunSS = @(t,x) BMR4_AB_h2o_ode(t,x,inputVectorSS,AC); 

% determine steady state as initial value for simulation: 
tSpanSS = [0,tSS]; 
x0Soeren = [0.091,0.508,0.944,956.97,3.26,0.956,0.413,2.569,1,0.315,0.78];  % Sörens GitHub, xAsh0 selbst gewählt
x0Soeren(4) = x0Soeren(4)/1000; % unit change for water g/l --> kg/l 
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0Soeren); 

x0 = xSS(end,:)';   % start in steady state  

%% Solve ODE in a loop: iterative solution between all feeding points (both on & off)
xSim = zeros(length(tOverall), nStates);% placeholder
tSim = zeros(length(tOverall),1);       % placeholder

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
    modelOdeFun = @(t, x) BMR4_AB_h2o_ode(t,x,inputVector,AC);                
        
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
    end % if
    
    % update initial values for next interval
    x0 = solVec(end, :);
end % for
toc

%% Call adm1_r4_out to get output variables
% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSol = xSim(idxGrid,:);
yClean = BMR4_AB_mgl_h2o_mat(xSol,AC.c);

%% add noise to measurements acc to noise covariances: 
volFlowClean = yClean(:,1);
pCh4Clean = yClean(:,2); 
pCo2Clean = yClean(:,3);
SINClean = yClean(:,4); 
TSClean = yClean(:,5); 
VSClean = yClean(:,6); 

% define std. deviations of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2;       % Trommelgaszähler Labor [L/h]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.12;    % NH4-N [g/L]
sigmaTS = 1.73/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.31/100; % [%] -> [-]

buffer = 1.5;         % conservative estimation of measurement accuracy
C = buffer*diag([sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS].^2);

nMeas = length(pCh4Clean);  % number of measurements
rng('default');     % fix seed for random number generation (for replicable results)

pCh4Meas = pCh4Clean + sigmaCh4*randn(nMeas,1); 
pCo2Meas = pCo2Clean + sigmaCo2*randn(nMeas,1); 
volFlowMeas = volFlowClean + sigmaV*randn(nMeas,1);
SINMeas = SINClean + sigmaSIN*randn(nMeas,1);
TSMeas = TSClean + sigmaTS*randn(nMeas,1);
VSMeas = VSClean + sigmaVS*randn(nMeas,1);

%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = xSS(end,:); 
MESS.x = xSol; 
% MESS.xSim = [tSim,xSim]'; 
MESS.u = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yClean = [volFlowClean, pCh4Clean, pCo2Clean, SINClean, TSClean, VSClean]; 
MESS.yMeas = [volFlowMeas, pCh4Meas, pCo2Meas, SINMeas, TSMeas, VSMeas]; 
MESS.R = C; 

save('SimonsMessung_ADM1_R4.mat', 'MESS', 'AC', 'pFix')

%% Plot clean and noisy measurements (partial pressures, gas volume flow, SIN and feedings)

figure()

% gas volume flow: 
subplot(3,2,1)
plot(tGrid,volFlowMeas,'b.','DisplayName','noisy')
hold on 
plot(tGrid,yClean(:,1),'r','DisplayName','clean','LineWidth',1.2)
ylabel('gas vol flow in L/h')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow in L/h')
legend('Location','NorthEast'); 

% pch4
subplot(3,2,2)
plot(tGrid,pCh4Meas,'b.','DisplayName','noisy')
hold on 
plot(tGrid,yClean(:,2),'r','DisplayName','clean','LineWidth',1.2)
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
set(gca, "YColor", 'k')     % make right y-axis black 
ylabel('feed vol flow in L/h')
legend('Location','NorthEast'); 

% pco2:
subplot(3,2,3)
plot(tGrid,pCo2Meas,'b.','DisplayName','noisy')
hold on     
plot(tGrid,yClean(:,3),'r','DisplayName','clean','LineWidth',1.2)
ylabel('p_{co2} in bar')
legend('Location','SouthEast'); 

% SIN:  
subplot(3,2,4)
plot(tGrid,SINMeas,'b.','DisplayName','noisy')
hold on     
plot(tGrid,yClean(:,4),'r','DisplayName','clean','LineWidth',1.2)
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
legend('Location','SouthEast'); 

% TS:  
subplot(3,2,5)
plot(tGrid,TSMeas,'b.',...
     'LineWidth',1.5,'DisplayName','noisy')
hold on; 
plot(tGrid,yClean(:,5),'r','DisplayName','clean','LineWidth',1.2)
ylabel('total solids [-]')
xlabel('time [d]')
% legend('Location','NorthEast'); 
legend()

% VS:  
subplot(3,2,6)
plot(tGrid,VSMeas,'b.',...
     'LineWidth',1.5,'DisplayName','noisy')
hold on; 
plot(tGrid,yClean(:,6),'r','DisplayName','clean','LineWidth',1.2)
ylabel('volatile solids [-]')
xlabel('time [d]')
% legend('Location','NorthEast'); 
legend()

sgtitle('Simulationen mit den Modellgleichungen aus arXiv')