%% Version
% (R2022b) Update 2
% Erstelldatum: Oktober 2022
% Autor: Simon Hellmann

% create the MESS-vector for parameter identification with MLE

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
Vl = 100;           % liquid volume, aus Sörens GitHub[l]
Vg = 10;            % gas volume, aus Sörens GitHub [l]
rho = 1000;         % density of digestate [kg/l]

modParams    = [Kch4, Kco2, R, T, kla, kch, kdec, kli, kp, kpr, ph2o]; 
physParams   = [Vl, Vg, patm];
% combine all fixed parameters in a separate vector: 
pFix = [Kch4, Kco2, R, T, kla, kp, ph2o, Vl, Vg, patm, rho];
 
% inlet concentrations [GitHub Sören]
%      S_ch4, S_IC, S_IN,  S_h2o,  X_ch,   X_pr, X_li, X_bac,  X_ash, S_ch4,g, S_co2,g
xIn =  [0,     0,   0.592, 960.512,23.398, 4.75, 1.381,  0,     10,     0,      0]; % [g/l], xAshIn = 10 selbst gewählt

% times: 
tEnd = 7;   % [d] End of Simulation
tSS = 100;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.5;     % [d]
intLong = 2.5;      % [d]
intNorm = 1;        % [d]
% start with feeding right after steady state: 
ints = [0,intShort,intNorm,intLong,intNorm]';
% ... and transform intervals to absolute times:
tFeedOn = cumsum(ints);    
durationFeed = 8/60/24;             % [min], converted to [d]
tFeedOff = tFeedOn + durationFeed; 
tEvents = sort([tFeedOn;tFeedOff]); 
tGridOn = (0:1/24:tEnd)';   % online time grid every 0.5h [d]. online model outputs will be evaluated here later on
tGridOff = (0:1/2:tEnd)';   % offline time grid every 12h
tOverall = unique([tGridOn; tGridOff; tEvents]); % Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state)
nIntervals = length(tEvents);   % # intervals, during which feeding =on or =off (both constantly)
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [l/d]
feedFactors = [70,30,100,50,60]'/100; 
portions = feedFactors*feedMax; 
feedSS = 0.5*feedMax;    % steady state feed volume flow 
feedVolFlow = zeros(nIntervals,1);  % placeholder
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nStates = length(xIn);
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % placeholder 
% at times when there is feeding, introduce the input concentrations:
xInMat(idxFeedOn,:) = repmat(xIn,nFeedings,1); 
inputMat = [feedVolFlow,xInMat]; % structure required by right-hand side of ODE: [u,xIn]
inputMatSS = [feedSS,xIn]; 

% define function handle to compute the steady state:
odeFunSS = @(t,x) BMR4_AB_ode(t, x, physParams, inputMatSS, modParams);    

% determine steady state as initial value for simulation: 
tSpanSS = [0,tSS]; 
x0Soeren = [0.091,0.508,0.944,956.97,3.26,0.956,0.413,2.569,1,0.315,0.78];  % Sörens GitHub, xAsh0 selbst gewählt
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0Soeren); 

x0Init = xSS(end,:)';   % start with actual simulation in steady state
x0 = x0Init;            % save a copy

%% Solve ODE in a loop: iterative solution between all feeding points (both on & off)
% placeholders:
xSim = zeros(length(tOverall), nStates);
tSim = zeros(length(tOverall),1);   

% integrate system ODEs interval-wise (each region with feeding=on or =off
% is one interval):
tic
for cI = 1:nIntervals
    if cI == nIntervals   % last interval goes from last feeding impuls (off) till end of simulation ...
        tCurrent   = tEvents(end);
        tNext      = tEnd;
    else    % ... all other feeding intervals:
        tCurrent   = tEvents(cI);
        tNext      = tEvents(cI + 1);
    end
    
    % Get current feeding volume flow and inlet concentrations:
    inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
    modelOdeFun = @(t, x) BMR4_AB_ode(t, x, physParams, inputVector, modParams);
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval  = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % the solver would interpret t_ode as a time 
        % span otherwise in this case. but we want precisely 3 evaluations:
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % set 3 equally-spaced integration points
        [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0);
        xSim(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode has more than 2 elements. These are the evaluation points for ode15s to integrate between 
        [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0); % Returns >= 3 values
        xSim(idxTimeInterval, :) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0 = solVec(end, :); % update initial value for next interval
end
toc

%% Call adm1_r4_out to get output variables
% Evaluate xSim either in the online or offline time grid, discard the rest
idxGridOn = ismember(tOverall, tGridOn); 
idxGridOff = ismember(tOverall, tGridOff); 
xSolOn = xSim(idxGridOn, :);
xSolOff = xSim(idxGridOff, :);

% evaluate system output for online- and offline measurements
yCleanOnTemp = biogasmodell_mgl(xSolOn,pFix);
yCleanOffTemp = biogasmodell_mgl(xSolOff,pFix);
yCleanOn = yCleanOnTemp(:,1:3); 
yCleanOff = yCleanOffTemp(:,4:6); 

%% Plot artificial measurements and feedings
figure()

% gas volume flow: 
subplot(3,2,1)
plot(tGridOn,yCleanOn(:,1),'c-')%'DisplayName','biogasmodell_mgl')
ylabel('gas vol flow in L/d')
% xlabel('time [d]')
% legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
plot(tGridOn,yCleanOn(:,2),'r-','DisplayName','biogasmodell-mgl')
ylabel('p_ch4 in bar')
yyaxis right
stairs(tEvents, feedVolFlow/24,'g', ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
% legend('Location','NorthEast'); 

% pco2:
subplot(3,2,3)
plot(tGridOn,yCleanOn(:,3),'b-')%,'DisplayName','biogasmodell_mgl')
ylabel('p_co2 in bar')
% legend('Location','SouthEast'); 

% SIN:  
subplot(3,2,4)
plot(tGridOff,yCleanOff(:,1),'rx')%,'DisplayName','biogasmodell_mgl')
ylabel('inorg. nitrogen in g/L')
% xlabel('time [d]')
% legend('Location','SouthEast'); 

% TS: 
subplot(3,2,5)
plot(tGridOff,yCleanOff(:,2),'bx')%,'DisplayName','biogasmodell_mgl')
ylabel('TS')
xlabel('time [d]')
% legend('Location','SouthEast'); 

% oTS: 
subplot(3,2,6)
plot(tGridOff,yCleanOff(:,3),'cx')%,'DisplayName','biogasmodell_mgl')
ylabel('oTS')
xlabel('time [d]')
% legend('Location','SouthEast'); 

sgtitle('Simulationen mit p0 aus 2 vers. Fkts. für die Ausgangsgleichugen')

%% corrupt measurements with noise acc. to noise covariances: 
 
volFlowClean = yCleanOn(:,1); 
pCh4Clean = yCleanOn(:,2);
pCo2Clean = yCleanOn(:,3);
SINClean = yCleanOff(:,1); 
TSClean = yCleanOff(:,2); 
VSClean = yCleanOff(:,3); 

% define variance of sensors assuming zero mean (see Übersicht_Messrauschen.xlsx):
% std. deviations: 
% sigmaV = 0.08*1000; % Trommelgaszähler FBGA [m³/h] -> [L/h]
sigmaV = 0.2;       % Trommelgaszähler Labor [L/h]
sigmaCh4 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaCo2 = 0.2/100; % [Vol-%] -> [-]; für ideale Gase und p0 ungefähr 1 bar: -> [bar]
sigmaSIN = 0.13;    % NH4-N [g/L]
sigmaTS = 1.82/100; % Trockenschrank großer Tiegel [%] -> [-]
sigmaVS = 0.33/100; % [%] -> [-]

buffer = 1.5;         % conservative estimation of measurement accuracy
C = buffer*diag([sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS].^2);

nMeasOn = length(pCh4Clean);  % # of online measurements
nMeasOff = length(SINClean);  % # of offline measurements

% corrupt online measurements: 
volFlowMeas = volFlowClean + sigmaV*randn(nMeasOn,1);
pCh4Meas = pCh4Clean + sigmaCh4*randn(nMeasOn,1); 
pCo2Meas = pCo2Clean + sigmaCo2*randn(nMeasOn,1); 
% corrupt offline measurements: 
SINMeas = SINClean + sigmaSIN*randn(nMeasOff,1);
TSMeas = TSClean + sigmaTS*randn(nMeasOff,1);
VSMeas = VSClean + sigmaVS*randn(nMeasOff,1);

%% save results in struct MESS: 
MESS.tOn = tGridOn; 
MESS.tOff = tGridOff;
MESS.tSim = tSim'; 
MESS.x0 = x0Init; 
MESS.x = xSolOn'; 
MESS.xSim = xSim'; 
MESS.u = [tEvents, feedVolFlow, xInMat];    % u in [L/d]
MESS.yCleanOn = [volFlowClean, pCh4Clean, pCo2Clean]; 
MESS.yMeasOn = [volFlowMeas, pCh4Meas, pCo2Meas]; 
MESS.yCleanOff = [SINClean, TSClean, VSClean]; 
MESS.yMeasOff = [SINMeas, TSMeas, VSMeas]; 
MESS.C = C; 

save('Messung_ADM1_R4.mat', 'MESS', 'pFix')
