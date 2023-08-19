% simulate a standardized feeding cycle for the BMR4+AB

clc
clear
close

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
p0 = 1.0133;    % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [l/bar/d]

T = 273.15;     % operating temperature [K]
Vl = 100;       % liquid volume, aus Sörens GitHub[l]
Vg = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % density of digestate [g/l]

modParams    = [Kch4, Kco2, R, T, kla, kch, kdec, kli, kp, kpr, ph2o, rho]; 
physParams   = [Vl, Vg, p0];
 
% inlet concentrations [GitHub Sören]
% XY: für welches Substrat gelten die Werte von Sören? 
% XY: hier brauchen wir noch einen zweiten Satz an Werten (anderes Substrat)
%      S_ch4, S_IC, S_IN,  S_h2o,  X_ch,   X_pr, X_li, X_bac,  X_ash, S_ch4,g, S_co2,g
xIn = [0,     0,   0.592, 960.512,23.398, 4.75, 1.381,  0,     10,     0,      0]; % [g/l], xAshIn = 10 selbst gewählt
xIn1 = [0,     0,   0.592, 960.512,23.398, 4.75,1.381,  0,     10,     0,      0]; % [g/l], xAshIn = 10 selbst gewählt
xIn2 = [0,     0,   0.592, 900.512, 50,    10,  2,      0,     10,     0,      0]; % [g/l], xAshIn = 10 selbst gewählt
% xInMaize = [0,      0,   0.764, 663,    314.5,  24.2, 18.1, 0,      10,     0,      0]';    % [g/l] 

% times: 
tEnd = 14;  % [d] End of Simulation
tSS = 100;  % [d] Dauer, bis steady state als erreicht gilt (~ 1/3 Jahr) 

% feeding intervals: 
intShort = 0.5;     % [d]
intLong = 2.5;      % [d]
intNorm = 1;        % [d]
% first week: start with feeding right after steady state: 
intsWeek1 = [0,intShort,intNorm,intLong,intNorm]';
% every other week: start with normal interval instead of 0:
intsWeek2 = [intNorm,intShort,intNorm,intLong,intNorm]';
% combine weeks to one vector...
ints = [intsWeek1;intsWeek2];
% ... and transform intervals to absolute times:
tFeedOn = cumsum(ints);    
durationFeed = 8/60/24;             % [min], converted to [d]
tFeedOff = tFeedOn + durationFeed; 
tEvents = sort([tFeedOn;tFeedOff]); 
tGrid = (0:1/48:tEnd)';             % time grid every 0.5h [d]. The model outputs will be evaluated here later on
tOverall = unique([tGrid; tEvents]);% Join and sort timestamps

% construct vector of feeding volume flows at times tFeedOn (we start in
% steady state)
nIntervals = length(tEvents);
[~,idxFeedOn] = ismember(tFeedOn,tEvents); 
feedMax = 60*24;  % max. feed volume flow [l/h] converted to [1/d]
feedFactorsWeek = [70,30,100,50,60]'/100; 
feedFactors = repmat(feedFactorsWeek,2,1);
portions = feedFactors*feedMax; 
feedSS = 0.5*feedMax;    % steady state feed volume flow 
feedVolFlow = zeros(nIntervals,1);  % placeholder
feedVolFlow(idxFeedOn) = portions;       

% construct matrix of inlet concentrations at tFeedOn: 
nStates = length(xIn);
nFeedings = length(tFeedOn); 
xInMat = zeros(nIntervals,nStates);         % placeholder
idxFeedOn1 = idxFeedOn(1:nFeedings/2);      % first week
idxFeedOn2 = idxFeedOn(nFeedings/2+1:end);  % second week
% feed one separate substrate per week:
xInMat(idxFeedOn1,:) = repmat(xIn1,nFeedings/2,1); 
xInMat(idxFeedOn2,:) = repmat(xIn2,nFeedings/2,1);
% XY: das muss 
% noch angepasst werden für die zweite Woche mit anderen Substraten!
inputMat = [feedVolFlow,xInMat]; % structure required by right-hand side of ODE 
inputVectorSS = [feedSS,xIn]; 

% define function handle to determine the steady state:
odeFunSS = @(t,x) BMR4_AB_ode(t, x, physParams, inputVectorSS, modParams);    

% determine steady state as initial value for simulation: 
tSpanSS = [0,tSS]; 
x0Soeren = [0.091,0.508,0.944,956.97,3.26,0.956,0.413,2.569,1,0.315,0.78];  % Sörens GitHub, xAsh0 selbst gewählt
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0Soeren); 

x0 = xSS(end,:)';   % start in steady state
% save('x0_BMR4_AB.mat','x0')

%% Solve ODE in a loop without interp1
xSol = zeros(length(tOverall), nStates); % placeholder

% integriere die System-DGLs abschnittsweise (jeder Bereich mit
% Fütterung =on oder =off ist ein Abschnitt):
tic
for cI = 1:nIntervals
    if cI == nIntervals   % letztes Intervall ...
        tCurrent   = tEvents(end);
        tNext      = tEnd;
    else    % ... alle anderen Fütterungsintervalle:
        tCurrent   = tEvents(cI);
        tNext      = tEvents(cI + 1);
    end
    
    % Get current feeding volume flow and inlet concentrations:
    inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
    modelOdeFun = @(t, x) BMR4_AB_ode(t, x, physParams, inputVector, modParams);
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval  = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode       = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % der Solver würde t_ode hier als Zeitspanne 
        % interpretieren und seine eigenen Integrationszeitpunkte wählen.
        % Wir fordern aber genau 3 Auswertungspunkte
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % lege 3 äquidistante Zeitpunkte zur Auswertung der Integration fest
        [~, solVec] = ode15s(modelOdeFun, t_ode, x0);
        xSol(idxTimeInterval,:) = solVec([1 end], :);   
    else    % t_ode hat mehr als zwei Zeitpunkte. Das sind die Integrations-Zeitpunkte für ode15s
        [~, solVec] = ode15s(modelOdeFun, t_ode, x0); % Returns >= 3 values
        xSol(idxTimeInterval, :) = solVec;
    end
    
    % update initial value for next interval
    x0 = solVec(end, :);
end
toc

%% Call adm1_r4_out to get output variables
% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSol = xSol(idxGrid, :);
simDigestionResults = adm1_r4_out(xSol, physParams, modParams);

%% Plot results (gas production)

figure;
plot(tGrid, simDigestionResults(:,end)/24, 'DisplayName', 'q_{gas} (l/h) BMR4+AB');
hold on; grid on;
plot(tGrid, simDigestionResults(:,5), 'DisplayName', 'X_{ch} (g/l) BMR4+AB');

xlim([min(tGrid) max(tGrid)]);
xlabel('time (d)')
legend('Location', 'SouthWest');

% plot the partial pressures: 
figure()
plot(tGrid,simDigestionResults(:,nStates+1:nStates+2))
legend('p_{ch4}', 'p_{co2}')
