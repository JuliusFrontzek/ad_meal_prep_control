%% Version
% (R2022b) Update 2
% Erstelldatum: November 2022
% Autor: Simon Hellmann

%% DAS Kalman Filter fürs ADM1-R4-frac

close all
clear all
clc

global counter
counter = 0; 

%% 1) Startwerte und Simulationsparameter

% Quelle: Sörens GitHub. Selbst angepasst: xAsh0, XCH je zu 50% auf XCHFast
% und XCHSlow aufgeteilt
% XY: warum ausgerechnet dieser Startwert? Können wir nicht im steady-state
% starten? 
x0Soeren = [0.091, 0.508, 0.944, 956.97, 0.5*3.26, 0.5*3.26, 0.956, 0.413, 2.569, 1, 0.315, 0.78]';  
x0SS = [0.016,0.555,0.563,958.4,1.263,2.816,2.654,0.972,2.894,10,0.374,0.948];
x0Init = x0SS;      % initial state value
x0 = x0Init;        % x0 will be replaced in every iteration later
nStates = length(x0); 

% Tuning Kalman Filter
% Falsche Anfangsbedingungen (verrauscht, nur positive Konzentrationen!)
xHat = x0Soeren.*abs(randn(nStates,1)); 
% xHat = x0Soeren;
% % force negative entries to zero: 
% x_hat(x_hat < 0) = 0; 

P0 = eye(nStates); % XY: sicher besseres Tuning möglich

%% 2) Load Measurement Data and initialize EKF
load Messung_ADM1_R4
tMeas = MESS.t;
nMeas = length(tMeas);  % number of measurements taken

% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMeas); 
dt = diff_t(1);             % sampling time
t = [tMeas,tMeas(end)+dt];  % add one time interval at end

% set up raw matrices for evaluation. all matrices which require
% initialization have 1 column more than the no of measurement instances:
MEAS = MESS.yMeas;
STATES = zeros(nStates,nMeas + 1);
ESTIMATES = zeros(nStates,nMeas + 1);
COVARIANCE = zeros(nStates,nMeas + 1);
GAIN = zeros(nStates,nMeas);

% Initialize Kalman Filter:
STATES(:,1) = x0;
ESTIMATES(:,1) = xHat;
COVARIANCE(:,1) = diag(P0);
POld = P0;
R = MESS.R; 
xMinus = xHat;

tEvents = MESS.u(1,:);          % times of feeding events (on/off)
nEvents = length(tEvents);      % number of feeding events 
inputMat = MESS.u;    % [tEvents;feed vol flow;inlet concentrations]

tic
% integrate across all (online) measurement intervals (like in reality):
for k = 1:nMeas 
    
    tSpan = [t(k) t(k+1)]; % measurement interval; in reality = distance 
    % between old and new measurements
        
    % simulate measurement ...
    xReal = MESS.x(:,k); 
    yMeas = MESS.yMeas(:,k); 
    % ... and save values in STATES and MEAS:
    STATES(:,k+1) = xReal;
    MEAS(:,k) = yMeas;
    
    %% Call EKF
    % pass only relevant feedings during the measurement interval, because
    % only those would be known in reality:
    idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
    tRelEvents = tEvents(idxRelEvents); % Auswertung anhand Index
    
    % Fall a: konstante Fütterung während Messintervall...
    if isempty(tRelEvents)
        
        % XY: statt der folgenden if-else-Abfrage könntest du hier 
        % vermutlich auch interp1 mit extrapolation 0 anwenden!
         
        % finde den Wert der letzten Fütterung (so diese existiert):
        idxCurrentFeed = find(tEvents < tSpan(1),1,'last');
        if isempty(idxCurrentFeed) 
            % Fall a1: wir sind ganz am Anfang und 
            % es gab noch gar keine Fütterung, übergib dann Nullen:
            feedInfo = zeros(1+1+nStates-2,1);
            % 1 for tFeed, 1 for VFeed, one X_in for all but the two gas
            % states
        else 
            % Fall a2: gib ihm den Wert der letzten verfügbaren Fütterung. 
            % Finde dazu den für das aktuelle Intervall relevanten 
            % Fütterungs-Zeitpunkt:
            feedInfo = inputMat(:,idxCurrentFeed);
        end  

    % Fall b: Fütterung während Messung veränderlich...
    else
        % ... finde die Fütterungszeitpunkte und den letzten Zeitpunkt
        % davor (denn dessen Wert ist maßgeblich für das anfängliche
        % Fütterungsregime):
        idxLastFeed = find(tEvents < tSpan(1),1,'last');
        idxCurrentFeed = idxRelEvents;
        feedInfo = inputMat(:,[idxLastFeed,idxCurrentFeed]);
    end

    [xPlus,PPlus,Kv] = extendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC,R); % standard EKF
%     [xPlus,PPlus,Kv] = constrainedExtendedKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC,R); % constrained EKF
    
    % Abspeichern der Ergebnisse:
    ESTIMATES(:,k+1) = xPlus;
    COVARIANCE(:,k+1) = diag(PPlus); % P sollte Diagonalmatrix sein -> diagonale enthält alle Infos
    GAIN(:,k) = Kv;     % Kalman Gain * Innovation
    
    % Update für nächste Iteration:  
    x0 = xReal;         % real state (which isn't available in reality)
    xMinus = xPlus;     % estimated state from Kalman Filter
    POld = PPlus;       % state error covariance matrix
end
toc
%% erzeuge aus den Zustands-Trajektorien Messwerte:
EKFOutput = BMR4_AB_frac_mgl_mat(ESTIMATES,AC.c);
yClean = MESS.yClean; 
feedVolFlow = inputMat(2,:);

%% berechne Bestimmtheitsmaß für alle Messgrößen
q = size(MESS.yMeas,1);     % number of measurement signals
RMSSE = zeros(q,1);         % allocate memory
num2 = zeros(q,1); 
denom2 = zeros(q,1); 

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = MESS.yMeas(kk,:); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = EKFOutput(kk,2:end);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk),num2(kk),denom2(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

MSSE_i = num2./denom2;  % quotients of squared numerators and denominators
RMSSETot = sqrt(sum(MSSE_i)); 

% note: mean of RMSSE is probably better suited for overall comparison than
% RMSSETot!

%% 3) Plots der Ergebnisse

% plotte den Modell-Outut auf Basis des EKF und vergleiche mit echtem
% Output:
figure

% gas volume flow: 
subplot(3,2,1)
plot(t,EKFOutput(1,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(1,:),'b--','DisplayName','clean model output')
ylim([0,50])
ylabel('gas vol flow in L/h')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
set(gca, "YColor", 'k')     % make right y-axis black 
xlabel('time [d]')
legend('Location','NorthEast'); 

% pch4: 
subplot(3,2,2)
plot(t,EKFOutput(2,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(2,:),'b--','DisplayName','clean model output')
ylim([0.4,0.7])
ylabel('p_{ch4} in bar')

% pco2:
subplot(3,2,3)
plot(t,EKFOutput(3,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(3,:),'b--','DisplayName','clean model output')
ylim([0.3,0.7])
ylabel('p_{co2} in bar')
legend('Location','NorthEast'); 

% SIN:  
subplot(3,2,4)
plot(t,EKFOutput(4,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(4,:),'b--','DisplayName','clean model output')
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
legend('Location','NorthEast'); 

% TS:  
subplot(3,2,5)
plot(t,EKFOutput(5,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(5,:),'b--','DisplayName','clean model output')
ylabel('total solids [-]')
xlabel('time [d]')
legend('Location','NorthEast'); 

% VS:  
subplot(3,2,6)
plot(t,EKFOutput(6,:),'r-',...
     'LineWidth',1.5,'DisplayName','EKF-Output')
hold on; 
plot(tMeas,yClean(6,:),'b--','DisplayName','clean model output')
ylabel('volatile solids [-]')
xlabel('time [d]')
ylim([0.7,0.8])
legend('Location','NorthEast'); 

sgtitle('Vergleich von EKF und ungestörtem Modell-Output')

% Biomasse X_bac:
figure
plot(t,STATES(8,:)','ok')
hold on
plot(t,ESTIMATES(8,:)','r')
plot(t,ESTIMATES(8,:)'+sqrt(COVARIANCE(8,:))','--r')
plot(t,ESTIMATES(8,:)'-sqrt(COVARIANCE(8,:))','--r')
hold off
ylabel('Biomasse X_{bac} [g/L]')
title('Simulierte und geschätzte Biomassekonzentration')
legend('Simulation', 'EKF', 'EKF +/- 1 \sigma')

%% Plotte Trajektorien markanter States: 
trueStates = MESS.x; 

figure()

% S_IN:
subplot(2,2,1)
plot(t,ESTIMATES(3,:),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,trueStates(3,:),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{IN} in g/l')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black 

% S_h2o:
subplot(2,2,2)
plot(t,ESTIMATES(4,:),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,trueStates(4,:),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{h2o} in g/l')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

% X_ch_fast:
subplot(2,2,3)
plot(t,ESTIMATES(5,:),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,trueStates(5,:),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('X_{ch,fast} in g/l')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

% S_ch4_gas:
subplot(2,2,4)
plot(t,ESTIMATES(11,:),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,trueStates(11,:),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{ch4}^{gas} in g/l')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

sgtitle('Vergleich der mit EKF geschätzten und wahren Zustände')

%% Kalman Gains for EKF
figure
subplot(321)
stairs(t(2:end),GAIN(1,:),'b')
ylabel('S_{ch4} [g/L]')
grid on
%
subplot(323)
stairs(t(2:end),GAIN(3,:),'r')
ylabel('S_{IN} [g/L]')
grid on
%
subplot(325)
stairs(t(2:end),GAIN(9,:),'g')
ylabel('X_{ash} [g/L]')
xlabel('Zeit [h]')
grid on
%%%%%%%%%%%%%%%%%%%%%%%
subplot(322)
stairs(t(2:end),GAIN(5,:),'b')
title('Korrektur der Zustände: K*v')
ylabel('X_{ch} [g/L]')
grid on
%
subplot(324)
stairs(t(2:end),GAIN(6,:),'r')
ylabel('X_{pr} [g/L]')
grid on
%
subplot(326)
stairs(t(2:end),GAIN(7,:),'g')
ylabel('X_{li} [g/L]')
xlabel('Zeit [h]')
grid on
sgtitle('Korrektur der Zustände: K*v')

%% create another subplot with only gas volume flow and pch4: 
% tL = tiledlayout(1,1,'Padding','tight');
% tL.Units = 'centimeters';
% tL.OuterPosition = [0.25 0.25 16 9];
% nexttile;

% myFig = figure('Name','First Structure');
figure
subplot(1,2,1)
plot(t,EKFOutput(1,:),'r-',...
     'LineWidth',1.5,'DisplayName','reconstructed with EKF')
hold on; 
plot(tMeas,yClean(1,:),'b--','DisplayName','clean model output')
ylim([0,50])
ylabel('gas volume flow in L/h')
xlabel('time [d]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
xlim([0,round(t(end))])
ylabel('feed volume flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black 
% 
subplot(1,2,2)
plot(t,EKFOutput(2,:),'r-',...
     'LineWidth',1.5,'DisplayName','reconstructed with EKF')
hold on; 
plot(tMeas,yClean(2,:),'b--','DisplayName','clean model output')
ylim([0.4,0.7])
ylabel('methane partial pressure in bar')
xlim([0,round(t(end))])
xlabel('time [d]')
legend('Location','NorthEast'); 

%% save as .png:
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0.25 0.25 32 12]);
print(gcf,'-dpng','EKF_simulation.png');
