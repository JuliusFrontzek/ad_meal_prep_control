%% main file for UKF of ADM1-R4 using Matlab's System Identification Toolbox

close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

nStates = 11; 

%% 2) Load Measurement Data and initialize EKF
% load SimonsMessung_ADM1_R4 % auf Basis meiner arXiv-Modelle
load Messung_ADM1_R4 % auf Basis von Sörens/Manuels Modellen
tMeas = MESS.t;
nMeas = length(tMeas);  % number of measurements taken

% % fester Parametersatz: (Hydrolysekonstanten)
% th = [0.25, 0.2, 0.1, 0.02];   % [kch, kpr, kli, kdec] 
% % füge den zum struct AC hinzu: 
% AC.th = th; 

% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMeas); 
dt = diff_t(1);             % sampling time
t = [tMeas;tMeas(end)+dt];  % add one time interval at end

% set up raw matrices for evaluation. all matrices which require
% initialization have 1 column more than the no of measurement instances:
MEAS = MESS.yMeas;
TrueSTATES = MESS.x; 
ESTIMATES = zeros(nMeas + 1,nStates);
COVARIANCE = zeros(nStates,nStates,nMeas + 1);
GAIN = zeros(nMeas,nStates);

% Tuning Kalman Filter:
xInitialGuess = TrueSTATES(1,:); 
% Falsche Anfangsbedingungen (verrauscht, nur positive Konzentrationen!)
% xInitialGuess = x0Soeren.*abs(randn(nStates,1)); 
PInitialGuess = eye(nStates);       % XY: sicherlich noch zu tunen!

% Initialize Kalman Filter:
ESTIMATES(1,:) = xInitialGuess;     % estimated states
COVARIANCE(:,:,1) = PInitialGuess;  % state error covariance matrix P
POld = PInitialGuess;
xMinus = xInitialGuess;

tEventsFromStart = MESS.u(:,1);  % times of feeding events (on/off) without
% nIntervals = length(tEventsFromStart);
fullInputMat = MESS.u;  % [tEventsFromStart,feed vol flow,inlet concentrations]

% define function handles for state transition and measurement functions:
% xTransFun = @(xOld,tSpan,fullInputMat,params) stateTransitionFunction(xOld,tSpan,fullInputMat,params); 
% measFun = @(x,pFix) biogasmodell_mgl(x,pFix);

% construct UKF object: 
ukf = unscentedKalmanFilter(...
    @stateTransitionFunction, ...  % state transisiton function (x, ...)
    @biogasmodell_mgl, ...    % Measurement function (x,v)
    xInitialGuess);
%     xTransFun,...
%     measFun,...

ukf.Alpha = 1E-3; % default: 1E-3
% choose beta such that wc0 = nu > 0 (assuming k = 0):
lambda = nStates *((ukf.Alpha)^2 - 1); 
nu = 1; 
ukf.Beta = (ukf.Alpha)^2 - 1 - lambda/(nStates + lambda) + nu; 

%  @(x, tSpan, feedInfo, AC)stateTransitionFunction(x, tSpan, feedInfo, AC),...     % State transition function (x)
%  @(x)BMR4_AB_mgl(x,AC.c),...       

% define measurement noise covariance: 
sigmaV = 0.2;       % Trommelgaszähler [L/h], converted to [L/s]
sigmaCh4 = 0.002;   % [%]
sigmaCo2 = 0.002;   % [%]
sigmaSIN = 0.12;    % [g/L]
sigmaTS = 0.01;     % Trocknungswaage
sigmaVS = 0.0035;   
sigma = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
R = diag(sigma.^2); % Kovarianzmatrix des Messrauschens (aus Sensordatenblättern)

Q = 1E-3*diag([0.0644, 0.6848, 1.6561, 348.1543, 2.0443, 0.4208, 0.6206, 0.5350, 1.5051, 0.5701, 0.0912]); % ad-hoc Ansatz: 
% Diagonalmatrix aus initialem Startwert xMinus
% XY: noch anzupassen!
% Q = 1E-3*eye(nStates); 

ukf.MeasurementNoise = R;
ukf.ProcessNoise = Q;

% integrate across all (online) measurement intervals (like in reality):
for k = 1:nMeas 
    
    tSpan = [t(k); t(k+1)]; % measurement interval; in reality = distance 
    % between old and new measurements
        
    yMeas = MESS.yMeas(k,:);  % measurement value
    
    %% Call UKF:
    
%     % prepare vector feedInfo for the current measurement interval: 
%     % we only want to pass relevant feedings during the measurement 
%     % interval, because only those would be known in reality:
%     idxFeedEventsInTSpan = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
%     tFeedEventsInTSpan = tEvents(idxFeedEventsInTSpan);
%     % Fall a: konstante Fütterung während Messintervall...
%     % XY: überprüfe, ob interp1 nicht effizienter wäre!
%     % XY: übergib als ersten Schuss vllt. einfach alle
%     % Fütterungs-Informationen und selektiere nur noch in
%     % stateTransitionFunction!
%     if isempty(tFeedEventsInTSpan)
%         % finde den Wert der letzten Fütterung (so diese existiert):
%         idxCurrentFeed = find(tEvents < tSpan(1),1,'last');
%         if isempty(idxCurrentFeed) 
%             % Fall a1: wir sind ganz am Anfang und 
%             % es gab noch gar keine Fütterung, übergib dann Nullen (1 für
%             % tFeed, 1 für feedVolFlow, nStates für inlet concentrations):
%             feedInfo = zeros(1,1+1+nStates);
%         else 
%             % Fall a2: gib ihm den Wert der letzten verfügbaren Fütterung. 
%             % Finde dazu den für das aktuelle Intervall relevanten 
%             % Fütterungs-Zeitpunkt:
%             feedInfo = fullInputMat(idxCurrentFeed,:);
%         end        
%     % Fall b: Fütterung während Messung veränderlich...
%     else
%         % ... finde die Fütterungszeitpunkte und den letzten Zeitpunkt
%         % davor (denn dessen Wert ist maßgeblich für das anfängliche
%         % Fütterungsregime):
%         idxLastFeed = find(tEvents < tSpan(1),1,'last');
%         idxCurrentFeed = idxFeedEventsInTSpan;
%         feedInfo = fullInputMat([idxLastFeed;idxCurrentFeed],:);
%     end
    
    % do time update: 
    % übergib vereinfacht sämtliche Fütterungs-Informationen!
    % XY: argumente prüfen!
    [xMinus,PMinus] = predict(ukf,tSpan,fullInputMat,params);
    
    % do measurement update: 
    [xPlus,PPlus] = correct(ukf,yMeas,pFix);  
    
    % XY: hard reset negativer Konzentrationen:
%     xPlus(xPlus<0) = 0; 

    % alt:
    % [xPlus,PPlus,Kv] = centralDiffKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC);
    
    % Abspeichern der Ergebnisse:
    ESTIMATES(k+1,:) = xPlus;
    COVARIANCE(:,:,k+1) = PPlus; 
    
end

%% erzeuge aus den Zustands-Trajektorien synthetische Messwerte:
UKFOutput = biogasmodell_mgl_mat(ESTIMATES,pFix);
yClean = MESS.yClean; 
feedVolFlow = fullInputMat(:,2);

%% 3) Plots der Ergebnisse

% plotte den Modell-Outut auf Basis des EKF und vergleiche mit echtem
% Output:
figure

% pch4: 
subplot(3,2,1)
plot(t,UKFOutput(:,2),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,2),'b--','DisplayName','clean model output')
% ylim([0.4,1])
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding');
set(gca, "YColor", 'k')             % make right y-axis black 
ylabel('feed vol flow in L/h')
legend('Location','NorthEast'); 

% pco2:
subplot(3,2,2)
plot(t,UKFOutput(:,3),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,3),'b--','DisplayName','clean model output')
ylabel('p_{co2} in bar')
legend('Location','SouthEast'); 

% gas volume flow: 
subplot(3,2,3)
plot(t,UKFOutput(:,1),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,1),'b--','DisplayName','clean model output')
% ylim([0,500])
ylabel('gas vol flow in L/s')
xlabel('time [d]')
legend('Location','SouthEast'); 

% SIN:  
subplot(3,2,4)
plot(t,UKFOutput(:,4),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,4),'b--','DisplayName','clean model output')
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
legend('Location','SouthWest'); 

% TS:  
subplot(3,2,5)
plot(t,UKFOutput(:,5),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,5),'b--','DisplayName','clean model output')
ylabel('TS [-]')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
set(gca, "YColor", 'k')             % make right y-axis black 
ylabel('feed vol flow in L/h')
xlabel('time [d]')
legend('Location','SouthWest'); 

% VS:  
subplot(3,2,6)
plot(t,UKFOutput(:,6),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas,yClean(:,6),'b--','DisplayName','clean model output')
ylabel('VS [-]')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
set(gca, "YColor", 'k')             % make right y-axis black 
ylabel('feed vol flow in L/h')
xlabel('time [d]')
legend('Location','NorthWest'); 

sgtitle('Vergleich von UKF und sauberem Modell-Output')

%% Plotte Trajektorien markanter States: 
figure()

% S_IN:
subplot(2,2,1)
plot(t,ESTIMATES(:,3),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,3),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{IN} in g/l')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black 

% S_h2o:
subplot(2,2,2)
plot(t,ESTIMATES(:,4),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,4),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{h2o} in g/l')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

% X_ch:
subplot(2,2,3)
plot(t,ESTIMATES(:,5),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,5),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('X_{ch} in g/l')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

% S_ch4_gas:
subplot(2,2,4)
plot(t,ESTIMATES(:,10),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,10),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{ch4}^{gas} in g/l')
yyaxis right
stairs(tEventsFromStart, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('vol flow in L/h')
legend('Location','NorthEast'); 
set(gca, "YColor", 'k')     % make right y-axis black

sgtitle('Vergleich der mit UKF geschätzten und wahren Zustände')

% XY: die Auswertung der Einträge im array COVARIANCE sind noch anzupassen:
% % Biomasse X_bac:
% figure
% plot(t,STATES(8,:)','ok')
% hold on
% plot(t,ESTIMATES(8,:)','r')
% plot(t,ESTIMATES(8,:)'+sqrt(COVARIANCE(8,:))','--r')
% plot(t,ESTIMATES(8,:)'-sqrt(COVARIANCE(8,:))','--r')
% hold off
% ylabel('Biomasse X_{bac} [g/L]')
% title('Simulierte und geschätzte Biomassekonzentration')
% legend('Simulation', 'CDKF', 'CDKF +/- 1 \sigma')
% 
% %% Kalman Gains for CDKF
% figure
% subplot(321)
% stairs(t(2:end),GAIN(1,:),'b')
% ylabel('S_{ch4} [g/L]')
% grid on
% %
% subplot(323)
% stairs(t(2:end),GAIN(3,:),'r')
% ylabel('S_{IN} [g/L]')
% grid on
% %
% subplot(325)
% stairs(t(2:end),GAIN(9,:),'g')
% ylabel('X_{ash} [g/L]')
% xlabel('Zeit [h]')
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%
% subplot(322)
% stairs(t(2:end),GAIN(5,:),'b')
% title('Korrektur der Zustände: K*v')
% ylabel('X_{ch} [g/L]')
% grid on
% %
% subplot(324)
% stairs(t(2:end),GAIN(6,:),'r')
% ylabel('X_{pr} [g/L]')
% grid on
% %
% subplot(326)
% stairs(t(2:end),GAIN(7,:),'g')
% ylabel('X_{li} [g/L]')
% xlabel('Zeit [h]')
% grid on
% sgtitle('Korrektur der Zustände: K*v')
