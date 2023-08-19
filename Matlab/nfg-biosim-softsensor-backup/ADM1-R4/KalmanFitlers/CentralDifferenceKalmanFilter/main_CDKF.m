%% DAS Central Difference Kalman Filter fürs ADM1-R4

close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

x0Soeren = [0.091,0.508,0.944,956.97,3.26,0.956,0.413,2.569,1,0.315,0.78]';  % Sörens GitHub, xAsh0 selbst gewählt
x0Init = x0Soeren;      % initial state value
x0 = x0Init;            % x0 will be replaced in every iteration later
nStates = length(x0); 

% Tuning Kalman Filter
% Falsche Anfangsbedingungen (verrauscht, nur positive Konzentrationen!)
xHat = x0Soeren.*abs(randn(nStates,1)); 
% xHat = x0Soeren;
% % force negative entries to zero: 
% x_hat(x_hat < 0) = 0; 
P0 = eye(nStates);      % XY: sicherlich noch zu tunen!

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
xMinus = xHat;

tEvents = MESS.u(1,:);          % times of feeding events (on/off)
nEvents = length(tEvents);      % number of feeding events 
inputMat = MESS.u;    % [tEvents;feed vol flow;inlet concentrations]

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
    
    %% Call CDKF:
    
    % pass only relevant feedings during the measurement interval, because
    % only those are known in reality:
    idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
    tRelEvents = tEvents(idxRelEvents); 
    % Fall a: konstante Fütterung während Messintervall...
    if isempty(tRelEvents)
        % finde den Wert der letzten Fütterung (so diese existiert):
        idxCurrentFeed = find(tEvents < tSpan(1),1,'last');
        if isempty(idxCurrentFeed) 
            % Fall a1: wir sind ganz am Anfang und 
            % es gab noch gar keine Fütterung, übergib dann Nullen:
            feedInfo = zeros(1+1+9,1);
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

    [xPlus,PPlus,Kv] = centralDiffKalmanFilter(xMinus,POld,tSpan,feedInfo,yMeas,AC);
    
    % Abspeichern der Ergebnisse:
    ESTIMATES(:,k+1) = xPlus;
    COVARIANCE(:,k+1) = diag(PPlus); % P sollte Diagonalmatrix sein -> diagonale enthält alle Infos
    GAIN(:,k) = Kv;     % Kalman Gain * Innovation
    
    % Update für nächste Iteration:  
    x0 = xReal;         % real state (which isn't available in reality)
    xMinus = xPlus;     % estimated state from Kalman Filter
    POld = PPlus;       % state error covariance matrix
end

%% erzeuge aus den Zustands-Trajektorien Messwerte:
CDKFOutput = BMR4_AB_mgl_mat(ESTIMATES,AC.c);
yClean = MESS.yClean; 
feedVolFlow = inputMat(2,:);

%% 3) Plots der Ergebnisse

% plotte den Modell-Outut auf Basis des EKF und vergleiche mit echtem
% Output:
figure

% pch4: 
subplot(2,2,1)
plot(t,CDKFOutput(2,:),'r-',...
     'LineWidth',1.5,'DisplayName','CDKF-Output')
hold on; 
plot(tMeas,yClean(2,:),'b--','DisplayName','clean model output')
ylim([0.4,1])
ylabel('p_{ch4} in bar')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
       'DisplayName','feeding'); 
ylabel('feed vol flow in L/h')
legend('Location','NorthEast'); 
% legend(["$p_{ch4}$", "$p_{co2}$" "$\dot{V}_{feed}$ [l/h]"],"Interpreter", "latex",'Location', 'NorthEast')
% xlim([min(tGrid) max(tGrid)]);
% xlabel('time [d]')

% pco2:
subplot(2,2,2)
plot(t,CDKFOutput(3,:),'r-',...
     'LineWidth',1.5,'DisplayName','CDKF-Output')
hold on; 
plot(tMeas,yClean(3,:),'b--','DisplayName','clean model output')
ylabel('p_{co2} in bar')
legend('Location','SouthEast'); 

% gas volume flow: 
subplot(2,2,3)
plot(t,CDKFOutput(1,:),'r-',...
     'LineWidth',1.5,'DisplayName','CDKF-Output')
hold on; 
plot(tMeas,yClean(1,:),'b--','DisplayName','clean model output')
ylim([0,500])
ylabel('gas vol flow in L/s')
xlabel('time [d]')
legend('Location','SouthEast'); 

% SIN:  
subplot(2,2,4)
plot(t,CDKFOutput(4,:),'r-',...
     'LineWidth',1.5,'DisplayName','CDKF-Output')
hold on; 
plot(tMeas,yClean(4,:),'b--','DisplayName','clean model output')
ylabel('inorg. nitrogen in g/L')
xlabel('time [d]')
legend('Location','NorthEast'); 

sgtitle('Vergleich von CDKF und sauberem Modell-Output')

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
legend('Simulation', 'CDKF', 'CDKF +/- 1 \sigma')

%% Kalman Gains for CDKF
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
