%% SPI - II - Übung 4 - DAS Kalman Filter

close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

t0 = 0;
dt = 0.5;
tend = 30;
t = t0:dt:tend;     % time vector

x0 = [10,75,5]';     % initial state value

% Tuning Kalman Filter
% x_hat = x0;         % Initialisierung mit echtem Startwert (1)-(2a)
% P0 = 1e-3*eye(3);   % für den Fall perfekt abgeschätzter Startwerte (2a)

% Falsche Anfangsbedingung
x_hat = [17 86 5.3]';   % (2b) ff.
P0 = diag([5,2,0.5]);   % (2d), selbst gewählt. Entsprechend rel. Abweichung des Anfangs-Schätzers \hat x_0 vom Anfangs-Zustand x0

%% 2) Durchführen des Prozesses
nEval = length(t); 
MESS = zeros(2,nEval);
STATES = zeros(3,nEval);
ESTIMATESEKF = zeros(3,nEval);
COVARIANCEEKF = zeros(3,nEval);
GAINEKF = zeros(3,nEval);
ESTIMATESUKF = zeros(3,nEval);
COVARIANCEUKF = zeros(3,nEval);
GAINUKF = zeros(3,nEval);

STATES(:,1) = x0;
MESS(1,1) = x0(1)/x0(3);
MESS(2,1) = x0(3);
ESTIMATESEKF(:,1) = x_hat;
COVARIANCEEKF(:,1) = diag(P0);
ESTIMATESUKF(:,1) = x_hat;
COVARIANCEUKF(:,1) = diag(P0);

%% Tuning
% Initialisierung
R = diag([1.15^2,0.25^2]);      % für Messrauschen (ist dank gegebener Sensordaten (Varianzen) fest)
% Q = diag([0.03,0.03,1]);      % für Prozessrauschen
Q = zeros(3);                   % (2a)
Q = diag([0.0527,0.3504,0.25]); % (2d) - mit Werten aus plainSimulation; bis 3% Fehler in x0 (bioprocess.m) gehen noch gut

Q = diag([0.0527,0.7,   0.25]); % (2d) - Simons beste Lösung bei 5% Ungenauigkeit
Q(1,3) = 0.03;                  % (2d) - gehört noch dazu (weil x1 und x3 am Ende stark in ihrer Unsicherheit korrelieren)

% Q = 0.01*diag([10,75,5]);       % (2d) - aus MuLö, funktioniert auch bei 5% Modell-Ungenauigkeit noch sehr gut!

% Parametersatz I
p_KF = [0.1,0.2,0.6];

% Parametersatz II
% p_KF = [0.11,0.205,0.59];

t(end+1) = t(end)+dt;   % Die Verlängerung um ein Interval braucht man, 
% damit man für time = tend immer noch einen zugehörigen Eintrag in u hat. 
% Grund ist auch der Index-Shift in Matlab, da man ja eig. bei 0 anfängt,
% Matlab aber bei 1 beginnt.
% Stellgrößenverlauf (Hut-kurve rechteckig):
nEvalNew = nEval + 1;
u = zeros(1,nEvalNew);
u(round(0.5*nEvalNew):round(0.6*nEvalNew)) = 0.25;
i = 2;

xMinusEKF = x_hat;
PMinusEKF = P0;
xMinusUKF = x_hat;
PMinusUKF = P0;

rng('default');     % fix seed for random number generation (for replicable results)
% integriere für jedes Zeitintervall separat:
for time = t0:dt:tend
    
    t_span = [time time+dt];
    
    %% Simulation und Messung
    % berechne tatsächliches x und verrauschtes y am Ende des Intervalls...
    [xReal,yMeas] = bioprocess(t_span,x0,u(i));  
    % ... und speichere diese Werte in STATES und MESS ab:
    STATES(:,i) = xReal;
    MESS(:,i) = yMeas;
    
    %% Aufruf des EKFs
    [xPlusEKF,PPlusEKF,KvEKF] = my_extended_kalman_filter(xMinusEKF,PMinusEKF,u(i),yMeas,t_span,p_KF,Q,R);
    [xPlusUKF,PPlusUKF,KvUKF] = my_UKF_additive(xMinusUKF,PMinusUKF,u(i),yMeas,t_span,p_KF,Q,R);
    [xPlusUKF,PPlusUKF] = my_UKF_fullyAugmented(xMinusUKF,PMinusUKF,u(i),yMeas,t_span,p_KF,Q,R);
%     [xPlus,PPlus,Kv] = extended_kalman_filter(xMinus,u(i),yMeas,t_span,POld);
    ESTIMATESEKF(:,i) = xPlusEKF;
    COVARIANCEEKF(:,i) = diag(PPlusEKF);
    GAINEKF(:,i) = KvEKF; 
    ESTIMATESUKF(:,i) = xPlusUKF;
    COVARIANCEUKF(:,i) = diag(PPlusUKF);
    GAINUKF(:,i) = KvUKF; 

    % Update für nächste Iteration:
    i = i+1;  
    x0 = xReal;
    xMinusEKF = xPlusEKF; 
    PMinusEKF = PPlusEKF; 
    xMinusUKF = xPlusUKF; 
    PMinusUKF = PPlusUKF;
end

%% 3) Plots der Ergebnisse

figure
% Biomasse
subplot(311)
plot(t,MESS(1,:)'.*STATES(3,:)','ok')
hold on
plot(t,STATES(1,:)','k')
plot(t,ESTIMATESEKF(1,:)','r')
plot(t,ESTIMATESUKF(1,:)','b')
plot(t,ESTIMATESEKF(1,:)'+sqrt(COVARIANCEEKF(1,:))',':r')
plot(t,ESTIMATESUKF(1,:)'+sqrt(COVARIANCEUKF(1,:))',':b')
plot(t,ESTIMATESEKF(1,:)'-sqrt(COVARIANCEEKF(1,:))','--r')
plot(t,ESTIMATESUKF(1,:)'-sqrt(COVARIANCEUKF(1,:))','--b')
hold off
ylabel('Biomasse m_X [g]')
xlim([0,t(end)])
title('Simulierte und geschätzte Zustände')
legend('Messung', 'Simulation', 'EKF', 'UKF', 'EKF +/- 1 \sigma', 'UKF +/- 1 \sigma')

% Substrat
subplot(312)
plot(t,STATES(2,:)','k')
hold on
plot(t,ESTIMATESEKF(2,:)','r')
plot(t,ESTIMATESUKF(2,:)','b')
plot(t,ESTIMATESEKF(2,:)'+sqrt(COVARIANCEEKF(2,:))',':r')
plot(t,ESTIMATESUKF(2,:)'+sqrt(COVARIANCEUKF(2,:))',':b')
plot(t,ESTIMATESEKF(2,:)'-sqrt(COVARIANCEEKF(2,:))','--r')
plot(t,ESTIMATESUKF(2,:)'-sqrt(COVARIANCEUKF(2,:))','--b')
hold off
ylabel('Substrat m_S [g]')
xlim([0,t(end)])

% Volumen
subplot(313)
plot(t,STATES(3,:)','k')
hold on
plot(t,MESS(2,:)','ok')
plot(t,ESTIMATESEKF(3,:)','r')
plot(t,ESTIMATESUKF(3,:)','b')
plot(t,ESTIMATESEKF(3,:)'+sqrt(COVARIANCEEKF(3,:))',':r')
plot(t,ESTIMATESUKF(3,:)'+sqrt(COVARIANCEUKF(3,:))',':b')
plot(t,ESTIMATESEKF(3,:)'-sqrt(COVARIANCEEKF(3,:))','--r')
plot(t,ESTIMATESUKF(3,:)'-sqrt(COVARIANCEUKF(3,:))','--b')
hold off
% ylim([0.9,1.1])
ylabel('Volumen [l]')
xlim([0,t(end)])
xlabel('Zeit [h]')

% Kalman Gains
figure
subplot(311)
stairs(t,GAINEKF(1,:),'b')
hold on
stairs(t,GAINUKF(1,:),'b','LineStyle','-.')
title('Korrektur der Zustände: K*v')
xlim([0,t(end)])
ylabel('Biomasse m_X [g]')
legend('EKF', 'UKF')
grid on
%
subplot(312)
stairs(t,GAINEKF(2,:),'r')
hold on
stairs(t,GAINUKF(2,:),'r','LineStyle','-.')
xlim([0,t(end)])
ylabel('Substrat m_S [g]')
grid on
%
subplot(313)
stairs(t,GAINEKF(3,:),'g')
hold on
stairs(t,GAINUKF(3,:),'g','LineStyle','-.')
ylabel('Substrat m_S [g]')
xlim([0,t(end)])
xlabel('Zeit [h]')
grid on