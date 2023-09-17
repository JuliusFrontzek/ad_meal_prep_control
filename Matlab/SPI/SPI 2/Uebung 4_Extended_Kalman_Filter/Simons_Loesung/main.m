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
ESTIMATES = zeros(3,nEval);
COVARIANCE = zeros(3,nEval);
GAIN = zeros(3,nEval);

STATES(:,1) = x0;
MESS(1,1) = x0(1)/x0(3);
MESS(2,1) = x0(3);
ESTIMATES(:,1) = x_hat;
COVARIANCE(:,1) = diag(P0);

t(end+1) = t(end)+dt;   % Die Verlängerung um ein Interval braucht man, 
% damit man für time = tend immer noch einen zugehörigen Eintrag in u hat. 
% Grund ist auch der Index-Shift in Matlab, da man ja eig. bei 0 anfängt,
% Matlab aber bei 1 beginnt.
% Stellgrößenverlauf (Hut-kurve rechteckig):
nEvalNew = nEval + 1;
u = zeros(1,nEvalNew);
u(round(0.5*nEvalNew):round(0.6*nEvalNew)) = 0.25;
i = 2;

POld = P0;
xMinus = x_hat;

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
    [xPlus,PPlus,Kv] = extended_kalman_filter(xMinus,u(i),yMeas,t_span,POld);
    ESTIMATES(:,i) = xPlus;
    COVARIANCE(:,i) = diag(PPlus);
    GAIN(:,i) = Kv; 
    
    % Update für nächste Iteration:
    i = i+1;  
    x0 = xReal;
    xMinus = xPlus; 
end

%% 3) Plots der Ergebnisse

figure
% Biomasse
subplot(311)
plot(t,MESS(1,:)'.*STATES(3,:)','ok')
hold on
plot(t,STATES(1,:)','k')
plot(t,ESTIMATES(1,:)','r')
plot(t,ESTIMATES(1,:)'+sqrt(COVARIANCE(1,:))','--r')
plot(t,ESTIMATES(1,:)'-sqrt(COVARIANCE(1,:))','--r')
hold off
ylabel('Biomasse m_X [g]')
xlim([0,t(end)])
title('Simulierte und geschätzte Zustände')
legend('Messung', 'Simulation', 'EKF', 'EKF +/- 1 \sigma')

% Substrat
subplot(312)
plot(t,STATES(2,:)','k')
hold on
plot(t,ESTIMATES(2,:)','r')
plot(t,ESTIMATES(2,:)'+sqrt(COVARIANCE(2,:))','--r')
plot(t,ESTIMATES(2,:)'-sqrt(COVARIANCE(2,:))','--r')
hold off
ylabel('Substrat m_S [g]')
xlim([0,t(end)])

% Volumen
subplot(313)
plot(t,STATES(3,:)','k')
hold on
plot(t,MESS(2,:)','ok')
plot(t,ESTIMATES(3,:)','r')
plot(t,ESTIMATES(3,:)'+sqrt(COVARIANCE(3,:))','--r')
plot(t,ESTIMATES(3,:)'-sqrt(COVARIANCE(3,:))','--r')
hold off
% ylim([0.9,1.1])
ylabel('Volumen [l]')
xlim([0,t(end)])
xlabel('Zeit [h]')

% Kalman Gain
figure
subplot(311)
stairs(t,GAIN(1,:),'b')
title('Korrektur der Zustände: K*v')
xlim([0,t(end)])
ylabel('Biomasse m_X [g]')
grid on
%
subplot(312)
stairs(t,GAIN(2,:),'r')
xlim([0,t(end)])
ylabel('Substrat m_S [g]')
grid on
%
subplot(313)
stairs(t,GAIN(3,:),'g')
ylabel('Substrat m_S [g]')
xlim([0,t(end)])
xlabel('Zeit [h]')
grid on