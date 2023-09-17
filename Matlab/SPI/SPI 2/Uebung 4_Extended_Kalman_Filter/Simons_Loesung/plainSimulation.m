%% SPI - II - Übung 4 - DAS Kalman Filter

close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

t0 = 0;
dt = 0.5;
tend = 30;
t = t0:dt:tend;     % time vector [h]

x0 = [10,75,5]';     % initial state value

%% 2) Durchführen des Prozesses
nEval = length(t); 
MESS = zeros(2,nEval);
STATES = zeros(3,nEval);

STATES(:,1) = x0;
MESS(:,1) = messgleichung(x0); 

t(end+1) = t(end)+dt;   % XY: wozu braucht man das? 
% Stellgrößenverlauf (Hut-kurve rechteckig):
nEvalNew = length(t);
u = zeros(1,nEvalNew);
u(round(0.5*nEvalNew):round(0.6*nEvalNew)) = 0.25;
i = 2;

% integriere für jedes Zeitintervall separat:
for time = t0:dt:tend
    
    t_span = [time time+dt];
    
    %% Simulation und Messung
    % berechne tatsächliches x und verrauschtes y am Ende des Intervalls...
    [xReal,yMeas] = bioprocess(t_span,x0,u(i));  
    % ... und speichere diese Werte in STATES und MESS ab:
    STATES(:,i) = xReal;
    MESS(:,i) = yMeas;
    
    % Update für nächste Iteration:
    i = i+1;  
    x0 = xReal;
end

%% compute derivatives of state trajectories: 
diffStates = diff(STATES,1,2); 
diffTime = diff(t,1,2); 

dotStates = diffStates./diffTime; 
maxDotStates = max(dotStates,[],2);
sigmaV = 0.25;  % Std.-Abweichung Messung von V, siehe Aufgabenstellung
Q = diag([maxDotStates(1:2)./10;sigmaV]);   % nimm an, dass Steigungen bis auf 10% genau bestimmt werden können, siehe Skript

%% 3) Plots der Ergebnisse

figure
% Biomasse
subplot(311)
plot(t,MESS(1,:)'.*STATES(3,:)','ok')
% hold on
plot(t,STATES(1,:)','k')
% plot(t,ESTIMATES(1,:)','r')
% plot(t,ESTIMATES(1,:)'+sqrt(COVARIANCE(1,:))','--r')
% plot(t,ESTIMATES(1,:)'-sqrt(COVARIANCE(1,:))','--r')
% hold off
ylabel('Biomasse m_X [g]')
xlim([0,t(end)])
title('Simulierte Zustände')
legend('Messung')

% Substrat
subplot(312)
plot(t,STATES(2,:)','k')
% hold on
% plot(t,ESTIMATES(2,:)','r')
% plot(t,ESTIMATES(2,:)'+sqrt(COVARIANCE(2,:))','--r')
% plot(t,ESTIMATES(2,:)'-sqrt(COVARIANCE(2,:))','--r')
% hold off
ylabel('Substrat m_S [g]')
xlim([0,t(end)])

% Volumen
subplot(313)
plot(t,STATES(3,:)','k')
% hold on
% plot(t,MESS(2,:)','ok')
% plot(t,ESTIMATES(3,:)','r')
% plot(t,ESTIMATES(3,:)'+sqrt(COVARIANCE(3,:))','--r')
% plot(t,ESTIMATES(3,:)'-sqrt(COVARIANCE(3,:))','--r')
% hold off
% ylim([0.9,1.1])
ylabel('Volumen [l]')
xlim([0,t(end)])
xlabel('Zeit [h]')

% % Kalman Gain
% figure
% subplot(311)
% stairs(t,GAIN(1,:),'b')
% title('Korrektur der Zustände: K*v')
% xlim([0,t(end)])
% ylabel('Biomasse m_X [g]')
% grid on
% subplot(312)
% stairs(t,GAIN(2,:),'b')
% xlim([0,t(end)])
% ylabel('Substrat m_S [g]')
% grid on
% subplot(313)
% stairs(t,GAIN(3,:),'b')
% ylabel('Substrat m_S [g]')
% xlim([0,t(end)])
% xlabel('Zeit [h]')
% grid on