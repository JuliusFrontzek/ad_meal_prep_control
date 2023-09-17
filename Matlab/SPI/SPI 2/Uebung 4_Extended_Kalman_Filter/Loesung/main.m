%% SPI II - Übung 4 - DAS Kalman Filter

% close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

t0 = 0;
dt = 0.5;
tend = 30;
t = t0:dt:tend;

x0 = [10,75,5];

% Tuning Kalman Filter
x_hat = x0;
P0 = 1e-3*eye(3);
% P0 = 0.1*diag(x0);

% Falsche Anfangsbedingung
x_hat = [17 86 5.3];
% P0 = eye(3);
% P0(1,1) = 20;
% P0(2,2) = 20;
% P0(3,3) = 0.5;
% P0(1,2) = 30;
% P0(2,1) = 30;

%% 2) Durchführen des Prozesses

MESS = zeros(2,length(t));
STATES = zeros(3,length(t));
ESTIMATES = zeros(3,length(t));
COVARIANCE = zeros(3,length(t));
GAIN = zeros(3,length(t));

STATES(:,1) = x0;
MESS(1,1) = x0(1)/x0(3);
MESS(2,1) = x0(3);
ESTIMATES(:,1) = x_hat;
COVARIANCE(:,1) = diag(P0);

t(end+1) = t(end)+dt;
u = zeros(1,length(t));
u(round(0.5*length(t)):round(0.6*length(t))) = 0.25;
i = 2;

for time = t0:dt:tend
    
    t_span = [time time+dt];
    
    %% Simulation und Messung
    [x_real,y_noise] = bioprocess(t_span,x0,u(i));
    STATES(:,i) = x_real;
    MESS(:,i) = y_noise;
    
    %% Aufruf des EKFs
    [x_hat,P0,Kv] = extended_kalman_filter(x_hat,y_noise,t_span,u(i),P0);
    ESTIMATES(:,i) = x_hat;
    COVARIANCE(:,i) = diag(P0);
    GAIN(:,i) = Kv;
    
    i = i+1;  
    x0 = x_real;
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
subplot(312)
stairs(t,GAIN(2,:),'b')
xlim([0,t(end)])
ylabel('Substrat m_S [g]')
grid on
subplot(313)
stairs(t,GAIN(3,:),'b')
ylabel('Substrat m_S [g]')
xlim([0,t(end)])
xlabel('Zeit [h]')
grid on