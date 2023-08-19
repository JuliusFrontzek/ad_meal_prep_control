%% main file for UKF of ADM1-R4-ess using Matlab's System Identification Toolbox

close all
clear all
clc

%% 1) Startwerte und Simulationsparameter

nStates = 8; 

%% 2) Load Measurement Data and initialize EKF
% load SimonsMessung_ADM1_R4 % auf Basis meiner arXiv-Modelle
load Messung_ADM1_R4_ess_Abklinger % auf Basis von Sörens/Manuels Modellen
tMeas = MESS.t;
nMeas = length(tMeas);  % number of measurements taken
inputVector = MESS.u;   % [volFlow, inlet concentrations]

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

% Tuning Kalman Filter:
xInitialTrue = TrueSTATES(1,:);
xInitialGuess = 1.01 * xInitialTrue;     % just a little off
% Falsche Anfangsbedingungen (verrauscht, nur positive Konzentrationen!)
% xInitialGuess = x0Soeren.*abs(randn(nStates,1)); 

% define function handles for state transition and measurement functions:
% xTransFun = @(xOld,tSpan,fullInputMat,params) stateTransitionFunction(xOld,tSpan,fullInputMat,params); 
% measFun = @(x,pFix) biogasmodell_mgl(x,pFix);

% construct UKF object: 
ukf = unscentedKalmanFilter(...
    @stateTransitionFunAbklinger, ...  % state transisiton function (x, ...)
    @biogasmodell_mgl_ess, ...    % Measurement function (x,v)
    xInitialGuess);
%     xTransFun,...
%     measFun,...

% ukf.Alpha = 1E-5; % default: 1E-3

%  @(x, tSpan, feedInfo, AC)stateTransitionFunction(x, tSpan, feedInfo, AC),...     % State transition function (x)
%  @(x)BMR4_AB_mgl(x,AC.c),...       

% apply Schneider & Georgakis' method to construct R, Q and P:
sigma = MESS.sigma;  % measurement noise levels
kR = 1;     % Tuning parameter for increased uncertainty w.r.t. measurement errors
R = kR*diag(sigma.^2); % Kovarianzmatrix des Messrauschens (aus Sensordatenblättern)

Jp = dfdp(xInitialGuess,AC.a);      % jacobian of state ode w.r.t. parameters
Cp = 1E-6*eye(length(AC.th));       % ad-hoc choice of cov. matrix of parameter error XY: hier fehlt noch ein Bootstrap/eine Fisher-Analyse, um C besser zu bestimmen!
kQ = 1;     % Tuning parameter for increased uncertainty w.r.t. process errors
Q = kQ*Jp*Cp*Jp'; 

PInitialGuess = diag((xInitialGuess - xInitialTrue).^2); 

% Initialize Kalman Filter:
ESTIMATES(1,:) = xInitialGuess;     % estimated states
COVARIANCE(:,:,1) = PInitialGuess;  % state error covariance matrix P

% Q = 1E-3*diag([0.0644, 0.6848, 2.0443, 0.4208, 0.6206, 0.5350, 0.5701, 0.0912]); % ad-hoc Ansatz: 
% Diagonalmatrix aus initialem Startwert xMinus
% XY: noch anzupassen!
% Q = 1E-3*eye(nStates); 
% Q = zeros(nStates);

ukf.MeasurementNoise = R;
ukf.ProcessNoise = Q;

% global counterX counterY
% counterX = 0; 
% counterY = 0; 

% integrate across all (online) measurement intervals (like in reality):
for k = 1:nMeas 
    
    tSpan = [t(k); t(k+1)]; % measurement interval; in reality = distance 
    % between old and new measurements
        
    yMeas = MEAS(k,:);  % measurement value
    
    %% Call UKF:
       
    % do time update: 
    [xMinus,PMinus] = predict(ukf,tSpan,inputVector,params);
    
    % do measurement update: 
    [xPlus,PPlus] = correct(ukf,yMeas,pFix);  
    
    % Abspeichern der Ergebnisse:
    ESTIMATES(k+1,:) = xPlus;
    COVARIANCE(:,:,k+1) = PPlus; 
    
end

%% erzeuge aus den Zustands-Trajektorien synthetische Messwerte:
UKFOutput = biogasmodell_mgl_ess_mat(ESTIMATES,pFix);
yClean = MESS.yClean; 

%% 3) Plots der Ergebnisse

% plotte den Modell-Outut auf Basis des EKF und vergleiche mit echtem
% Output:
figure

% gas volume flow: 
subplot(3,1,1)
plot(t(51:end),UKFOutput(51:end,1),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas(51:end),yClean(51:end,1),'b--','DisplayName','clean model output')
% ylim([0,500])
ylabel('gas vol flow in L/s')

% pch4: 
subplot(3,1,2)
plot(t(51:end),UKFOutput(51:end,2),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas(51:end),yClean(51:end,2),'b--','DisplayName','clean model output')
% ylim([0.4,1])
ylabel('p_{ch4} in bar')

% pco2:
subplot(3,1,3)
plot(t(51:end),UKFOutput(51:end,3),'r-',...
     'LineWidth',1.5,'DisplayName','UKF-Output')
hold on; 
plot(tMeas(51:end),yClean(51:end,3),'b--','DisplayName','clean model output')
ylabel('p_{co2} in bar')

sgtitle('Vergleich von UKF und sauberem Modell-Output')

%% Plotte Trajektorien markanter States: 
figure()

% X_ch:
subplot(2,1,1)
plot(t,ESTIMATES(:,3),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,3),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('X_{ch} in g/l')

% S_ch4_gas:
subplot(2,1,2)
plot(t,ESTIMATES(:,7),'r-',...
     'LineWidth',1.5,'DisplayName','estimate')
hold on; 
plot(tMeas,TrueSTATES(:,7),'b--','DisplayName','true state')
% ylim([0.4,0.7])
ylabel('S_{ch4}^{gas} in g/l')

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
