%% DAS Unscented Kalman Filter fürs ADM1-R4

% close all
clear all
clc

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX
counterSigmaInit = 0;
counterSigmaProp = 0;
counterSigmaX = 0; 
counterX = 0; 

%% Load Measurement Data
load SimonsMessung_ADM1_R4
tMeas = MESS.t;
nMeas = length(tMeas);  % number of measurements taken

% extend time vector by one instance to account for index shift (Matlab
% only starts counting at 1, but we start with t0, x0...)
diff_t = diff(tMeas); 
dt = diff_t(1);             % sampling time
t = [tMeas;tMeas(end)+dt];  % add one time interval at end

x0 = MESS.x0; 
xRef = x0'; 
yRef = BMR4_AB_mgl_h2o(x0,AC.c);    % reference output signal (for normalization of yMeas)
nStates = length(x0); 

%% Tuning Kalman Filter
rng('default');     % fix seed for random number generation (for replicable results)
s = rng;

% Falsche Anfangsbedingungen (verrauscht, nur positive Konzentrationen!)
xHat = x0.*abs(randn()); 
% % force negative entries to zero: 
% x_hat(x_hat < 0) = 0; 
e = xHat - x0;  % initial state error
P0 = diag(e.^2);% initial state error covariance 

% Messrauschen
R = MESS.R;     % Kovarianzmatrix des Messrauschens (aus Sensordatenblättern)
kR = 10; 
R = kR*R;
% R(1,1) = 0.5e1*R(1,1); % tuning on volFlow
% R(4,4) = 0.5e1*R(4,4); % tuning on SIN
R(5,5) = 1e1*R(5,5); % tuning on TS
R(6,6) = 1e1*R(6,6); % tuning on VS
Sr = chol(R,'lower'); 

% Prozessrauschen
kQ = 1;  
Q = kQ*diag([0.0644, 0.6848, 1.6561, 348.1543*1E-3, 2.0443, 0.4208, 0.6206, 0.5350, 1.5051, 0.5701, 0.0912]); % alter Wert (für cUKF)
% Q = kQ*diag([1,1,1e-1,1e-3,1,1,1,1,1e-2,1e1,1e1]); % trial and error (für alle UKFs)
% Q = kQ*diag([0.016,  0.555,  0.563,  958.4*1E-3,    2.816,  2.654,  0.972,  2.894,  1.5, 5*0.374, 5e-2*0.948]);   % Q aus EKF (R4)
Sq = chol(Q,'lower');
% R = 10*R; 
QNorm = Q;
% Diagonalmatrix aus initialem Startwert xMinus
% XY: noch anzupassen!

%% Allocate memory:

% set up raw matrices for evaluation. all matrices which require
% initialization have 1 column more than the no of measurement instances:
MEAS = MESS.yMeas;
nSignals = size(MESS.yMeas,2);  % number of measurement signals
STATES = MESS.x;
ESTIMATESUkf = zeros(nMeas + 1,nStates);
COVARIANCEUkf = zeros(nStates,nStates,nMeas + 1);
GAIN = zeros(nStates,nSignals,nMeas);

% Initialize Kalman Filter:
% STATES(1,:) = x0;
ESTIMATESUkf(1,:) = xHat;
COVARIANCEUkf(:,:,1) = P0;
POld = P0;
SOld = chol(POld,'lower'); 
xMinus = xHat;

tEvents = MESS.u(:,1);          % times of feeding events (on/off)
nEvents = length(tEvents);      % number of feeding events 
inputMat = MESS.u;    % [tEvents,feed vol flow,inlet concentrations]

tic
% integrate across all (online) measurement intervals:
for k = 1:nMeas 
    
    tSpan = [t(k);t(k+1)]; % measurement interval; in reality = distance 
    % between old and new measurements
        
    % simulate measurement ...
    xReal = MESS.x(k,:); 
    yMeas = MESS.yMeas(k,:); 
    % ... and save values in STATES and MEAS:
%     STATES(k+1,:) = xReal;
%     MEAS(k,:) = yMeas;
    
    %% Call UKF:
    
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
            feedInfo = zeros(1,1+1+9);
        else 
            % Fall a2: gib ihm den Wert der letzten verfügbaren Fütterung. 
            % Finde dazu den für das aktuelle Intervall relevanten 
            % Fütterungs-Zeitpunkt:
            feedInfo = inputMat(idxCurrentFeed,:);
        end        
    % Fall b: Fütterung während Messung veränderlich...
    else
        % ... finde die Fütterungszeitpunkte und den letzten Zeitpunkt
        % davor (denn dessen Wert ist maßgeblich für das anfängliche
        % Fütterungsregime):
        idxLastFeed = find(tEvents < tSpan(1),1,'last');
        idxCurrentFeed = idxRelEvents;
        feedInfo = inputMat([idxLastFeed;idxCurrentFeed],:);
    end

    [xPlus,PPlus,K] = unscKalmanFiltervdMerwe(xMinus',POld,tSpan,feedInfo,yMeas,AC,R,Q);
%     [xPlus,PPlus,K] = unscKalmanFilterKolasAdditive(xMinus',POld,tSpan,feedInfo,yMeas,AC,R,Q);
%     [xPlus,PPlus,K] = unscKalmanFilterKolasAdditiveNorm(xMinus',POld,tSpan,feedInfo,yMeas,yRef,xRef,AC,R,QNorm);    
%     [xPlus,PPlus,K] = unscKalmanFilterKolasFullyAug(xMinus',POld,tSpan,feedInfo,yMeas,AC,R,Q);
%     [xPlus,PPlus,K] = unscKalmanFilterVachhani(xMinus',POld,tSpan,feedInfo,yMeas,AC,R,Q);
%     [xPlus,PPlus,K] = SRunscKalmanFilterVachhani(xMinus',SOld,tSpan,feedInfo,yMeas,AC,Sr,Sq);
%     [xPlus,PPlus] = constrainedUnscKalmanFilter(xMinus',POld,tSpan,feedInfo,yMeas,AC,R,Q);
    
    % Abspeichern der Ergebnisse:
    ESTIMATESUkf(k+1,:) = xPlus';
    COVARIANCEUkf(:,:,k+1) = PPlus; 
%     GAIN(:,:,k) = K;     % Kalman Gain
    
    % Update für nächste Iteration:  
%     x0 = xReal;         % real state (which isn't available in reality)
    xMinus = xPlus';    % estimated state from Kalman Filter
    POld = PPlus;       % state error covariance matrix
end
toc

%% erzeuge aus den Zustands-Trajektorien Messwerte:
UKFOutput = BMR4_AB_mgl_h2o_mat(ESTIMATESUkf,AC.c);
% UKFOutput = biogasmodell_mgl_mat(ESTIMATES,pFix);
yClean = MESS.yClean; 
yMeas = MESS.yMeas; 
% yCleanNew = BMR4_AB_mgl_mat(STATES,AC.c);
feedVolFlow = inputMat(:,2);

%% berechne Bestimmtheitsmaß für alle Messgrößen
q = size(MESS.yMeas,2);     % number of measurement signals
RMSSE = zeros(q,1);         % allocate memory
num2 = zeros(q,1); 
denom2 = zeros(q,1); 

% compute RMSSE for each measurement signal:
for kk = 1:q
    measurements = yMeas(:,kk); 
    
    % ignore the first value because that's only the output of x0:
    estimatedMeasurements = UKFOutput(2:end,kk);    
    
    % get RMSSE and squared numerators and denominators:
    [RMSSE(kk),num2(kk),denom2(kk)] = computeRMSSE(measurements,estimatedMeasurements); 
end

MSSE_i = num2./denom2;  % quotients of squared numerators and denominators
RMSSETot = sqrt(sum(MSSE_i)); 
aRMSSE = mean(RMSSE); % averaged root mean squared scaled error

%% 3) Plots der Ergebnisse
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plotte den Modell-Outut auf Basis des UKF und vergleiche mit echtem
% Output:
figure

% gas volume flow: 
% subplot(1,2,1)
% plot(t,UKFOutput(:,1),'r-',...
%      'LineWidth',1.5,'DisplayName','UKF-Output')
% hold on; 
% plot(tMeas,yClean(:,1),'b--','DisplayName','clean model output')
% ylim([-5,30])
% ylabel('gas vol flow in l/h')
% yyaxis right
% stairs(tEvents, feedVolFlow/24, 'Color', [0.9290 0.6940 0.1250], ...
%        'DisplayName','feeding'); 
% set(gca, "YColor", 'k')     % make right y-axis black 
% ylabel('feed vol flow in l/h')
% xlabel('time [d]')
% % legend('Location','NorthEast'); 
% legend()

% gas volume flow: 
subplot(1,2,1)
scatter(tMeas,MESS.yMeas(:,1),'DisplayName','Messwerte',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(t,UKFOutput(:,1), 'DisplayName','Unscented Kalman Filter',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(tMeas,yClean(:,1),'DisplayName','Wahrer Wert',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylim([0,10])
ylabel('Gasvolumenstrom [l/h]')
yyaxis right 
stairs(tEvents, feedVolFlow/24, 'DisplayName','Fütterung', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel({'Fütterungs-','Volumenstrom [l/h]'})
xlabel('Zeit [d]')
legend('Location','NorthEast'); 

% SIN:  
subplot(1,2,2)
scatter(tMeas,MESS.yMeas(:,4),'DisplayName','Messwerte',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(t,UKFOutput(:,4),'DisplayName','Unscented Kalman Filter',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8) 
plot(tMeas,yClean(:,4),'DisplayName','Wahrer Wert',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5)
ylabel('Anorg. Stickstoff [g/l]')
yyaxis right 
stairs(tEvents, feedVolFlow/24, 'DisplayName','Fütterung', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel({'Fütterungs-','Volumenstrom [l/h]'})
xlabel('Zeit [d]')


%%
% % pco2:
% subplot(3,2,3)
% plot(t,UKFOutput(:,3),'r-',...
%      'LineWidth',1.5,'DisplayName','UKF-Output')
% hold on; 
% plot(tMeas,yClean(:,3),'b--','DisplayName','clean model output')
% ylabel('p_{co2} in bar')
% % legend('Location','NorthEast'); 
% legend()
% 
% % SIN:  
% subplot(3,2,4)
% plot(t,UKFOutput(:,4),'r-',...
%      'LineWidth',1.5,'DisplayName','UKF-Output')
% hold on; 
% plot(tMeas,yClean(:,4),'b--','DisplayName','clean model output')
% ylabel('inorg. nitrogen in g/l')
% xlabel('time [d]')
% % legend('Location','NorthEast'); 
% legend()
% 
% % TS:  
% subplot(3,2,5)
% plot(t,UKFOutput(:,5),'r-',...
%      'LineWidth',1.5,'DisplayName','UKF-Output')
% hold on; 
% plot(tMeas,yClean(:,5),'b--','DisplayName','clean model output')
% ylabel('total solids [-]')
% xlabel('time [d]')
% % legend('Location','NorthEast'); 
% legend()
% 
% % VS:  
% subplot(3,2,6)
% plot(t,UKFOutput(:,6),'r-',...
%      'LineWidth',1.5,'DisplayName','UKF-Output')
% hold on; 
% plot(tMeas,yClean(:,6),'b--','DisplayName','clean model output')
% ylabel('volatile solids [-]')
% xlabel('time [d]')
% ylim([0,1.5])
% % legend('Location','NorthEast'); 
% legend()
% 
% sgtitle('Comparison of UKF estimates and clean model outputs')

%% make plots for CMP presentation with nice color palette:
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 

% plotte Modell-Oututs auf Basis des UKF und vergleiche mit echten Outputs:
figure()

% gas volume flow: 
subplot(2,2,1)
scatter(tMeas,yMeas(:,1),'DisplayName','noisy measurements',...
    'Marker','.','MarkerEdgeColor', colorPaletteHex(1));
hold on
plot(t,UKFOutput(:,1), 'DisplayName','UKF estimate',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.5);
% plot(t,EKFOutput(1,:)/24, 'DisplayName','EKF estimate',...
%     'LineStyle','-', 'Color', colorPaletteHex(5), 'LineWidth',0.8);
plot(tMeas,yClean(:,1), 'DisplayName','clean model outputs',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylim([0,8])
ylabel('gas vol flow [l/h]')
fontsize(gca,15,'points'); 
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
fontsize(gca,15,'points'); 
% legend('Location','NorthEast'); 

% pch4: 
subplot(2,2,2)
scatter(tMeas,yMeas(:,2),'DisplayName','noisy measurements',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(t,UKFOutput(:,2), 'DisplayName','UKF estimate',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
plot(tMeas,yClean(:,2), 'DisplayName','clean model outputs',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylim([0.4,0.7])
ylabel('p_{ch4} [bar]')
fontsize(gca,15,'points'); 
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
fontsize(gca,15,'points'); 
lgg = legend('Location','NorthEast');
lgd.FontSize = 15; 

% SIN:  
subplot(2,2,3)
scatter(tMeas,yMeas(:,4),'DisplayName','noisy measurements',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(t,UKFOutput(:,4), 'DisplayName','UKF estimate',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
plot(tMeas,yClean(:,4), 'DisplayName','clean model outputs',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
ylim([0,2])
ylabel('inorg. nitrogen [g/l]')
fontsize(gca,15,'points'); 
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
fontsize(gca,15,'points'); 

% oTS:  
subplot(2,2,4)
scatter(tMeas,yMeas(:,6),'DisplayName','noisy measurements',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(t,UKFOutput(:,6), 'DisplayName','UKF estimate',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
plot(tMeas,yClean(:,6), 'DisplayName','clean model outputs',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
% ylim([0.7,0.8])   % close up
ylim([0,1.5])
ylabel('volatile solids [-]')
fontsize(gca,15,'points'); 
% scatter(tMeas,yMeas(:,5),'DisplayName','noisy measurements',...
%     'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
% hold on
% plot(t,UKFOutput(:,5), 'DisplayName','UKF estimate',...
%     'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
% plot(tMeas,yClean(:,5), 'DisplayName','clean model outputs',...
%     'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
% % ylim([0.7,0.8])
% ylim([0,0.2])
% ylabel('total solids [-]')
% fontsize(gca,15,'points'); 
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
fontsize(gca,15,'points'); 

sgtitle('Comparison of UKF estimates and clean model outputs')

%% water and ash (responsible for TS/VS)
figure
% reshape ndarrays to vectors:
sigmaWaterArray = sqrt(COVARIANCEUkf(4,4,:)); 
sigmaAshArray = sqrt(COVARIANCEUkf(9,9,:)); 
nSamplings = size(sigmaWaterArray,3); 
sigmaWater = reshape(sigmaWaterArray,nSamplings,1);
sigmaAsh = reshape(sigmaAshArray,nSamplings,1);

% Water Sh2o:
subplot(2,1,1)
plot(tMeas,STATES(:,4)','DisplayName','true',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.2)
hold on
plot(t,ESTIMATESUkf(:,4)','DisplayName','estimated',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,4)+sigmaWater,'DisplayName','\pm 1 \sigma boundary', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,4)-sigmaWater,'HandleVisibility','off', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
ylabel('water concentration [kg/l]')
ylim([0.65,1.2])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
% xlim([0,7])
legend('location','Southeast')

% Ash XAsh:
subplot(2,1,2)
plot(tMeas,STATES(:,9),'DisplayName','true',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.2)
hold on
plot(t,ESTIMATESUkf(:,9),'DisplayName','estimated',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,9)+sigmaAsh,'DisplayName','\pm 1 \sigma boundary', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,9)-sigmaAsh,'HandleVisibility','off', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
ylabel('ash concentration [g/l]')
ylim([0,20])
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
% xlim([0,7])
% legend('location','Southeast')

sgtitle('UKF estimation of water and ash')

%% non-measurable states: 
% reshape ndarrays to vectors:
sigmaBiomassArray = sqrt(COVARIANCEUkf(8,8,:)); 
sigmaCHArray = sqrt(COVARIANCEUkf(5,5,:)); 
sigmaBiomass = reshape(sigmaBiomassArray,nSamplings,1);
sigmaCH = reshape(sigmaCHArray,nSamplings,1);

figure

% Biomasse X_bac:
subplot(2,1,1)
plot(tMeas,STATES(:,8),'DisplayName','true',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.2)
hold on
plot(t,ESTIMATESUkf(:,8),'DisplayName','estimated',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
xlabel('time [d]')
plot(t,ESTIMATESUkf(:,8)+sigmaBiomass,'DisplayName','\pm 1 \sigma boundary', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,8)-sigmaBiomass,'HandleVisibility','off', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
ylabel('biomass concentration [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
% xlim([0,7])
legend('location','Southeast')

% Carbohydrates Xch:
subplot(2,1,2)
plot(tMeas,STATES(:,5),'DisplayName','true',...
    'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.2)
hold on
plot(t,ESTIMATESUkf(:,5),'DisplayName','estimated',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,5)+sigmaCH,'DisplayName','\pm 1 \sigma boundary', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
plot(t,ESTIMATESUkf(:,5)-sigmaCH,'HandleVisibility','off', ...
    'LineStyle',':', 'Color', colorPaletteHex(5), 'LineWidth',0.8)
ylabel('carbohydrate concentration [g/l]')
yyaxis right
stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding pulses', ...
       'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
set(gca,'YColor','k')
ylabel('feed volume flow [l/h]')
ylim([0,70])
xlabel('time [d]')
% xlim([0,7])
% legend('location','Southeast')

sgtitle('Estimation of non-measurable states')

%% Kalman Gains for UKF
% figure
% subplot(321)
% stairs(t(2:end),GAIN(:,1),'b')
% ylabel('S_{ch4} [g/l]')
% grid on
% %
% subplot(323)
% stairs(t(2:end),GAIN(:,3),'r')
% ylabel('S_{IN} [g/l]')
% grid on
% %
% subplot(325)
% stairs(t(2:end),GAIN(:,9),'g')
% ylabel('X_{ash} [g/l]')
% xlabel('Zeit [h]')
% grid on
% %%%%%%%%%%%%%%%%%%%%%%%
% subplot(322)
% stairs(t(2:end),GAIN(:,5),'b')
% title('Korrektur der Zustände: K*v')
% ylabel('X_{ch} [g/l]')
% grid on
% %
% subplot(324)
% stairs(t(2:end),GAIN(:,6),'r')
% ylabel('X_{pr} [g/l]')
% grid on
% %
% subplot(326)
% stairs(t(2:end),GAIN(:,7),'g')
% ylabel('X_{li} [g/l]')
% xlabel('Zeit [h]')
% grid on
% sgtitle('Korrektur der Zustände: K*v')
