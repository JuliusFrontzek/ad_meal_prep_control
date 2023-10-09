%% SPI - II - Übung 4 - DAS Kalman Filter

% close all
clear all
clc

global counterSigmaInit
global counterSigmaProp
global counterSigmaX
global counterX
% define those separately for fully-augmented (FA) case
global counterSigmaInitFA
global counterSigmaPropFA
global counterSigmaXFA
global counterXFA
% define separately for cUKF to check if constraints are respected:
global counterSigmaXcUKF

counterSigmaInit = 0; 
counterSigmaProp = 0; 
counterSigmaX = 0; 
counterX = 0; 
% separate for FA case:
counterSigmaInitFA = 0; 
counterSigmaPropFA = 0; 
counterSigmaXFA = 0; 
counterXFA = 0; 
counterSigmaXcUKF = 0;

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
x_hat_CKF = 1.0*x_hat;
P0 = diag([5,2,0.5]);   % (2d), selbst gewählt. Entsprechend rel. Abweichung des Anfangs-Schätzers \hat x_0 vom Anfangs-Zustand x0

nEval = length(t); 

t(end+1) = t(end) + dt;   
% alle Schätzwerte müssen bei t=0 erstmal initialisiert werden. Matlab
% beginnt aber nicht bei 0, sondern bei 1 zu zählen. Daher brauchen wir
% einen Zeitschritt mehr, um bis zu tEnd simulieren und schätzen zu können.

% Stellgrößenverlauf (Hut-kurve rechteckig):
u = zeros(1,nEval+1);
u(round(0.5*(nEval+1)):round(0.6*(nEval+1))) = 0.25;
%% 2) Durchführen des Prozesses
% allocate memory:
MESS = zeros(2,nEval);
STATES = zeros(3,nEval);
% alle Matrizen für x und P werden initialisiert, brauchen also einen
% Eintrag mehr als die wahren Zustände und die Messwerte:
ESTIMATESEKF = zeros(3,nEval+1);      	% EKF
COVARIANCEEKF = zeros(3,nEval+1);
GAINEKF = zeros(3,nEval+1);
ESTIMATESUKFAdd = zeros(3,nEval+1);     % UKF additive noise
COVARIANCEUKFAdd = zeros(3,nEval+1);
GAINUKFAdd = zeros(3,nEval+1);
ESTIMATESSRUKFAdd = zeros(3,nEval+1);   % SQ-UKF additive noise
% SRCOVARIANCESRUKFAdd = zeros(3,nEval+1);
ESTIMATEScUKFAdd = zeros(3,nEval+1);    % constrained UKF with additive noise
COVARIANCEcUKFAdd = zeros(3,nEval+1);
ESTIMATESUKFAug = zeros(3,nEval+1);     % UKF with augmented system noise
COVARIANCEUKFAug = zeros(3,nEval+1);
ESTIMATEScUKFAug = zeros(3,nEval+1);    % constrainted UKF with augmented system noise
COVARIANCEcUKFAug = zeros(3,nEval+1);
ESTIMATESUKFFullyAug = zeros(3,nEval+1);% fully augmented UKF
COVARIANCEUKFFullyAug = zeros(3,nEval+1);
ESTIMATEScUKFFullyAug = zeros(3,nEval+1);% constrained fully augmented UKF
COVARIANCEcUKFFullyAug = zeros(3,nEval+1);
ESTIMATESCKF = zeros(3,nEval+1);    % Cubature Kalman Filter
COVARIANCECKF = zeros(3,nEval+1);

% initialize all matrices:
STATES(:,1) = x0;
MESS(1,1) = x0(1)/x0(3);
MESS(2,1) = x0(3);
ESTIMATESEKF(:,1) = x_hat;
COVARIANCEEKF(:,1) = diag(P0);
ESTIMATESUKFAdd(:,1) = x_hat;
COVARIANCEUKFAdd(:,1) = diag(P0);
ESTIMATESUKF_sysID(:,1) = x_hat;
COVARIANCEUKF_sysID(:,1) = diag(P0);
ESTIMATESSRUKFAdd(:,1) = x_hat;
% SRCOVARIANCESRUKFAdd(:,1) = diag(S0);
ESTIMATEScUKFAdd(:,1) = x_hat;
COVARIANCEcUKFAdd(:,1) = diag(P0);
ESTIMATESUKFAug(:,1) = x_hat;
COVARIANCEUKFAug(:,1) = diag(P0);
ESTIMATEScUKFAug(:,1) = x_hat;
COVARIANCEcUKFAug(:,1) = diag(P0);
ESTIMATESUKFFullyAug(:,1) = x_hat;
COVARIANCEUKFFullyAug(:,1) = diag(P0);
ESTIMATEScUKFFullyAug(:,1) = x_hat;
COVARIANCEcUKFFullyAug(:,1) = diag(P0);
% ESTIMATESCKF(:,1) = x_hat_CKF;
% COVARIANCECKF(:,1) = diag(P0);

xMinusEKF = x_hat;
PMinusEKF = P0;
xMinusUKFAdd = x_hat;
PMinusUKFAdd = P0;
xMinusSRUKFAdd = x_hat; % define corr. Covariance further down!
xMinuscUKFAdd = x_hat;
PMinuscUKFAdd = P0;
xMinusUKFAug = x_hat;
PMinusUKFAug = P0;
xMinuscUKFAug = x_hat;
PMinuscUKFAug = P0;
xMinusUKFFullyAug = x_hat;
PMinusUKFFullyAug = P0;
xMinuscUKFFullyAug = x_hat;
PMinuscUKFFullyAug = P0;
xMinusCKF = x_hat_CKF;
PMinusCKF = P0;

%% Tuning
R = diag([1.15^2,0.25^2]);      % für Messrauschen (ist dank gegebener Sensordaten (Varianzen) fest)
% Q = diag([0.03,0.03,1]);      % für Prozessrauschen
% Q = zeros(3);                   % (2a)
% Q = diag([0.0527,0.3504,0.25]); % (2d) - mit Werten aus plainSimulation; bis 3% Fehler in x0 (bioprocess.m) gehen noch gut
Q = diag([0.0527,0.7,   0.25]); % (2d) - Simons beste Lösung bei 5% Ungenauigkeit
Q(1,3) = 0.03;                  % (2d) - gehört noch dazu (weil x1 und x3 am Ende stark in ihrer Unsicherheit korrelieren)
% make Q symmetric: 
Q = Q + triu(Q,1)';

% Parametersatz I for model underlying the Kalman Filters
% p_KF = [0.1,0.2,0.6];

% Parametersatz II
p_KF = [0.11,0.205,0.59];

%% CKF Tuning
% tuning matrices for CKF:
R_CKF = R;
Q_CKF = Q;

% model parameters for CKF
p_CKF = [0.11,0.205,0.59]; %[0.18,0.405,0.309]; % model parameters for CKF
testParamValue = 0;     % nur um zu sehen, ob man der Messgleichung auf zusätzliche Werte übergeben kann

% Set up Cubature Kalman Filter
% ckf = trackingCKF(@StateTransitionFcn,@MeasurementFcn,x_hat_CKF, ...
%     'HasAdditiveProcessNoise',true, 'HasAdditiveMeasurementNoise',true, ...
%     'MeasurementNoise',R_CKF, 'ProcessNoise',Q_CKF);

%% UKF acc. to Matlab Toolbox

ukf = unscentedKalmanFilter(@StateTransitionFcn,@MeasurementFcn,xMinusUKFAdd, ...
    'HasAdditiveProcessNoise',true, 'HasAdditiveMeasurementNoise',true, ...
    'MeasurementNoise',R, 'ProcessNoise',Q);

%% SR-UKF van der Merwe (2001)

% funktioniert:
S0 = chol(P0,'upper'); 
SMinusSRUKFAdd = S0;
SQ = chol(Q,'lower'); 
SR = chol(R,'lower'); 

% % Probe:
% S0 = chol(P0,'upper'); 
% SMinusSRUKFAdd = S0;
% SQ = chol(Q,'lower'); 
% SR = chol(R,'lower'); 

%% call Kalman Filters iteratively for fixed time grid
i = 2;

rng('default');     % fix seed for random number generation (for replicable results)
for time = t0:dt:tend
    
    tSpan = [time time+dt];
    
    %% Simulation und Messung
    % berechne tatsächliches x und verrauschtes y am Ende des Intervalls...
    [xReal,yMeas] = bioprocess(tSpan,x0,u(i));  
    % ... und speichere diese Werte in STATES und MESS ab:
    STATES(:,i) = xReal;
    MESS(:,i) = yMeas;
    
    %% Aufruf der EKF/UKFs
    %     [xPlus,PPlus,Kv] = extended_kalman_filter(xMinus,u(i),yMeas,t_span,POld);
    [xPlusEKF,PPlusEKF,KvEKF] = my_extended_kalman_filter(xMinusEKF,PMinusEKF,u(i),yMeas,tSpan,p_KF,Q,R);
    [xPlusUKFAdd,PPlusUKFAdd,KvUKFAdd] = my_UKF_additive(xMinusUKFAdd,PMinusUKFAdd,u(i),yMeas,tSpan,p_KF,Q,R);
    [xPlusSRUKFAdd,SPlusSRUKFAdd] = my_SR_UKF_additive(xMinusSRUKFAdd,SMinusSRUKFAdd,u(i),yMeas,tSpan,p_KF,SQ,SR);
%     [xPlusCKF,PPlusCKF] = my_CKF(ckf,u(i),yMeas,tSpan,p_CKF,testParamValue);
    [xPlusUKF_sysID,PPlusUKF_sysID] = my_UKF(ukf,u(i),yMeas,tSpan,p_KF,testParamValue);
%     [xPluscUKFAdd,PPluscUKFAdd] = my_cUKF_additive(xMinuscUKFAdd,PMinuscUKFAdd,u(i),yMeas,tSpan,p_KF,Q,R);
    [xPlusUKFAug,PPlusUKFAug] = my_UKF_augmented(xMinusUKFAug,PMinusUKFAug,u(i),yMeas,tSpan,p_KF,Q,R);
%     [xPluscUKFAug,PPluscUKFAug] = my_cUKF_augmented(xMinuscUKFAug,PMinuscUKFAug,u(i),yMeas,tSpan,p_KF,Q,R);
    [xPlusUKFFullyAug,PPlusUKFFullyAug] = my_UKF_fullyAugmented(xMinusUKFFullyAug,PMinusUKFFullyAug,u(i),yMeas,tSpan,p_KF,Q,R);
%     [xPluscUKFFullyAug,PPluscUKFFullyAug] = my_cUKF_fullyAugmented(xMinuscUKFFullyAug,PMinuscUKFFullyAug,u(i),yMeas,tSpan,p_KF,Q,R);

% packe das beides muss in separate Funktion, genau wie für ekf und ukf auch:
%     [xPredCKF,PPredCKF] = predict(ckf, u(i), tSpan, p_CKF);
%     [xCorrCKF,pCorrCKF] = correct(ckf,yMeas);
%     [xPlusCKF,PPlusCKF] = my_CKF(ckf,u(i),yMeas,tSpan,p_CKF,testParamValue);

    ESTIMATESEKF(:,i) = xPlusEKF;
    COVARIANCEEKF(:,i) = diag(PPlusEKF);
    GAINEKF(:,i) = KvEKF; 
    ESTIMATESUKFAdd(:,i) = xPlusUKFAdd;
    COVARIANCEUKFAdd(:,i) = diag(PPlusUKFAdd);
    ESTIMATESUKF_sysID(:,i) = xPlusUKF_sysID;
    COVARIANCEUKF_sysID(:,i) = diag(PPlusUKF_sysID);
    ESTIMATESSRUKFAdd(:,i) = xPlusSRUKFAdd;
%     SRCOVARIANCESRUKFAdd(:,i) = diag(SPlusSRUKFAdd);
    GAINUKFAdd(:,i) = KvUKFAdd; 
%     ESTIMATEScUKFAdd(:,i) = xPluscUKFAdd;
%     COVARIANCEcUKFAdd(:,i) = diag(PPluscUKFAdd);
    ESTIMATESUKFAug(:,i) = xPlusUKFAug;
    COVARIANCEUKFAug(:,i) = diag(PPlusUKFAug);
%     ESTIMATEScUKFAug(:,i) = xPluscUKFAug;
%     COVARIANCEcUKFAug(:,i) = diag(PPluscUKFAug);
    ESTIMATESUKFFullyAug(:,i) = xPlusUKFFullyAug;
    COVARIANCEUKFFullyAug(:,i) = diag(PPlusUKFFullyAug);
%     ESTIMATEScUKFFullyAug(:,i) = xPluscUKFFullyAug;
%     COVARIANCEcUKFFullyAug(:,i) = diag(PPluscUKFFullyAug);
%     ESTIMATESCKF(:,i) = xPlusCKF;
%     COVARIANCECKF(:,i) = diag(PPlusCKF);

    % Update für nächste Iteration:
    i = i+1;  
    x0 = xReal;
    xMinusEKF = xPlusEKF; 
    PMinusEKF = PPlusEKF; 
    xMinusUKFAdd = xPlusUKFAdd; 
    PMinusUKFAdd = PPlusUKFAdd;
%     xMinusUKF_sysID = xPlusUKF_sysID; % is done automatically by ukf
%     object
%     PMinusUKF_sysID = PPlusUKF_sysID;
    xMinusSRUKFAdd = xPlusSRUKFAdd; 
    SMinusSRUKFAdd = SPlusSRUKFAdd;
%     xMinuscUKFAdd = xPluscUKFAdd; 
%     PMinuscUKFAdd = PPluscUKFAdd;
    xMinusUKFAug = xPlusUKFAug; 
    PMinusUKFAug = PPlusUKFAug;
%     xMinuscUKFAug = xPluscUKFAug; 
%     PMinuscUKFAug = PPluscUKFAug;
    xMinusUKFFullyAug = xPlusUKFFullyAug; 
    PMinusUKFFullyAug = PPlusUKFFullyAug;
%     xMinuscUKFFullyAug = xPluscUKFFullyAug; 
%     PMinuscUKFFullyAug = PPluscUKFFullyAug;
%     xMinusCKF = xPlusCKF; % is done automatically by ckf object
%     PMinusCKF = PPlusCKF;
end

%% 3) Plots der Ergebnisse

% figure
% % Biomasse
% subplot(311)
% plot(t,MESS(1,:)'.*STATES(3,:)','ok')
% hold on
% plot(t,STATES(1,:)','k')
% plot(t,ESTIMATESEKF(1,:)','r')
% plot(t,ESTIMATESUKFAdd(1,:)','b')
% plot(t,ESTIMATESEKF(1,:)'+sqrt(COVARIANCEEKF(1,:))',':r')
% plot(t,ESTIMATESUKFAdd(1,:)'+sqrt(COVARIANCEUKFAdd(1,:))',':b')
% plot(t,ESTIMATESEKF(1,:)'-sqrt(COVARIANCEEKF(1,:))','--r')
% plot(t,ESTIMATESUKFAdd(1,:)'-sqrt(COVARIANCEUKFAdd(1,:))','--b')
% hold off
% ylabel('Biomasse m_X [g]')
% xlim([0,t(end)])
% title('Simulierte und geschätzte Zustände')
% legend('Messung', 'Simulation', 'EKF', 'UKF', 'EKF +/- 1 \sigma', 'UKF +/- 1 \sigma')
% 
% % Substrat
% subplot(312)
% plot(t,STATES(2,:)','k')
% hold on
% plot(t,ESTIMATESEKF(2,:)','r')
% plot(t,ESTIMATESUKFAdd(2,:)','b')
% plot(t,ESTIMATESEKF(2,:)'+sqrt(COVARIANCEEKF(2,:))',':r')
% plot(t,ESTIMATESUKFAdd(2,:)'+sqrt(COVARIANCEUKFAdd(2,:))',':b')
% plot(t,ESTIMATESEKF(2,:)'-sqrt(COVARIANCEEKF(2,:))','--r')
% plot(t,ESTIMATESUKFAdd(2,:)'-sqrt(COVARIANCEUKFAdd(2,:))','--b')
% hold off
% ylabel('Substrat m_S [g]')
% xlim([0,t(end)])
% 
% % Volumen
% subplot(313)
% plot(t,STATES(3,:)','k')
% hold on
% plot(t,MESS(2,:)','ok')
% plot(t,ESTIMATESEKF(3,:)','r')
% plot(t,ESTIMATESUKFAdd(3,:)','b')
% plot(t,ESTIMATESEKF(3,:)'+sqrt(COVARIANCEEKF(3,:))',':r')
% plot(t,ESTIMATESUKFAdd(3,:)'+sqrt(COVARIANCEUKFAdd(3,:))',':b')
% plot(t,ESTIMATESEKF(3,:)'-sqrt(COVARIANCEEKF(3,:))','--r')
% plot(t,ESTIMATESUKFAdd(3,:)'-sqrt(COVARIANCEUKFAdd(3,:))','--b')
% hold off
% % ylim([0.9,1.1])
% ylabel('Volumen [l]')
% xlim([0,t(end)])
% xlabel('Zeit [h]')

%% compare EKF and all 6 different UKF implementations: 
viridisColorPaletteHex = ["#f0f921", "#fdb42f", "#ed7953", "#cc4778", "#9c179e", "#5c01a6", "#0d0887"]; 
% created with 7 colors and plasma color palette at 
% https://waldyrious.net/viridis-palette-generator/

figure()
% Biomasse
subplot(311)
plot(t,MESS(1,:)'.*STATES(3,:)','ok', 'DisplayName','Messung')
hold on
plot(t,STATES(1,:)','k', 'DisplayName','Simulation')
plot(t,ESTIMATESEKF(1,:)', 'Color','red', 'DisplayName','EKF')
plot(t,ESTIMATESUKFAdd(1,:)', 'Color',viridisColorPaletteHex(1), ...
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','UKF-add')
plot(t,ESTIMATESUKF_sysID(1,:)', 'Color',viridisColorPaletteHex(6), ...
    'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','UKF-sysID')
plot(t,ESTIMATESSRUKFAdd(1,:)', 'Color','green', ... % viridisColorPaletteHex(1)
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','SR-UKF-add')
% plot(t,ESTIMATEScUKFAdd(1,:)', 'Color',viridisColorPaletteHex(2), ...
%     'LineStyle', ':', 'LineWidth',1.6, 'DisplayName','cUKF-add')
plot(t,ESTIMATESUKFAug(1,:)', 'Color',viridisColorPaletteHex(3), ...
    'LineStyle', '-.', 'LineWidth',0.8, 'DisplayName','UKF-aug')
% plot(t,ESTIMATEScUKFAug(1,:)', 'Color',viridisColorPaletteHex(4), ...
%     'LineStyle', '-.', 'LineWidth',1.6, 'DisplayName','cUKF-aug')
plot(t,ESTIMATESUKFFullyAug(1,:)', 'Color',viridisColorPaletteHex(5), ...
    'LineStyle', '--', 'LineWidth',0.8, 'DisplayName','UKF-fully-aug')
% plot(t,ESTIMATEScUKFFullyAug(1,:)', 'Color',viridisColorPaletteHex(5), ...
%     'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','cUKF-fully-aug')
% plot(t,ESTIMATESCKF(1,:)', 'Color',viridisColorPaletteHex(6), ...
%     'LineStyle', '--', 'LineWidth',1.6)
ylabel('Biomasse m_X [g]')
xlim([0,t(end)])
title('Simulierte und geschätzte Zustände')
legend()

% Substrat
subplot(312)
plot(t,STATES(2,:)','k', 'DisplayName','Simulation')
hold on
plot(t,ESTIMATESEKF(2,:)', 'Color','red', 'DisplayName','EKF')
plot(t,ESTIMATESUKFAdd(2,:)', 'Color', viridisColorPaletteHex(1), ...
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','UKF-add')
plot(t,ESTIMATESUKF_sysID(2,:)', 'Color',viridisColorPaletteHex(6), ...
    'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','UKF-sysID')
plot(t,ESTIMATESSRUKFAdd(2,:)', 'Color','green', ... % viridisColorPaletteHex(1)
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','SR-UKF-add')
% plot(t,ESTIMATEScUKFAdd(2,:)', 'Color',viridisColorPaletteHex(2), ...
%     'LineStyle', ':', 'LineWidth',1.6, 'DisplayName','cUKF-add')
plot(t,ESTIMATESUKFAug(2,:)', 'Color',viridisColorPaletteHex(3), ...
    'LineStyle', '-.', 'LineWidth',0.8, 'DisplayName','UKF-aug')
% plot(t,ESTIMATEScUKFAug(2,:)', 'Color',viridisColorPaletteHex(4), ...
%     'LineStyle', '-.', 'LineWidth',1.6, 'DisplayName','cUKF-aug')
plot(t,ESTIMATESUKFFullyAug(2,:)', 'Color',viridisColorPaletteHex(5), ...
    'LineStyle', '--', 'LineWidth',0.8, 'DisplayName','UKF-fully-aug')
% plot(t,ESTIMATEScUKFFullyAug(2,:)', 'Color',viridisColorPaletteHex(5), ...
%     'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','cUKF-fully-aug')
% plot(t,ESTIMATESCKF(2,:)', 'Color',viridisColorPaletteHex(6), ...
%     'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','CKF')
ylabel('Substrat m_S [g]')
title('Simulierte und geschätzte Zustände')
legend()
xlim([0,t(end)])

% Volumen
subplot(313)
plot(t,STATES(3,:)','k', 'DisplayName','Simulation')
hold on
plot(t,MESS(2,:)','ok', 'DisplayName','Messung')
plot(t,ESTIMATESEKF(3,:)', 'Color','red', 'DisplayName','EKF')
plot(t,ESTIMATESUKFAdd(3,:)', 'Color', viridisColorPaletteHex(1), ...
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','UKF-add')
plot(t,ESTIMATESUKF_sysID(3,:)', 'Color',viridisColorPaletteHex(6), ...
    'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','UKF-sysID')
plot(t,ESTIMATESSRUKFAdd(3,:)', 'Color','green', ... % viridisColorPaletteHex(1)
    'LineStyle', ':', 'LineWidth',0.8, 'DisplayName','SR-UKF-add')
% plot(t,ESTIMATEScUKFAdd(3,:)', 'Color',viridisColorPaletteHex(2), ...
%     'LineStyle', ':', 'LineWidth',1.6, 'DisplayName','cUKF-add')
plot(t,ESTIMATESUKFAug(3,:)', 'Color',viridisColorPaletteHex(3), ...
    'LineStyle', '-.', 'LineWidth',0.8, 'DisplayName','UKF-aug')
% plot(t,ESTIMATEScUKFAug(3,:)', 'Color',viridisColorPaletteHex(4), ...
%     'LineStyle', '-.', 'LineWidth',1.6, 'DisplayName','cUKF-aug')
plot(t,ESTIMATESUKFFullyAug(3,:)', 'Color',viridisColorPaletteHex(5), ...
    'LineStyle', '--', 'LineWidth',0.8, 'DisplayName','UKF-fully-aug')
% plot(t,ESTIMATEScUKFFullyAug(3,:)', 'Color',viridisColorPaletteHex(5), ...
%     'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','cUKF-fully-aug')
% plot(t,ESTIMATESCKF(3,:)', 'Color',viridisColorPaletteHex(6), ...
%     'LineStyle', '--', 'LineWidth',1.6, 'DisplayName','CKF')
ylabel('Volumen [l]')
legend()
xlim([0,t(end)])
xlabel('Zeit [h]')

%% Kalman Gains
% figure
% subplot(311)
% stairs(t,GAINEKF(1,:),'b')
% hold on
% stairs(t,GAINUKFAdd(1,:),'b','LineStyle','-.')
% title('Korrektur der Zustände: K*v')
% xlim([0,t(end)])
% ylabel('Biomasse m_X [g]')
% legend('EKF', 'UKF')
% grid on
% %
% subplot(312)
% stairs(t,GAINEKF(2,:),'r')
% hold on
% stairs(t,GAINUKFAdd(2,:),'r','LineStyle','-.')
% xlim([0,t(end)])
% ylabel('Substrat m_S [g]')
% grid on
% %
% subplot(313)
% stairs(t,GAINEKF(3,:),'g')
% hold on
% stairs(t,GAINUKFAdd(3,:),'g','LineStyle','-.')
% ylabel('Substrat m_S [g]')
% xlim([0,t(end)])
% xlabel('Zeit [h]')
% grid on
