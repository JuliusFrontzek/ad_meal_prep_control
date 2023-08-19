% create the MESS-vector for Kalman Filtering based on Sörens/Manuels model 

clc
clear
close all

%% define parameters: 
% kinetic constants [1/d] (Tab. B.7 in Sörens Diss): 
kch = 0.25; 
kpr = 0.2; 
kli = 0.1; 
kdec = 0.02; 
k = [kch, kpr, kli, kdec];

% Henry coefficients: [mol/l/bar] (Tab. B.7 in Sörens Diss)
Kch4 = 0.0011;      
Kco2 = 0.025; 

% miscellaneous parameters (Weinrich, 2017, Tab. B.7): 
R = 0.08315;    % id. gas constant [bar l/mol/K]
ph2o = 0;       % 0.0657 partial pressure of water in gas phase (saturated) [bar]
patm = 1.0133;  % atmospheric pressure [bar]
kla = 200;      % mass transfer coefficient [1/d]
kp = 5e4;       % friction parameter [l/bar/d]

T = 273.15;     % operating temperature [K]
Vl = 100;       % liquid volume, aus Sörens GitHub[l]
Vg = 10;        % gas volume, aus Sörens GitHub [l]
rho = 1000;     % density of digestate [g/l]

modParams    = [Kch4, Kco2, R, T, kla, kch, kdec, kli, kp, kpr, ph2o]; 
physParams   = [Vl, Vg, patm];
% combine all fixed parameters in a separate vector: 
pFix = [Kch4, Kco2, R, T, kla, kp, ph2o, Vl, Vg, patm, rho];

% combine Sörens params in struct: 
params.phyParams = physParams; 
params.modParams = modParams; 

% inlet concentrations [GitHub Sören]
%      S_ch4, S_IC, X_ch,   X_pr, X_li,     X_bac, S_ch4,g, S_co2,g
xIn = [0,     0,    23.398, 4.75, 1.381,    0,     0,       0]; % [g/l]

%% combine model parameters for compact control notation: 
c1 = 1/Vl;  
c2 = kla; 
c3 = kla*Kch4*R*T; 
c4 = kla*Kco2*R*T; 
c5 = kla*Vl/Vg; 
c6 = kp/patm*(R*T/16)^2;
c7 = 2*kp/patm*(R*T)^2/16/44;
c8 = kp/patm*(R*T/44)^2;
c9 = kp/patm*(R*T)/16*(2*ph2o - patm);
c10 = kp/patm*(R*T)/44*(2*ph2o - patm);
c11 = kp/patm*(ph2o - patm)*ph2o; 
c12 = R*T/16;
c13 = R*T/44;
c14 = rho; 
c15 = -kp/patm/Vg*(R*T/16)^2; 
c16 = -2*kp/patm/Vg*(R*T)^2/16/44; 
c17 = -kp/patm/Vg*(R*T/44)^2; 
c18 = -kp/patm/Vg*(R*T)/16*(2*ph2o - patm);
c19 = -kp/patm/Vg*(R*T)/44*(2*ph2o - patm);
c20 = -kla*Vl/Vg*Kch4*R*T - kp/patm/Vg*(ph2o - patm)*ph2o; 
c21 = -kla*Vl/Vg*Kco2*R*T - kp/patm/Vg*(ph2o - patm)*ph2o; 
c22 = Vl/Vg; 
c23 = -kp/patm/Vg*(ph2o - patm)*ph2o; 
c = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23];

a = [0.2482,  0.6809,  -0.0207,    -0.0456,    -1,      0,      0,       0.1372,    0,      0,      0; 
     0.3221,  0.7954,   0.1689,    -0.4588,     0,     -1,      0,       0.1723,    0,      0,      0; 
     0.6393,  0.5817,  -0.0344,    -0.4152,     0,      0,     -1,       0.2286,    0,      0,      0; 
     0,       0,        0,          0,          0.18,   0.77,   0.05,   -1,         0,      0,      0;
    -1,       0,        0,          0,          0,      0,      0,       0,         0,      c22,    0; 
     0,      -1,        0,          0,          0,      0,      0,       0,         0,      0,      c22]';

th = [kch, kpr, kli, kdec];     % kinetic parameters

% combine constant parameters in struct: 
AC.a = a; 
AC.c = c; 
AC.th = th;

% times: 
t0 = 0;     % starting time
tEnd = 7;   % [d] End of Simulation
tSS = 365;  % [d] Dauer, bis steady state auf jeden fall erreicht

dt = 15/60/24;    % sampling time [min], converted to [d]
tGrid = (t0:dt:tEnd)';  % time grid every 0.5h [d]. The model outputs will be evaluated here later on

totalFeed = 10;     % [L]

% construct vector of inlet concentrations 
nStates = length(xIn);
inputVectorSS = [totalFeed/tEnd,xIn]; 

% define function handle to determine the steady state:
odeFunSS = @(t,x) BMR4_AB_ess_ode(t, x, physParams, inputVectorSS, modParams);    

% determine steady state as initial value for simulation: 
tSpanSS = [t0,tSS]; 
x0Soeren = [0.091,0.508,3.26,0.956,0.413,2.569,0.315,0.78];  % Sörens GitHub
[tVecSS,xSS] = ode15s(odeFunSS,tSpanSS,x0Soeren); 

x0 = xSS(end,:);   % start in steady state 

%% Solve ODE at once:

% we dont feed anything anymore during the Abklinger:
inputVector = zeros(size(inputVectorSS)); 
modelOdeFun = @(t, x) BMR4_AB_ess_ode(t, x, physParams, inputVector, modParams);
[~, xSol] = ode15s(modelOdeFun, tGrid, x0);

%% get output variables

yClean = biogasmodell_mgl_ess_mat(xSol,pFix);

%% add noise to measurements acc to noise covariances: 

volFlowClean = yClean(:,1);
pCh4Clean = yClean(:,2); 
pCo2Clean = yClean(:,3);

% define properties of normal distributions of sensors (assume zero mean), 
% see Übersicht_Messrauschen.xlsx):
% std. deviations: 
sigmaV = 0.2; % Trommelgaszähler [L/h]
sigmaCh4 = 0.002;   % [%]
sigmaCo2 = 0.002;   % [%]

sigmas = [sigmaV, sigmaCh4, sigmaCo2]; 

nMeas = length(pCh4Clean);  % number of measurements
volFlowMeas = volFlowClean + sigmaV*randn(nMeas,1);
pCh4Meas = pCh4Clean + sigmaCh4*randn(nMeas,1); 
pCo2Meas = pCo2Clean + sigmaCo2*randn(nMeas,1); 

%% Plot results (partial pressures, gas volume flow, SIN and feedings)

% plot and compare the results from biogasmodell_mgl and BMR4_AB_mgl: 
figure()

% gas volume flow: 
subplot(3,1,1)
plot(tGrid,yClean(:,1),'r','LineWidth',2,'DisplayName','clean')
hold on
plot(tGrid, volFlowMeas, 'b-.','LineWidth',0.5,'DisplayName','noisy')
ylabel('gas vol flow in L/h')
xlabel('time [d]')
% legend('Location','NorthEast'); 

% pch4: 
subplot(3,1,2)
plot(tGrid,yClean(:,2),'r','LineWidth',2,'DisplayName','clean')
hold on
plot(tGrid, pCh4Meas, 'b-.','LineWidth',0.5,'DisplayName','noisy')
ylabel('p_{ch4} in bar')

% pco2:
subplot(3,1,3)
plot(tGrid,yClean(:,3),'r','LineWidth',2,'DisplayName','clean')
hold on
plot(tGrid, pCo2Meas, 'b-.','LineWidth',0.5,'DisplayName','noisy')
ylabel('p_{co2} in bar')
legend('Location','SouthEast'); 

sgtitle('Simulation eines Abklingers mit Sörens/Manuels Modell')


%% save results in struct MESS: 
MESS.t = tGrid; 
MESS.x0 = x0; 
MESS.x = xSol; 
% MESS.xSim = [tSim,xSim]; 
MESS.u = inputVector;    % u in [L/d]
MESS.sigma = sigmas; 
MESS.yClean = [volFlowClean, pCh4Clean, pCo2Clean]; 
MESS.yMeas = [volFlowMeas, pCh4Meas, pCo2Meas]; 

save('Messung_ADM1_R4_ess_Abklinger.mat', 'MESS', 'pFix', 'params', 'AC')
