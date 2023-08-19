%% Version
% (R2022b) Update 5
% Erstelldatum: 28.02.23
% Autor: Simon Hellmann

%% Parameteridentifikation mit Max. Likelihood Estimation
close all
clear all
clc

warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

%% Laden der künstlichen Messdaten
% Struct MESS besteht aus den Feldern
% -- tOn: online Messzeitpunkte
% -- tOff: offline Messzeitpunkte
% -- tSim: Zeitvektor über kompletten Horizont
% -- x0: Anfangsbedingungen 
% -- x: Werte der Zustände zu den online Messzeitpunkten
% -- xSim: Zustand über kompletten Horizont x_i_Sim
% -- u: Stellgröße (Reihenfolge: t_u, u, x_i_in)
% -- yCleanOn: online Ausgangsgrößen zu den Messzeitpunkten (ohne Rauschen)
% -- yCleanOff: offline Ausgangsgrößen zu den Messzeitpunkten (ohne Rauschen)
% -- yMeasOn: online Messungen zu den Messzeitpunkten (mit Rauschen)
% -- yMeasOff: offline Messungen zu den Messzeitpunkten (mit Rauschen)

load Messung_ADM1_R4

% plot all 6 measurements on the basis of the current set of parameters:
fh1 = figure();
tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact');
% biogas flow:
nexttile; %subplot(3,2,1)
plot(MESS.tOn,MESS.yMeasOn(:,1),'k.', 'DisplayName','Messung');
hold on
plot(MESS.tOn,MESS.yCleanOn(:,1),'r-', 'LineWidth',1.5, 'DisplayName','ungestörter Wert');
ylabel('biogas flow [L/h]');
ylim([80,520]);

% p_ch4:
nexttile; %subplot(3,2,2)
yyaxis right    % wechsel über zur rechten y-Achse
stairs(MESS.u(:,1),MESS.u(:,2),'g', 'DisplayName','Fütterung')
ylabel('feed volume flow [L/h]')
yyaxis left     % wechsel über zur linken y-Achse
plot(MESS.tOn,MESS.yMeasOn(:,2),'k.', 'DisplayName','online Messung')
hold on
plot(MESS.tOn,MESS.yCleanOn(:,2),'r-','DisplayName','ungestörter Wert')
ylabel('p_{ch4} [bar]')
legend('Location','NorthEast')

% p_co2:
nexttile; %subplot(3,2,3)
plot(MESS.tOn,MESS.yMeasOn(:,3),'k.', 'DisplayName','Messung');
hold on
plot(MESS.tOn,MESS.yCleanOn(:,3),'b-', 'LineWidth',1.5, 'DisplayName','ungestörter Wert');
ylabel('p_{co2} [bar]');
ylim([0.3,0.55]);

% S_IN: 
nexttile; %subplot(3,2,4)
plot(MESS.tOff,MESS.yMeasOff(:,1),'kx', 'DisplayName','offline Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,1),'b-', 'LineWidth',1.5, 'DisplayName','ungestörter Wert')
ylabel('Inorganic Nitrogen [g/L]')
ylim([0.2,1.1])
legend('Location','SouthEast')

% TS: 
nexttile; %subplot(3,2,5)
plot(MESS.tOff,MESS.yMeasOff(:,2),'kx', 'DisplayName','Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,2),'r-', 'LineWidth',1.5, 'DisplayName','ungestörter Wert')
ylabel('TS [kg/kg FM]');
xlabel('time [d]');

% oTS: 
nexttile; %subplot(3,2,6)
plot(MESS.tOff,MESS.yMeasOff(:,3),'kx', 'DisplayName','Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,3),'b-', 'LineWidth',1.5, 'DisplayName','ungestörter Wert')
ylabel('oTS [kg/kg TS]');
xlabel('time [d]');

sgtitle('Messwerte und Simulation mit p0 aus MESS')

%% Teste biogasmodell_lsg_zdgl mit den Standard-Werten der Paramter: 
% (kann man überspringen!)

% Startwerte der Optimierung und Beschränkungen
% Reihenfolge: kch, kpr, kli, kdec
p0  = [0.25, 0.2, 0.1, 0.02];       % Startwert Parametervektor
pLB = [0.02, 0.01, 0.01, 0.001];    % untere Grenze Parametervektor
pUB = [2.88, 2.688, 0.76, 1];       % obere Grenze Parametervektor

%% Parameteridentifikation mit MaxLikelihoodEstimation (unbekannte Messunsicherheit)

% Optimierung 
costFun = @(p) myGuete_pi_MLEOnOffline(p,MESS,pFix); 
% costFun = @(p) extractCostFunctionMLEOnOffline(p,MESS,pFix);     % die Funktion computeCostFunctionMLE 
% gibt nur die erste von zwei output-Variablen von guete_pi_MLE wieder (skalarer Wert des Gütefunktionals,
% damit daraus ein function handle werden kann.

% use incorrect (off) initial values of parameters and check, what MLE converges to:
pOff1 = [0.1, 1,       0.5,    0.2];
pOff2 = [1,   0.1,     0.01,   0.02];
pOff3 = [0.01,0.01,    0.1,    0.002];
% pOpt1 = fmincon(costFun, pOff1,[],[],[],[], pLB, pUB);
tic;
pOpt2 = fmincon(costFun, pOff2,[],[],[],[], pLB,pUB);
toc
% pOpt3 = fmincon(costFun, pOff3,[],[],[],[], pLB, pUB);

% erhalte die Kovarianzmatrix durch erneuten Aufruf der Gütefunktion: 
% [~,CEst] = myGuete_pi_MLEOnOffline(pOpt1,MESS,pFix);

% Simulation with newly identified parameters:
[SIMOPT1a] = biogasmodell_lsg_zdglOnOffline(MESS,pOpt2,pFix);

%% Parameteridentifikation mit MaxLikelihoodEstimation (bekannte Messunsicherheit)

% Optimierung 
C = MESS.C;     % covariance matrix of measurement noise
costFunMarkov = @(p) costFunMarkovOnOffline(p,MESS,pFix,C); 

% use incorrect (off) initial values of parameters and check, what MLE converges to:
pOffMarkov = [0.1, 1,       0.5,    0.2];
pOptMarkov = fmincon(costFunMarkov, pOffMarkov,[],[],[],[], pLB, pUB);

% erhalte die Kovarianzmatrix durch erneuten Aufruf der Gütefunktion: 
% [~,CEst] = myGuete_pi_MLEOnOffline(pOpt1,MESS,pFix);

% Simulation with newly identified parameters:
[SIMOPT1b] = biogasmodell_lsg_zdglOnOffline(MESS,pOptMarkov,pFix);

%% plot all 6 measurements with a Markov-based MLE:
fh2 = figure();
tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact');
% biogas flow:
nexttile; %subplot(3,2,1)
plot(MESS.tOn,MESS.yMeasOn(:,1),'k.', 'DisplayName','Messung')
hold on
plot(SIMOPT1b.tOn,SIMOPT1b.yOn(:,1),'c-', 'LineWidth',1.5, 'DisplayName','Simulation')
ylabel('biogas flow [L/h]');
ylim([80,520])

% p_ch4:
nexttile; %subplot(3,2,2)
yyaxis right    % wechsel über zur rechten y-Achse
stairs(MESS.u(:,1),MESS.u(:,2),'g', 'DisplayName','Fütterung')
ylabel('feed volume flow [L/h]')
yyaxis left     % wechsel über zur linken y-Achse
plot(MESS.tOn,MESS.yMeasOn(:,2),'k.', 'DisplayName','Messung')
hold on
plot(SIMOPT1b.tOn,SIMOPT1b.yOn(:,2),'r-', 'DisplayName','Simulation')
ylabel('p_{ch4} [bar]')
legend('Location','NorthEast')

% p_co2:
nexttile; %subplot(3,2,3)
plot(MESS.tOn,MESS.yMeasOn(:,3),'k.', 'DisplayName','Messung')
hold on
plot(SIMOPT1b.tOn,SIMOPT1b.yOn(:,3),'b-', 'LineWidth',1.5, 'DisplayName','Simulation')
ylabel('p_{co2} [bar]')
ylim([0.3,0.55])
% legend('Location','SouthEast')

% S_IN: 
nexttile; %subplot(3,2,4)
plot(MESS.tOff,MESS.yMeasOff(:,1),'kx', 'DisplayName','Messung')
hold on
plot(SIMOPT1b.tOff,SIMOPT1b.yOff(:,1),'b-', 'LineWidth',1.5, 'DisplayName','Simulation')
ylabel('Inorganic Nitrogen [g/L]')
ylim([0.2,1.1])

% TS: 
nexttile; %subplot(3,2,5)
plot(MESS.tOff,MESS.yMeasOff(:,2),'kx', 'DisplayName','Messung')
hold on
plot(SIMOPT1b.tOff,SIMOPT1b.yOff(:,2),'r-', 'LineWidth',1.5, 'DisplayName','Simulation')
ylabel('TS [kg/kg FM]');
xlabel('time [d]');

% oTS: 
nexttile; %subplot(3,2,6)
plot(MESS.tOff,MESS.yMeasOff(:,3),'kx','DisplayName','Messung')
hold on
plot(SIMOPT1b.tOff,SIMOPT1b.yOff(:,3),'c-', 'LineWidth',1.5, 'DisplayName','Simulation')
ylabel('oTS [kg/kg TS]');
xlabel('time [d]');

sgtitle('Messwerte und Simulation mit optimierten Parametern')

disp('Parameter-Schätzungen mit UNbekanntem C')
% disp(pOpt1)
disp(pOpt2)
% disp(pOpt3)

disp('Parameter-Schätzungen mit BEkanntem C')
disp(pOptMarkov)