%% Version
% (R2022b) Update 2
% Erstelldatum: 22.02.23
% Autor: Simon Hellmann

%% Bootstrap Analysis with Max. Likelihood Estimation

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
figure

% volume flow:
subplot(3,2,1)
plot(MESS.tOn,MESS.yMeasOn(:,1),'k.','DisplayName','Messung')
hold on
plot(MESS.tOn,MESS.yCleanOn(:,1),'r-','LineWidth',1.5,'DisplayName','clean online value')
ylabel('volume flow in l/h')
ylim([80,520])
legend('Location','NorthEast')

% p_ch4:
subplot(3,2,2)
yyaxis right    % wechsel über zur rechten y-Achse
stairs(MESS.u(:,1),MESS.u(:,2),'g','DisplayName','Fütterung')
ylabel('feed volume flow in l/h')
yyaxis left     % wechsel über zur linken y-Achse
plot(MESS.tOn,MESS.yMeasOn(:,2),'k.','DisplayName','online Messung')
hold on
plot(MESS.tOn,MESS.yCleanOn(:,2),'r-','DisplayName','ungestörter Wert')
ylabel('p_{ch4} in bar')

% p_co2:
subplot(3,2,3)
plot(MESS.tOn,MESS.yMeasOn(:,3),'k.','DisplayName','Messung')
hold on
plot(MESS.tOn,MESS.yCleanOn(:,3),'r-','LineWidth',1.5,'DisplayName','ungestörter Wert')
ylabel('p_{co2} in bar')
ylim([0.3,0.55])

% S_IN: 
subplot(3,2,4)
plot(MESS.tOff,MESS.yMeasOff(:,1),'kx','DisplayName','offline Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,1),'b-','LineWidth',1.5,'DisplayName','clean offline value')
ylabel('Inorganic Nitrogen in g/L')
ylim([0.2,1.1])
legend('Location','SouthEast')

% TS: 
subplot(3,2,5)
plot(MESS.tOff,MESS.yMeasOff(:,2),'kx','DisplayName','Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,2),'b-','LineWidth',1.5,'DisplayName','ungestörter Wert')
ylabel('TS')
xlabel('time in d')

% oTS: 
subplot(3,2,6)
plot(MESS.tOff,MESS.yMeasOff(:,3),'kx','DisplayName','Messung')
hold on
plot(MESS.tOff,MESS.yCleanOff(:,3),'b-','LineWidth',1.5,'DisplayName','ungestörter Wert')
ylabel('oTS')
xlabel('time in d')

sgtitle('Messwerte und Simulation mit p0 aus MESS')

%% Initiale Parameteridentifikation mit Maximum Likelihood Estimation 
% Markov-Schätzung (bekannte Messunsicherheit C)

C = MESS.C;     % cov. matrix of measurement noise

% Startwerte der Optimierung und Beschränkungen (Diss Sören, Tab. 3.6a)
% Reihenfolge: kch, kpr, kli, kdec
pInit = [0.1, 1,       0.5,    0.2];% initial guess for parameter vector
pLB = [0.02, 0.01, 0.01, 0.001];    % untere Grenze Parametervektor
pUB = [2.88, 2.688, 0.76, 1];       % obere Grenze Parametervektor

% Optimierung 
costFun = @(p) costFunMarkovOnOffline(p,MESS,pFix,C);

% use incorrect (off) initial values of parameters and check, what MLE converges to:
pOpt0 = fmincon(costFun, pInit,[],[],[],[], pLB, pUB);

%% Validate initial parameter estimation
% Simulation with newly identified parameters:
[SIMOPT1a] = biogasmodell_lsg_zdglOnOffline(MESS,pOpt0,pFix);

%% plot all 6 measurements on the basis of the current set of parameters:
figure

% volume flow:
subplot(3,2,1)
plot(MESS.tOn,MESS.yMeasOn(:,1),'b.','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOn,SIMOPT1a.yOn(:,1),'r-','LineWidth',1.5,'DisplayName','Simulation')
ylabel('volume flow in l/h')
ylim([80,520])

% p_ch4:
subplot(3,2,2)
yyaxis right    % wechsel über zur rechten y-Achse
stairs(MESS.u(:,1),MESS.u(:,2),'g','DisplayName','Fütterung')
ylabel('feed volume flow in l/h')
yyaxis left     % wechsel über zur linken y-Achse
plot(MESS.tOn,MESS.yMeasOn(:,2),'b.','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOn,SIMOPT1a.yOn(:,2),'r-','DisplayName','Simulation')
ylabel('p_{ch4} in bar')
legend('Location','NorthEast')

% p_co2:
subplot(3,2,3)
plot(MESS.tOn,MESS.yMeasOn(:,3),'k.','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOn,SIMOPT1a.yOn(:,3),'b-','LineWidth',1.5,'DisplayName','Simulation')
ylabel('p_{co2} in bar')
ylim([0.3,0.55])
% legend('Location','SouthEast')

% S_IN: 
subplot(3,2,4)
plot(MESS.tOff,MESS.yMeasOff(:,1),'kx','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOff,SIMOPT1a.yOff(:,1),'b-','LineWidth',1.5,'DisplayName','Simulation')
ylabel('Inorganic Nitrogen in g/L')
ylim([0.2,1.1])

% TS: 
subplot(3,2,5)
plot(MESS.tOff,MESS.yMeasOff(:,2),'kx','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOff,SIMOPT1a.yOff(:,2),'r-','LineWidth',1.5,'DisplayName','Simulation')
ylabel('TS')

% oTS: 
subplot(3,2,6)
plot(MESS.tOff,MESS.yMeasOff(:,3),'kx','DisplayName','Messung')
hold on
plot(SIMOPT1a.tOff,SIMOPT1a.yOff(:,3),'c-','LineWidth',1.5,'DisplayName','Simulation')
ylabel('oTS')

sgtitle('Messwerte und Simulation mit optimierten Parametern')

disp('Parameter aus der Identifikation')
disp(pOpt0)

%% actual bootstrap: 
N = 1000;           % number of synthetic experiments
m = numel(pInit);   % number of parameters
pMat = nan(N,m);    % allocate memory for estimated parameters

rng('default');     % fix seed for random number generation (for replicable results)

% create new noisy measurements once:
[yCleanOn, yCleanOff] = computeSyntheticCleanMeasurements(MESS,pOpt0,pFix); 

% add to synthetic measurements random noise (N times) and perform N times
% MLE estimations
tic
for k = 1:N
    
    % add random, normally distributed noise to measurements:  
    [yMeasOn, yMeasOff] = addNoiseToCleanMeasurements(yCleanOn,yCleanOff,C);

    % overwrite measurements in struct: 
    MESS.yMeasOn = yMeasOn; 
    MESS.yMeasOff = yMeasOff;

    newCostFun = @(p) costFunMarkovOnOffline(p,MESS,pFix,C);

    % perform new max. likelihood estimation with initial optimal parameter 
    % sets as starting point:
    pMat(k,:) = fmincon(newCostFun, pOpt0,[],[],[],[], pLB, pUB);

end
toc

% compute mean and covariance of parameter estimate: 
pMean = mean(pMat); 
diffpFromMean = pMat - pMean;   % standing matrix (N,m)
pCov = 1/(N-m) * (diffpFromMean')*diffpFromMean; 

% compute ideal normal distributions of parameters based on pMat and pCov
pSigma = sqrt(diag(pCov))';
pEval = nan(N,m); % evaluation points of normal distribution (pdf)
pNormal = nan(N,m); 
mypNormal = nan(N,m); 

% compute 1D normal distributions of individual parameters:
for k = 1:m
    pMin = pMean(k) - 2*pSigma(k); 
    pMax = pMean(k) + 2*pSigma(k); 
    pEval(:,k) = linspace(pMin, pMax, N); 
    pNormal(:,k) = normpdf(pEval(:,k),pMean(k),pSigma(k)); 
    mypNormal(:,k) = normalPDF(pEval(:,k),pMean(k),pSigma(k));
end

%% plot rel frequencies (pdf) of parameters
figure

% kch
subplot(2,2,1)
nBins = 20; 
h1 = histogram(pMat(:,1),nBins);        
% h1.Normalization = 'probability'; % realized distribution (rel. frequency)  
% hold on 
% plot(pEval(:,1), pNormal(:,1))
% legend('hist', 'Gauß')
xlabel('k_{ch}')

% kpr
subplot(2,2,2)
nBins = 20; 
histogram(pMat(:,2),nBins) % realized distribution
xlabel('k_{pr}')

% kli
subplot(2,2,3)
nBins = 20; 
histogram(pMat(:,3),nBins) % realized distribution
xlabel('k_{li}')

% kdec
subplot(2,2,4)
nBins = 20; 
histogram(pMat(:,4),nBins) % realized distribution
xlabel('k_{dec}')

%% re-simulate the system with the mean parameter values: 
[SIMOPTFinal] = biogasmodell_lsg_zdglOnOffline(MESS,pMean,pFix);

%% plot all 6 measurements on the basis of the current set of parameters:
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
figure

% volume flow:
subplot(3,2,1)
scatter(MESS.tOn,MESS.yMeasOn(:,1)/24,'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(SIMOPT1a.tOn,SIMOPTFinal.yOn(:,1)/24,'DisplayName','Simulation',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
ylabel('Gasvolumenstrom [l/h]')
yyaxis right
stairs(MESS.u(:,1),MESS.u(:,2)/24,'g','DisplayName','Fütterung',...
    'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5);
set(gca,'YColor','k')
ylabel({'Fütterungs-','Volumenstrom [l/h]'})

% p_ch4:
subplot(3,2,2)
yyaxis right    % wechsel über zur rechten y-Achse
stairs(MESS.u(:,1),MESS.u(:,2)/24,'g','DisplayName','Fütterung',...
    'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5);
ylabel({'Fütterungs-','Volumenstrom [l/h]'})
set(gca,'YColor','k')
yyaxis left     % wechsel über zur linken y-Achse
scatter(MESS.tOn,MESS.yMeasOn(:,2),'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
hold on
plot(SIMOPTFinal.tOn,SIMOPTFinal.yOn(:,2),'DisplayName','Simulation', ...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8); 
ylabel('p_{ch4} [bar]')
legend('Location','NorthEast')

% p_co2:
subplot(3,2,3)
scatter(MESS.tOn,MESS.yMeasOn(:,3),'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on
plot(SIMOPTFinal.tOn,SIMOPTFinal.yOn(:,3),'DisplayName','Simulation',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8);
ylabel('p_{co2} [bar]')
ylim([0.3,0.55])
% legend('Location','SouthEast')

% S_IN: 
subplot(3,2,4)
scatter(MESS.tOff,MESS.yMeasOff(:,1),'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on
plot(SIMOPTFinal.tOff,SIMOPTFinal.yOff(:,1),'DisplayName','Simulation',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('Anorg. Stickstoff [g/l]')
ylim([0.2,1.1])

% TS: 
subplot(3,2,5)
scatter(MESS.tOff,MESS.yMeasOff(:,2),'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on
plot(SIMOPTFinal.tOff,SIMOPTFinal.yOff(:,2),'DisplayName','Simulation',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('TS [-]')

% oTS: 
subplot(3,2,6)
scatter(MESS.tOff,MESS.yMeasOff(:,3), 'DisplayName','Messung',...
    'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
hold on
plot(SIMOPTFinal.tOff,SIMOPTFinal.yOff(:,3),'DisplayName','Simulation',...
    'LineStyle','-', 'Color', colorPaletteHex(3), 'LineWidth',0.8)
ylabel('oTS [-]')

sgtitle('Messwerte und Simulation mit Mittelwerten der Parametern nach Bootstrap')

disp('Mean Parameters after bootstrap:')
disp(pMean)


%% save results
save('bootstrapResults', 'pCov', 'pMean')

