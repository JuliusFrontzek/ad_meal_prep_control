%% SPI 1 Übung 5 - Parameteridentifikation und OVP am einfachen biologischen Beispiel

close all
clear all
clc

warning('on','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

load Messung_Biomodell
addpath('derivatives','functions')

% Laden der Messdaten des Pulsexperiments
% Struct MESS besteht aus den Feldern
% -- t:	Messzeitpunkte
% -- x0: Anfangsbedingungen (Reihenfolge der Zeilen: m_X, m_S, V)
% -- x: Werte der Zustände zu den Messzeitpunkten
% -- x_sim: Zustand über kompletten Horizont (Reihenfolge: t_sim, m_X, m_S, V)
% -- u: Stellgröße (Reihenfolge: t_u, u, c_Sin)
% -- y: Messungen zu den Messzeitpunkten

%% 1a,1b,1c
% siehe 
% biomodell_zdgl_XP.m 
% biomodell_FIM.m 
% Parameteranalyse.m 

%% 1d Optimiertes Modell aus Übung 5 bei bekannter Messunsicherheit

% % Kovarianzmatrix des Messrauschens
C = 0.05 * eye(2);
invC=C\eye(length(C));
% 
% % % Startwerte der Optimierung und Beschränkungen
% % % Reihenfolge: mu_max, K_S, Y_XS
% p0 = [0.2;0.1;0.4];
pLB = [1e-3;1e-3;1e-1];
pUB = [0.7;0.7;1];

% Ergebnis aus Optimierung war 
pOpt1 = [0.1182,0.4587,0.4570]';

% Berechnung der Fisher'schen Informationsmatrix und Simulation
% -- Funktionsaufruf siehe biomodell_FIM.m
% -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
[FM1,SIMOPT] = biomodell_FIM(MESS,pOpt1,invC);

figure
subplot(3,1,1)
plot(...
	MESS.t,MESS.y(1,:),'ko',...
	SIMOPT.t,SIMOPT.y(1,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
	SIMOPT.x_sim(1,:),SIMOPT.x_sim(2,:)./SIMOPT.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_X in g/l')
legend('Messung','Simulation optimierte Parameter','Location','NorthWest')
title(['1) Simulation mit optimierten Werten der Parameter p0=[', num2str(pOpt1(1)),';',num2str(pOpt1(2)),';',num2str(pOpt1(3)),']' ])

subplot(3,1,2)
plot(...
	MESS.t,MESS.y(2,:),'ko',...
	SIMOPT.t,SIMOPT.y(2,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
	SIMOPT.x_sim(1,:),SIMOPT.x_sim(3,:)./SIMOPT.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_S in g/l')

subplot(3,1,3)
stairs(MESS.u(1,:),MESS.u(2,:),'k')
set(gca,'XLim',[0 50],'YLim',[0 0.3])
xlabel('t in h')
ylabel('u in l/h')

drawnow;

disp('Parameter aus der Identifikation')
disp(['mu_{max} ',' K_S ',' Y_{XS}'])
disp(pOpt1)

%% 1e) Parameterunsicherheit

% Approximation von CV durch inverse Fishermatrix
CVopt1=inv(FM1);

% Berechnung wichtiger Kenngrößen über externe Funktion
[Corr1,StdDev1,relStddev1,EW1,EV1,CN1]=Parameteranalyse(CVopt1,pOpt1)

%´Grafische Darstellung der Parameterunsicherheit
figure
SD=1; %Anzahl Standardabweichungen 1= ca 65% 2= ca 95% 3= ca 99%
h = plot_gaussian_ellipsoid(pOpt1, CVopt1, SD);
xlabel('\mu_{max}');ylabel('K_S');zlabel('Y_{XS}')
title('1e) Parameterunsicherheit für Optimierung aus Übung 5')
grid on

%% 2a
% siehe
% OVP.m
% guete_ovp

%% 2b) Optimale Versuchsplanung

% Zusammenbau eines Structs als Übergabeelement der OVP
% Felder:
% -- t: Messzeitpunkte
VP.t = 0:5:50;
% -- x0: Anfangsbedingungen
VP.x0 = MESS.x0;
% -- u: Stellgrößen, Aufbau wie in MESS
tu = 0:5:50;
u0 = 0.02 * ones(size(tu));
VP.u = [tu; u0; 100 * ones(size(tu))];

% -- p: Parameterwerte
VP.p = pOpt1;

%% 2b: Definieren Sie die Beschränkungen
% -- CONS: Struct für Beschränkungen mit den Feldern
% ---- umin/umax: der Stellgrößen
VP.CONS.umin = 
VP.CONS.umax =

% ---- xmin/xmax: der Zustände
VP.CONS.xmin =
VP.CONS.xmax =

% ---- ymin/ymax: der Messgrößen
VP.CONS.ymin = 
VP.CONS.ymax = 

% -- mit Parametersatz aus der PI

%% 2b Durchführung der OVP
% -- Aufruf siehe OVP.m
% -- Ausgabe: Struct mit optimalem Profil und dazu passenden Mess- und
%	   Zustandsgrößen
SIMOPT1(2) = OVP(VP,FM1,invC);

%% 2b) Theoretische Ergebnisse nach der OVP 

% Ausgabe der berechneten Stellgrößen nach der OVP 
uVP=SIMOPT1(2).u(2,:);

% Ausrechnen der Fishermatrix nach der OVP
[FM2] = guete_ovp_FIM(VP,uVP,FM1,invC); 
disp('fishermatrix nach OVP')
disp(FM2)

% Approximation von CV durch inverse Fishermatrix
CVopt2=inv(FM2);

% Berechnung wichtiger Kenngrößen über externe Funktion
[Corr2,StdDev2,relStddev2,EW2,EV2,CN2]=Parameteranalyse(CVopt2,pOpt1)

%´Grafische Darstellung der Parameterunsicherheit
figure
SD=1; %Anzahl Standardabweichungen 1= ca 65% 2= ca 95% 3= ca 99%
h = plot_gaussian_ellipsoid(pOpt1, CVopt2, SD);
xlabel('\mu_{max}');ylabel('K_S');zlabel('Y_{XS}')
title('5) Erwartete Parameterunsicherheit nach 2. Experiment')
grid on

%% 2c) Durchführung des neuen Experiments und Vergleich mit dem alten Modell

% Erweitern des MESS-Structs
MESS(2).t = SIMOPT1(2).t;
MESS(2).x0 = SIMOPT1(2).x0;
MESS(2).u = SIMOPT1(2).u;

%% 2c Berechnung der Messwerte mit biomodell_messdaten.m
MESS(2) = 

% Plot des neu durchgeführten Experiments und Vergleich mit altem Modell
figure
subplot(3,1,1)
plot(...
	MESS(2).t,MESS(2).y(1,:),'ko',...
	SIMOPT1(2).t,SIMOPT1(2).y(1,:),'ro',...
	MESS(2).x_sim(1,:),MESS(2).x_sim(2,:)./MESS(2).x_sim(4,:),'k',...
	SIMOPT1(2).x_sim(1,:),SIMOPT1(2).x_sim(2,:)./SIMOPT1(2).x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_X in g/l')
legend('Messung','Simulation mit altem Modell','Location','NorthWest')
title('6) Messung und Simulation des neuen Versuchs mit altem Modell')

subplot(3,1,2)
plot(...
	MESS(2).t,MESS(2).y(2,:),'ko',...
	SIMOPT1(2).t,SIMOPT1(2).y(2,:),'ro',...
	MESS(2).x_sim(1,:),MESS(2).x_sim(3,:)./MESS(2).x_sim(4,:),'k',...
	SIMOPT1(2).x_sim(1,:),SIMOPT1(2).x_sim(3,:)./SIMOPT1(2).x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_S in g/l')

subplot(3,1,3)
stairs(MESS(2).u(1,:),MESS(2).u(2,:),'k')
set(gca,'XLim',[0 50])
xlabel('t in h')
ylabel('u in l/h')

drawnow;

%% 2d) Erneute Parameteridentifikation mit neuem Experiment

%% 2d: Aufruf von guete_pi_WLS
pOpt2 = fmincon(...)

% Analyse des neuen Modells

[FM3,SIMOPT2] = biomodell_FIM(MESS,pOpt2,invC);

% Approximation von CV durch inverse Fishermatrix
CVopt3=inv(FM3);

% Berechnung wichtiger Kenngrößen über externe Funktion
[Corr3,StdDev3,relStddev3,EW3,EV3,CN3]=Parameteranalyse(CVopt3,pOpt2)

%´Grafische Darstellung der Parameterunsicherheit
figure
SD=1; %Anzahl Standardabweichungen 1= ca 65% 2= ca 95% 3= ca 99%
h = plot_gaussian_ellipsoid(pOpt1, CVopt3, SD);
xlabel('\mu_{max}');ylabel('K_S');zlabel('Y_{XS}')
title('8) Parameterunsicherheit nach erneuter Identifikation')
grid on

%% Plot des neuen Modells

figure
for idxMess = 1:length(MESS)
% 	figure
	subplot(3,length(MESS),idxMess)
	plot(...
		MESS(idxMess).t,MESS(idxMess).y(1,:),'ko',...
		SIMOPT2(idxMess).t,SIMOPT2(idxMess).y(1,:),'ro',...
		MESS(idxMess).x_sim(1,:),MESS(idxMess).x_sim(2,:)./MESS(idxMess).x_sim(4,:),'k',...
		SIMOPT2(idxMess).x_sim(1,:),SIMOPT2(idxMess).x_sim(2,:)./SIMOPT2(idxMess).x_sim(4,:),'r'...
		)
	set(gca,'XLim',[0 50])
	ylabel('c_X in g/l')
	legend('Messung','Simulation','Location','NorthWest')
    title(['Experiment ', num2str(idxMess), ': Neues Modell'])
    
	subplot(3,length(MESS),idxMess+2)
	plot(...
		MESS(idxMess).t,MESS(idxMess).y(2,:),'ko',...
		SIMOPT2(idxMess).t,SIMOPT2(idxMess).y(2,:),'ro',...
		MESS(idxMess).x_sim(1,:),MESS(idxMess).x_sim(3,:)./MESS(idxMess).x_sim(4,:),'k',...
		SIMOPT2(idxMess).x_sim(1,:),SIMOPT2(idxMess).x_sim(3,:)./SIMOPT2(idxMess).x_sim(4,:),'r'...
		)
	set(gca,'XLim',[0 50])
	ylabel('c_S in g/l')
	subplot(3,length(MESS),idxMess+4)
	stairs(MESS(idxMess).u(1,:),MESS(idxMess).u(2,:),'k')
	set(gca,'XLim',[0 50])
	xlabel('t in h')
	ylabel('u in l/h')
	drawnow;
end % for idxMess


% for idxMess = 1:length(MESS)
% 	figure
% 	subplot(3,1,1)
% 	plot(...
% 		MESS(idxMess).t,MESS(idxMess).y(1,:),'ko',...
% 		SIMOPT2(idxMess).t,SIMOPT2(idxMess).y(1,:),'ro',...
% 		MESS(idxMess).x_sim(1,:),MESS(idxMess).x_sim(2,:)./MESS(idxMess).x_sim(4,:),'k',...
% 		SIMOPT2(idxMess).x_sim(1,:),SIMOPT2(idxMess).x_sim(2,:)./SIMOPT2(idxMess).x_sim(4,:),'r'...
% 		)
% 	set(gca,'XLim',[0 50])
% 	ylabel('c_X in g/l')
% 	legend('Messung','Simulation','Location','NorthWest')
%     title(['Experiment ', num2str(idxMess), ': Messung und Simulation mit neu identifiziertes Modell'])
%     
% 	subplot(3,1,2)
% 	plot(...
% 		MESS(idxMess).t,MESS(idxMess).y(2,:),'ko',...
% 		SIMOPT2(idxMess).t,SIMOPT2(idxMess).y(2,:),'ro',...
% 		MESS(idxMess).x_sim(1,:),MESS(idxMess).x_sim(3,:)./MESS(idxMess).x_sim(4,:),'k',...
% 		SIMOPT2(idxMess).x_sim(1,:),SIMOPT2(idxMess).x_sim(3,:)./SIMOPT2(idxMess).x_sim(4,:),'r'...
% 		)
% 	set(gca,'XLim',[0 50])
% 	ylabel('c_S in g/l')
% 	subplot(3,1,3)
% 	stairs(MESS(idxMess).u(1,:),MESS(idxMess).u(2,:),'k')
% 	set(gca,'XLim',[0 50])
% 	xlabel('t in h')
% 	ylabel('u in l/h')
% 	drawnow;
% end % for idxMess