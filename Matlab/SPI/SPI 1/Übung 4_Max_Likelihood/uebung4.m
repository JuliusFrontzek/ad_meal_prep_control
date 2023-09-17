%% SPI Übung 4 - Parameteridentifikation

close all
clear all
clc

warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

%% 1) Pulsexperiment

% Laden der Messdaten des Pulsexperiments
% Struct MESS besteht aus den Feldern
% -- t:	Messzeitpunkte
% -- x0: Anfangsbedingungen (Reihenfolge der Zeilen: m_X, m_S, V)
% -- x: Werte der Zustände zu den Messzeitpunkten
% -- x_sim: Zustand über kompletten Horizont (Reihenfolge: t_sim, m_X, m_S, V)
% -- u: Stellgröße (Reihenfolge: t_u, u, c_Szu)
% -- y: Messungen zu den Messzeitpunkten

load Messung_Biomodell

figure
subplot(3,1,1)
plot(MESS.t,MESS.y(1,:),'ko',MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k')
title('Messwerte und wahrer Verlauf')
ylabel('c_X in g/l')
legend('Messung','wahrer Wert','Location','SouthEast')

subplot(3,1,2)
plot(MESS.t,MESS.y(2,:),'ko',MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k')
ylabel('c_S in g/l')

subplot(3,1,3)
stairs(MESS.u(1,:),MESS.u(2,:),'k')
set(gca,'XLim',[0 50],'YLim',[0 0.3])
xlabel('t in h')
ylabel('u in l/h')

drawnow;

%% 1a:  Simulation mit den Startwerten der Optimierung
% Gütefunktional definiert in 'güte_pi_WLS.m'

% % Startwerte der Optimierung und Beschränkungen
% % Reihenfolge: mu_max, K_S, Y_XS
p0 = [0.2, 0.05, 0.4];      % Startwert Parametervektor
pLB = [1E-3, 1E-3, 1E-1];   % untere Grenze Parametervektor
pUB = [0.5, 0.5, 1];        % obere Grenze Parametervektor
% 
% Kovarianzmatrix des Messrauschens
C = 0.05 * eye(2);
invC = C\eye(length(C)); % berechnen der Inversen gleich einmal zu Beginn

%  Simulation
% -- Funktionsaufruf siehe biomodell_lsg_zdgl.m
[SIMp0] = biomodell_lsg_zdgl(MESS,p0);

figure
subplot(3,1,1)
plot(...
	MESS.t,MESS.y(1,:),'ko',...
	SIMp0.x_sim(1,:),SIMp0.x_sim(2,:)./SIMp0.x_sim(4,:),'r--'...
	)
set(gca,'XLim',[0 50])
ylabel('c_X in g/l')
legend('Messung','Simulation','Location','NorthWest')
title(['Simulation mit den Startwerten p0=[', num2str(p0(1)),';',num2str(p0(2)),';',num2str(p0(3)),']' ])

subplot(3,1,2)
plot(...
	MESS.t,MESS.y(2,:),'ko',...
	SIMp0.x_sim(1,:),SIMp0.x_sim(3,:)./SIMp0.x_sim(4,:),'r--'...
	)
set(gca,'XLim',[0 50])
ylabel('c_S in g/l')

subplot(3,1,3)
stairs(MESS.u(1,:),MESS.u(2,:),'r')
set(gca,'XLim',[0 50],'YLim',[0 0.3])
xlabel('t in h')
ylabel('u in l/h')

drawnow;

%% 1b:  Parameteridentifikation mit bekannter Messunsicherheit
% Optimierung
% -- Optimierer: fmincon -- Argumente in der Hilfe nachschauen
% -- zusätzliche Argumente der Gütefunktion (nicht die Parameter) werden
%    dahinter angehängt. (Inverse wird direkt übergeben)

% costFun = @(p) guete_pi_WLS(p,MESS,invC); 
% np = length(p0);    % number of parameters
% pOpt1 = fmincon(costFun, p0,[],[],[],[], pLB, pUB);
% % 
% % % Berechnung der Fisher'schen Informationsmatrix und Simulation des Systems:
% % % -- Funktionsaufruf siehe biomodell_FIM.m
% % % -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
% [SIMOPT1] = biomodell_lsg_zdgl(MESS,pOpt1);
% % 
% 
% % Plot
% figure
% subplot(3,1,1)
% plot(...
% 	MESS.t,MESS.y(1,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(1,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(2,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_X in g/l')
% legend('Messung','Simulation','Location','NorthWest')
% 
% subplot(3,1,2)
% plot(...
% 	MESS.t,MESS.y(2,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(2,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(3,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_S in g/l')
% 
% subplot(3,1,3)
% stairs(MESS.u(1,:),MESS.u(2,:),'k')
% set(gca,'XLim',[0 50],'YLim',[0 0.3])
% xlabel('t in h')
% ylabel('u in l/h')
% 
% set(gcf,'NextPlot','add');
% axes;
% h = title(['1a) Angepasst, bekannte Unsicherheit, p=[', num2str(pOpt1(1)),';',num2str(pOpt1(2)),';',num2str(pOpt1(3)),']']);
% set(gca,'Visible','off');
% set(h,'Visible','on'); 
% 
% drawnow;
% 
% disp('Parameter aus der Identifikation')
% disp(pOpt1)

%% 1c: Parameteridentifikation mit unbekannter Messunsicherheit

% Optimierung 
% costFun = @(p) computeCostFunctionMLE(p,MESS);     % die Funktion computeCostFunctionMLE 
% % gibt nur die erste von zwei output-Variablen von guete_pi_MLE wieder (skalarer Wert des Gütefunktionals,
% % damit daraus ein function handle werden kann.
% pOpt1a = fmincon(costFun, p0,[],[],[],[], pLB, pUB);
% 
% % erhalte die Kovarianzmatrix durch erneuten Aufruf der Gütefunktion: 
% [dummy,Ca] = guete_pi_MLE(pOpt1a,MESS);
% 
% % Simulation
% [SIMOPT1a] = biomodell_lsg_zdgl(MESS,pOpt1a);
% 
% figure
% subplot(3,1,1)
% plot(...
% 	MESS.t,MESS.y(1,:),'ko',...
% 	SIMOPT1a.t,SIMOPT1a.y(1,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1a.x_sim(1,:),SIMOPT1a.x_sim(2,:)./SIMOPT1a.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_X in g/l')
% legend('Messung','Simulation','Location','NorthWest')
% 
% subplot(3,1,2)
% plot(...
% 	MESS.t,MESS.y(2,:),'ko',...
% 	SIMOPT1a.t,SIMOPT1a.y(2,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1a.x_sim(1,:),SIMOPT1a.x_sim(3,:)./SIMOPT1a.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_S in g/l')
% 
% subplot(3,1,3)
% stairs(MESS.u(1,:),MESS.u(2,:),'k')
% set(gca,'XLim',[0 50],'YLim',[0 0.3])
% xlabel('t in h')
% ylabel('u in l/h')
% 
% set(gcf,'NextPlot','add');
% axes;
% h = title(['1c) Angepasst, unbekannte Unsicherheit, p=[', num2str(pOpt1a(1)),';',num2str(pOpt1a(2)),';',num2str(pOpt1a(3)),']' ]);
% set(gca,'Visible','off');
% set(h,'Visible','on'); 
% legend('pOpt1a','pOpt1a \pm \sigma_p','Location','SouthEast')
% drawnow;
% 
% disp('Parameter aus der Identifikation')
% disp(pOpt1a)

%% 2a) 2b) Mit einfachen Erhaltungsterm (Maintenance) erweitertes Modell
% in biomodell_zdgl ein kommentieren 

% Startwerte der Optimierung und Beschränkungen
% Reihenfolge: mu_max, K_S, Y_XS, muMmax (,K_M)
p0 = [0.2, 0.05, 0.4, 0.02]; 
pLB = [1E-3, 1E-3, 1E-1, 1E-3]; 
pUB = [0.5, 0.5, 1, 0.5];

% % Optimierung
% % -- Optimierer: fmincon -- Argumente in der Hilfe nachschauen
% % -- zusätzliche Argumente der Gütefunktion (nicht die Parameter) werden dahinter angehängt.
% costFun2b = @(p) guete_pi_WLS(p,MESS,invC);
% pOpt1 = fmincon(costFun2b,p0,[],[],[],[],pLB,pUB);
% % 
% % Berechnung der Fisher'schen Informationsmatrix und Simulation
% % -- Funktionsaufruf siehe biomodell_FIM.m
% % -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
% [SIMOPT1] = biomodell_lsg_zdgl(MESS,pOpt1);
% 
% figure
% subplot(3,1,1)
% plot(...
% 	MESS.t,MESS.y(1,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(1,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(2,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_X in g/l')
% legend('Messung','Simulation','Location','NorthWest')
% 
% subplot(3,1,2)
% plot(...
% 	MESS.t,MESS.y(2,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(2,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(3,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_S in g/l')
% 
% subplot(3,1,3)
% stairs(MESS.u(1,:),MESS.u(2,:),'k')
% set(gca,'XLim',[0 50],'YLim',[0 0.3])
% xlabel('t in h')
% ylabel('u in l/h')
% 
% set(gcf,'NextPlot','add');
% axes;
% h = title(['2b) Einfacher Maintenanceterm p=[', num2str(pOpt1(1)),';',num2str(pOpt1(2)),';',num2str(pOpt1(3)),';',num2str(pOpt1(4)),']' ]);
% set(gca,'Visible','off');
% set(h,'Visible','on'); 
% 
% drawnow;
% 
% disp('Parameter aus der Identifikation')
% disp(pOpt1)

%% 2c) Mit komplexen Erhaltungsterm erweitertes Modell
% in biomodell_zdgl einkommentieren

% Startwerte der Optimierung und Beschränkungen
% Reihenfolge: mu_max, K_S, Y_XS,mu_M_max,K_M
% p0 = [0.2, 0.05, 0.4, 0.02, 0.05]; 
% pLB = [1E-3, 1E-3, 1E-1, 1E-3, 1E-3];
% pUB = [0.5, 0.5, 1, 0.5, 0.5];
% 
% % Optimierung
% % -- Optimierer: fmincon -- Argumente in der Hilfe nachschauen
% % -- zusätzliche Argumente der Gütefunktion (nicht die Parameter) werden dahinter angehängt.
% costFun2c = @(p) guete_pi_WLS(p,MESS,invC);
% pOpt1 = fmincon(costFun2c,p0,[],[],[],[],pLB,pUB);
% 
% % Berechnung der Fisher'schen Informationsmatrix und Simulation
% % -- Funktionsaufruf siehe biomodell_FIM.m
% % -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
% [SIMOPT1] = biomodell_lsg_zdgl(MESS,pOpt1);
% 
% figure
% subplot(3,1,1)
% plot(...
% 	MESS.t,MESS.y(1,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(1,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(2,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_X in g/l')
% legend('Messung','Simulation','Location','NorthWest')
% 
% subplot(3,1,2)
% plot(...
% 	MESS.t,MESS.y(2,:),'ko',...
% 	SIMOPT1.t,SIMOPT1.y(2,:),'ro',...
% 	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
% 	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(3,:)./SIMOPT1.x_sim(4,:),'r'...
% 	)
% set(gca,'XLim',[0 50])
% ylabel('c_S in g/l')
% 
% subplot(3,1,3)
% stairs(MESS.u(1,:),MESS.u(2,:),'k')
% set(gca,'XLim',[0 50],'YLim',[0 0.3])
% xlabel('t in h')
% ylabel('u in l/h')
% 
% set(gcf,'NextPlot','add');
% axes;
% h = title(['2c) komplexer Maintenanceterm p=[', num2str(pOpt1(1)),';',num2str(pOpt1(2)),';',num2str(pOpt1(3)),';',num2str(pOpt1(4)),';',num2str(pOpt1(5)),']' ]);
% set(gca,'Visible','off');
% set(h,'Visible','on'); 
% 
% drawnow;
% disp('Parameter aus der Identifikation')
% disp(pOpt1)

%% 2d) (Erweiterung von c): Modell mit Monod-Maintenenca und Schätzung von C
% in biomodell_zdgl einkommentieren

% Startwerte der Optimierung und Beschränkungen
% Reihenfolge: mu_max, K_S, Y_XS,mu_M_max,K_M
p0 = [0.2, 0.05, 0.4, 0.02, 0.05]; 
pLB = [1E-3, 1E-3, 1E-1, 1E-3, 1E-3];
pUB = [0.5, 0.5, 1, 0.5, 0.5];

% Optimierung
% -- Optimierer: fmincon -- Argumente in der Hilfe nachschauen
% -- zusätzliche Argumente der Gütefunktion (nicht die Parameter) werden dahinter angehängt.
costFun2d = @(p) computeCostFunctionMLE(p,MESS);     % die Funktion computeCostFunctionMLE 
% gibt nur die erste von zwei output-Variablen von guete_pi_MLE wieder (skalarer Wert des Gütefunktionals,
% damit daraus ein function handle werden kann.
pOpt2d = fmincon(costFun2d,p0,[],[],[],[],pLB,pUB);

% Kovarianzmatrix durch erneuten Aufruf der Gütefunktion bestimmen:
[dummy,Cd] = guete_pi_MLE(pOpt2d,MESS);

% Berechnung der Fisher'schen Informationsmatrix und Simulation
% -- Funktionsaufruf siehe biomodell_FIM.m
% -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
[SIMOPT1] = biomodell_lsg_zdgl(MESS,pOpt2d);

figure
subplot(3,1,1)
plot(...
	MESS.t,MESS.y(1,:),'ko',...
	SIMOPT1.t,SIMOPT1.y(1,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(2,:)./SIMOPT1.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_X in g/l')
legend('Messung','Simulation','Location','NorthWest')

subplot(3,1,2)
plot(...
	MESS.t,MESS.y(2,:),'ko',...
	SIMOPT1.t,SIMOPT1.y(2,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(3,:)./SIMOPT1.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_S in g/l')

subplot(3,1,3)
stairs(MESS.u(1,:),MESS.u(2,:),'k')
set(gca,'XLim',[0 50],'YLim',[0 0.3])
xlabel('t in h')
ylabel('u in l/h')

set(gcf,'NextPlot','add');
axes;
h = title(['2d) komplexer Maintenanceterm, C schätzen, p=[', num2str(pOpt2d(1)),';',num2str(pOpt2d(2)),';',num2str(pOpt2d(3)),';',num2str(pOpt2d(4)),';',num2str(pOpt2d(5)),']' ]);
set(gca,'Visible','off');
set(h,'Visible','on'); 

drawnow;
disp('Parameter aus der Identifikation')
disp(pOpt2d)