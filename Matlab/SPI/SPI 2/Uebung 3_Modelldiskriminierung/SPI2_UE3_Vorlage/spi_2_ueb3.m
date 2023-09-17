%% spi_2_ueb3_LSG - Modelldiskriminierung mittels Akaike Kriterium
% Modelldiskriminierung für die Modellkandidaten 1, 2, 3

close all
clear all
clc

warning('off','optim:fmincon:SwitchingToMediumScale');
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

%% 0) Pulsexperiment erzeugen aus Simulationsdateien

% Laden der Messdaten des Pulsexperiments
% Struct MESS besteht aus den Feldern
% -- t:	Messzeitpunkte
% -- x0: Anfangsbedingungen (Reihenfolge der Zeilen: m_X, m_S, V)
% -- x: Werte der Zustände zu den Messzeitpunkten
% -- x_sim: Zustand über kompletten Horizont (Reihenfolge: t_sim, m_X, m_S, V)
% -- u: Stellgröße (Reihenfolge: t_u, u, c_Sin)
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

%drawnow;

%% Parameteridentifikation und AIC Berechnung

% Berechnung Akaike in einer Schleife
for iM = 1:3
    Modell = iM;
    
    if iM == 1       
        pLB = [1e-3;1e-3;1e-2];
%         pUB = [0.5;0.5;1]; 
        pUB = [0.2;0.05;0.4];
    elseif iM == 2
        pLB = [1e-3;1e-3;1e-1;1e-3];
%         pUB = [0.5;0.5;1;10];
        pUB = [0.2;0.1;0.4;0.01];
    elseif iM == 3   
        pLB = [1e-3;1e-3;1e-1;1e-3;1e-3];
%         pUB = [0.5;0.5;1;10;10];         
        pUB = [0.2;0.1;0.4;0.01;5]; 
    end
    
    %% Aufgabe: 1c
    p0 = 
    pOpt = fmincon(@guete_pi_MLE,p0,[],[],[],[],pLB,pUB,[],[],MESS,Modell);
       
    % AIC durch Aufruf der Funktion
    [I,C,AIC_c(iM)] = aic_ueb3(pOpt,MESS,Modell);
    
    % Speicherung der optimierten Parameter in einer Cell-Struct
    p_Opt(iM).p = pOpt;
end

AIC_min = min(AIC_c)
% AIC_min = find(AIC == min(AIC))


%% Aufgabe 2a: Berechnung der Gewichtung

aic_distance(1) = 
w_aic(1) =

%% Anwendung auf MMTP

% Vorbereitung der Prozessoptimierung (Trajektorienplanung: TP)
% Bei einer TP mit unterschiedlichen Modellen (Multi-Modell) bezeichnet man
% die TP als MMTP (Multimodell-Trajektorienplanung)
% Zusammenbau eines Structs als Übergabeelement der TP
% Felder:
% -- t: Messzeitpunkte
TP.t = 0:5:50;
% -- x0: Anfangsbedingungen
TP.x0 = MESS.x0;
TP.x0(2) = TP.x0(2) / 5;
% -- u: Stellgrößen, Aufbau wie in MESS
tu = 0:1:50;
u0 = 0.02 * ones(size(tu));
TP.u = [tu; u0; 100 * ones(size(tu))];

% -- CONS: Struct für Beschränkungen mit den Feldern
% ---- umin/umax: der Stellgrößen
TP.CONS.umin = 0 * ones(size(tu));
TP.CONS.umax = 0.04 * ones(size(tu));

% ---- xmin/xmax: der Zustände
TP.CONS.xmin = [0;0;2];
TP.CONS.xmax = [35;25;4];

% ---- ymin/ymax: der Messgrößen
TP.CONS.ymin = [-Inf;-Inf];
TP.CONS.ymax = [Inf;Inf];


% -- p: Parameterwerte
TP.p = p_Opt;

optTP = MMTP(TP,w_aic)

% Plotten der Optimierungsergebnisse (Trajektorienplanung)
figure
stairs(optTP,'k')
set(gca,'XLim',[0 50],'YLim',[0 0.05])
xlabel('t in h')
ylabel('u in l/h')
