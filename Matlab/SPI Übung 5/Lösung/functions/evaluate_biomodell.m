
MESS = struct('t',[],'x0',[],'x',[],'x_sim',[],'u',[],'y',[]);

MESS.t = [0:5:30 30.1 35:5:50];
MESS.x0 = [1;10;2];
MESS.u = [0 30 30.1;0 0.2 0;100 700 100];

MESS = biomodell_messdaten(MESS);

save Messung_Biomodell MESS

%% -- Bekannte Messunsicherheit

% Kovarianzmatrix des Messrauschens
C = 0.05 * eye(2);
invC=C\eye(length(C));

% Startwerte der Optimierung und Beschränkungen
% Reihenfolge: mu_max, K_S, Y_XS
p = [0.1;0.2;0.6];
pLB = [1e-3;1e-3;1e-1];
pUB = [0.5;0.5;1];

% Optimierung
% -- Optimierer: fmincon -- Argumente in der Hilfe nachschauen
% -- zusätzliche Argumente der Gütefunktion (nicht die Parameter) werden
%    dahinter angehängt.
pOpt1 = p;

% Berechnung der Fisher'schen Informationsmatrix und Simulation
% -- Funktionsaufruf siehe biomodell_FIM.m
% -- Ausgabe: FIM und Simulation als Struct, analog zu MESS
[FM1,SIMOPT1] = biomodell_FIM(MESS,pOpt1,invC);

% Relative Parameterunsicherheit
Gam1 = diag(sqrt(inv(diag(pOpt1)*FM1*diag(pOpt1))));

figure
subplot(3,2,1)
plot(...
	MESS.t,MESS.y(1,:),'ko',...
	SIMOPT1.t,SIMOPT1.y(1,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(2,:)./MESS.x_sim(4,:),'k',...
	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(2,:)./SIMOPT1.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_X in g/l')
legend('Messung','Simulation','Location','NorthWest')
title(['Simulation mit optimierten Werten der Parameter p0=[', num2str(pOpt1(1)),';',num2str(pOpt1(2)),';',num2str(pOpt1(3)),']' ])


subplot(3,2,3)
plot(...
	MESS.t,MESS.y(2,:),'ko',...
	SIMOPT1.t,SIMOPT1.y(2,:),'ro',...
	MESS.x_sim(1,:),MESS.x_sim(3,:)./MESS.x_sim(4,:),'k',...
	SIMOPT1.x_sim(1,:),SIMOPT1.x_sim(3,:)./SIMOPT1.x_sim(4,:),'r'...
	)
set(gca,'XLim',[0 50])
ylabel('c_S in g/l')

subplot(3,2,5)
stairs(MESS.u(1,:),MESS.u(2,:),'k')
set(gca,'XLim',[0 50],'YLim',[0 0.3])
xlabel('t in h')
ylabel('u in l/h')

subplot(3,2,[4 6])
y = MESS.y(2,:);
y(y < 0) = 0;
cSmax = max(y);
cS = linspace(0,cSmax,50);
plot(...
	cS,pOpt1(1)*cS./(pOpt1(2)+cS),'k',...
	cS,(1+Gam1(1))*pOpt1(1)*cS./((1+Gam1(2))*pOpt1(2)+cS),'k--',...
	cS,(1-Gam1(1))*pOpt1(1)*cS./((1-min([Gam1(2) 0.9]))*pOpt1(2)+cS),'k--',...
	y,pOpt1(1)*y./(pOpt1(2)+y),'ko')
xlabel('c_S in g/l')
ylabel('\mu(c_S) in 1/h')
legend('pOpt1','pOpt1 \pm \sigma_p','Location','SouthEast')

drawnow;