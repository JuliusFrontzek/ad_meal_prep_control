function [ERGEBNIS] = biomodell_lsg_zdgl(MESS,p)
% integriere das Systemverhalten auf bei gegebenen Parametern und
% F�tterungsplan

% Anzahl der Zust�nde
n = length(MESS(1).x0);

np = length(p);   % Anzahl Parameter
% hier noch allgemeine routine schreiben zum finden von y
ny = 2;           % Anzahl Messgr��en

% Iteration �ber alle Experimente (sp�ter in der OVP wird in MESS mehr als
% ein einziges Experiment enthalten sein; dann erst wird diese Iteration
% relevant):
for idxMess = 1:length(MESS)
	
	% Anfangsbedingungen f�r DGLn mit Sensitivit�ten
	x0 = [MESS(idxMess).x0];	% Anfangszustand
    x_sim = [];                 % wird einmal der Eintrag eines Structs
	y = biomodell_mgl(x0);
    x = x0;
    
    % Iteration �ber alle Messzeitpunkte innerhalb eines Experiments: 
	for idxT = 2:length(MESS(idxMess).t)    % MESS(idxMess).t ist der Vektor der Messzeitpunkte
		
%% Simulation des Systems ohne Sensitivit�ten
		
		% Simulation zwischen zwei Messpunkten
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
        
        % Folgendes ist die alte Notation zur �bergabe von Parametern an
        % die odefun. Die aktuelle Notation ist unterhalb angegeben:
        [T,X] = ode45(@biomodell_zdgl,tspan,x0,[],MESS(idxMess).u,p); % eine 
        % Zeitspanne tspan hochintegrieren: liefert viele Punkte f�r X und T
        % das ist die alte Notation (wie auch in guete_pi_WLS bei ode15s.
        % Alternativ ginge auch diese (modernere) Notation:
%         odefun = @(t,x) biomodell_zdgl(t,x,MESS(idxMess).u,p); 
%         [T,X] = ode45(odefun, tspan, x0); 
		x0 = X(end,:)';     % neuer Anfangszustand = alter Endzustand
		x(:,idxT) = x0;     % Zustands-Trajektorie zu den Messzeitpunkten
		x_sim = [x_sim [T'; X(:,1:n)']];    % erweitere x_sim mit jeder Iteration um einige neue Punkte
        % erste Zeile: Zeitpunkte, weitere drei Zeilen: Trajektorie des
        % Zustandsvektors zu den Simulations-Zeitpunkten (x_sim enth�lt
        % viel mehr Punkte als x, welches nur die Zust�nde w�hrend der
        % Messzeitpunkte beschreibt)
		y(:,idxT) = biomodell_mgl(x0);
		
	end % for idxT
	
	% Ausgabe des Simulationsstructs
	ERGEBNIS(idxMess).t = MESS(idxMess).t;
	ERGEBNIS(idxMess).x0 = MESS(idxMess).x0;
	ERGEBNIS(idxMess).x = x;
	ERGEBNIS(idxMess).x_sim = x_sim;
	ERGEBNIS(idxMess).u = MESS(idxMess).u;
	ERGEBNIS(idxMess).y = y;
	
end % for idxMess

