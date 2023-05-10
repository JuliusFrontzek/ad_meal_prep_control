% Simuliere ein ganzes Experiments (integriere dazu DGLs für x und X auf)
% und berechne die Fisher'sche Informationsmatrix FM 
% (einfaches biologisches Modell)

function [FM,ERGEBNIS] = biomodell_FIM(MESS,p,invC)
% FM - Fisher-Matrix
% ERGEBNIS - Struct mit allen Simulationsergebnissen (x,y,u,t,x0...)
% p - Parameter
% invC - Inverse der Kovarianzmatrix des Messrauschens
% MESS - struct mit folgenden Feldern:
% -- t:	Messzeitpunkte
% -- x0: Anfangsbedingungen (Reihenfolge der Zeilen: m_X, m_S, V)
% -- x: Werte der Zustände zu den Messzeitpunkten
% -- x_sim: Zustand über kompletten Horizont (Reihenfolge: t_sim, m_X, m_S, V)
% -- u: Stellgröße (Reihenfolge: t_u, u, c_Sin)
% -- y: Messungen zu den Messzeitpunkten

n = length(MESS(1).x0); % Anzahl der Zustände
m = length(p);          % Anzahl der Parameter

% Anfangsbedingungen der Sensitivitäten
XP0 = biomodell_dzdgl0dp(MESS(1).x0,p);     % liefert Nullmatrix mit Dim. (n x m)

% Initialisierung der FIM
FM = zeros(m,m);

% Iteration über alle Experimente (diese Iteration wird erst ab Aufgabe 2c mehrfach durchgeführt)
for idxMess = 1:length(MESS)  
	
	% Anfangsbedingungen für DGLn mit Sensitivitäten
	x0_ = [MESS(idxMess).x0;XP0(1:n*m)'];	% erweiterter Zustandsvektor x_... 
    % (aus x und X); die Spalten der Matrix werden untereinander gepackt, alle unter x
	
	x_ = x0_;               % erweiterter Zustandsvektor (x & X) zu den Messzeitpunkten
	x = MESS(idxMess).x0;   % reiner Zustandsvektor zu den Messzeitpunkten 
	x_sim = [];             % Platzhalter für vollständige Simulationsergebnisse 
    % mit vom Integrator gewählten Schrittweite: 
    % oberste Zeile: Simulations-Zeitpunkte T, alle weiteren 12 Zeilen: 
    % erst 3 Zeilen für Zustände zu Zeitpunkten T, dann 9 Zeilen für die
    % Zustands-Sensitivitäten Xp zu Zeitpunkten T
    
	y = biomodell_mgl(x0_); % eigentlich unnötig, den erweiterten Zustandsvektor zu ...
    % übergeben (nur x -also ohne Xp- würde reichen), denn die übrigen
    % Einträge werden überhaupt nicht verwendet
    
    % Iteration über alle Messzeitpunkte innerhalb eines Experiments:
	for idxT = 2:length(MESS(idxMess).t)    
%% Simulation des Systems mit Sensitivitäten
		
		% Simulation zwischen je zwei Messpunkten:
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
        % ode45 wird entspr. der alten Notation aufgerufen:
		[T,X] = ode45(@biomodell_zdgl_XP,tspan,x0_,[],MESS(idxMess).u,p);
		x0_ = X(end,:)';    % neuer (erweiterter) Startwert ist Endwert des alten Messintervalls:
		x_(:,idxT) = x0_;   % Auswertung zu Messzeitpunkten
		x(:,idxT) = x0_(1:n);   % nur die ersten drei Einträge sind wirklich 
        % der Zustand, alle weiteren sind die Sensitivität des Zustands Xp
		x_sim = [x_sim [T'; X(:,1:n)']]; % mit jeder Iteration einen weiteren Block rechts anfügen
		% Ausgangsgleichungen bestimmen:
        y(:,idxT) = biomodell_mgl(x0_); % erneut unnötig (s.o.): x = x0_(1:3) würde reichen
		
%% Berechnung der Fisher'schen Informationsmatrix (Fisher-Matrix)
		
		% Sensitivitäten aus Simulation		
		XP = reshape(x_(n+1:end,idxT),n,m); % aus dem erweitertern Zustandsvektor...
        % nur noch Xp extrahieren, und daraus eine (n,m)-Matrix formen
		
		% dgdx: partielle Ableitung von g nach den Zuständen
		dgdx = biomodell_dmgldx(x(:,idxT)); % Dim. (q,n)
        
        % dgdp: partielle Ableitung von g nach den Parametern (0 in unserem
		% Fall, weil die Messgleichung nicht explizit von den Parametern abhängt.)
		dgdp = zeros(size(dgdx,1),size(XP,2));  % (q,m)-Nullmatrix

		%% Aufgabe 1b: Berechnung von dgdtheta (totales Differential) und FIM
		dgdtheta = dgdx * XP + dgdp;
        FM = FM + dgdtheta' * invC * dgdtheta;      % siehe (4.46) im ausgedruckten Skript (AS)    
        
	end % for idxT
	
	% Ausgabe des Simulationsstructs
	ERGEBNIS(idxMess).t = MESS(idxMess).t;
	ERGEBNIS(idxMess).x0 = MESS(idxMess).x0;
	ERGEBNIS(idxMess).x = x;
	ERGEBNIS(idxMess).x_sim = x_sim;
	ERGEBNIS(idxMess).u = MESS(idxMess).u;
	ERGEBNIS(idxMess).y = y;
	
end % for idxMess
