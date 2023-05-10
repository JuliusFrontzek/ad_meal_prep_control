% Funktion zur Optimalen Versuchsplanung
% Als Teil der OVP wird auch schon das Systemverhalten mit dem optimalen
% Stellgrößenverlauf simuliert. Es wird mit jeder Iteration von fmincon
% über die Berechnung der nichtlinearen Beschränkungen im Struct ERGEBNIS
% gespeichert (und dort mit jeder Iteration überschrieben)

% im struct VP sind alle Parameter der aktuellen Versuchsplanung enthalten:
% Zeit t, Anfangszustand x0, Stellgrößen u (samt Stellzeitpunkten und Eingangs-
% Konzentrationen), Parameter p und constraints CONS 
function ERGEBNIS = OVP(VP,FM_old,invC)

% finde Beschränkungen für x und y, wenn welche existieren, und merke dir 
% die Indizes der Vektoren (das brauchen wir später für die state und
% output constraints in cons_fcn); beachte: nur wenn hier +inf oder -inf
% angegeben wird, bleibt der jeweilige Eintrag leer! 
VP.CONS.xmin_idx = find(isfinite(VP.CONS.xmin));
VP.CONS.xmax_idx = find(isfinite(VP.CONS.xmax));
VP.CONS.ymin_idx = find(isfinite(VP.CONS.ymin));
VP.CONS.ymax_idx = find(isfinite(VP.CONS.ymax));

% initialer Stellgrößenverlauf (je eine Wert in u für jeden Messzeitpunkt):
u0 = VP.u(2,:); 
% optimiere Stellgrößenverlauf; berücksichtige die state und output
% constraints in den nichtlinearen inequalities cons_fcn:
uOpt = fmincon(@guete_ovp,u0,[],[],[],[],VP.CONS.umin,VP.CONS.umax,@cons_fcn,[]);
% beachte: uOpt muss deswegen nicht explizit Rückgabewert der Funktion
% OVP.m sein, weil in jeder Iteration von fmincon der aktuelle Zwischenwert
% des Stellgrößenverlaufs in VP.u abgespeichert wird!

%% Gütefunktion für PI (als local function)

	function I = guete_ovp(uVP)
		%      [Stellzeitpunkte; Fütterungsmengen; Eingangskonzentrationen]
		VP.u = [VP.u(1,:); uVP; VP.u(3,:)]; % so können die Fütterungsmengen 
        % im Struct VP mit jeder Iteration von fmincon überschrieben werden
		
        %% 2a: Normierung der FM und Definition des Gütekriteriums
% -- Fishermatrix

        FM = biomodell_FIM(VP,VP.p,invC) + FM_old; 
        % warum wird FM_old zur neuen Fisher-Matrix hinzuaddiert? 
        % --> weil die Fisher-Matrix ALLE bisher aufgenommenen Messwerte
        % enthalten muss, also auch die aus dem Ausgangs-Experiment (das
        % FM_old geliefert hat). Die neue Parameter-Identifikation müsste
        % ebenfalls ALLE Messwerte umfassen. Siehe Unterlagen zu dieser
        % Übung, S.10 unten.
        
        % Skalierung mit der aktuellen Schätzung der Parametervektoren:
        Lam = diag(VP.p); 
        FMTilde = Lam * FM * Lam; % normierte Fisher-Matrix;  
        
% -- A-Kriterium
        FMTildeInv = FMTilde\eye(size(FM)); 
        I = trace(FMTildeInv);    % minimiere spur von F^(-1)

	end % function guete_ovp

%% Funktion zur Berechnung der nichtlinearen Beschränkungen
% XY: warum muss das eine separate nested function sein? Kann man das nicht
% als anonymous function als Teil von OVP implementieren, sodass man sich
% den doppelten Aufruf von biomodell_FIM spart? 
	function [c,ceq] = cons_fcn(uVP)
		
		VP.u = [VP.u(1,:); uVP; VP.u(3,:)];
		
% -- Fishermatrix und Simulation
		
        % das ERGEBNIS muss für die Auswertung der constraints vorliegen:
		[FM,ERGEBNIS] = biomodell_FIM(VP,VP.p,invC);
		
        % inequality constraints: xmin < x < xmax; ymin < y < ymax alle
        % nach c(u) < 0 in Spaltenvektor umbauen (siehe doku für nonlcon in fmincon). 
        % jeder constraint ist eine Zeile im Vektor c(x).
        % Beachte, dass x und y zeitveränderlich sind, daher mit jedem 
        % Aufruf von fmincon eine neue Definition der Constraints. 
        % Durch die Verwendung von z.B. xmin_idx werden nur aktive Constraints 
        % (< inf) genutzt.
        % Die constraints sind gleich an allen Zeitpunkten, an denen u optimiert wird (hier 11 Stück). 
		c = reshape([  ...
            % xmin < x:
			VP.CONS.xmin(VP.CONS.xmin_idx) * ones(1,length(VP.t)) - ERGEBNIS.x(VP.CONS.xmin_idx,:) ; ...
            % x < xmax:
            ERGEBNIS.x(VP.CONS.xmax_idx,:) - VP.CONS.xmax(VP.CONS.xmax_idx) * ones(1,length(VP.t)) ; ...
            % ymin < y:
            VP.CONS.ymin(VP.CONS.ymin_idx) * ones(1,length(VP.t)) - ERGEBNIS.y(VP.CONS.ymin_idx,:) ; ...
            % y < ymax:
            ERGEBNIS.y(VP.CONS.ymax_idx,:) - VP.CONS.ymax(VP.CONS.ymax_idx) * ones(1,length(VP.t))], ...
            (length(VP.CONS.xmin_idx)+length(VP.CONS.xmax_idx)+length(VP.CONS.ymin_idx)+length(VP.CONS.ymax_idx))*length(VP.t),1);
		ceq = [];   % no equality constraints
		
	end % function cons_fcn

end % function OVP
