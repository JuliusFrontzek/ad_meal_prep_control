	function [FM,EW,invFM] = guete_ovp_FIM(VP,uVP,FM_old,invC)
		% berechne die Fisher-Matrix, die nach einem neuen Experiment mit 
        % dem Stellgrößenverlauf aus der OVP resultiert
    
        % überschreibe den Stellgrößenverlauf mit den aus der OVP
		% erhaltenen Werten für u (tu und cSF bleiben gleich):
        VP.u = [VP.u(1,:); uVP; VP.u(3,:)];
		
%% -- Fishermatrix
		FM = biomodell_FIM(VP,VP.p,invC) + FM_old;
		% beachte: die OVP liefert ein neues Experiment. Die Fisher-Matrix
		% ist zu berechnen unter Hinzunahme dieses neuen Experiments zu
		% ihrem vorherigen Wert (ohne die OVP)
        
%% -- A-Kriterium
		% Bilden der Inverse:
        B = eye(length(VP.p));
		%invFM = inv(FM);
        invFM=FM\B;
        % Achtung: die MuLö hat hier nicht das A-Kriterium implementiert,
        % sondern miniert die Summe der Eigenwerte von F^-1, maximiert also 
        % die Summe der EW von F, minimiert also die Summe der Hauptachsen.
        % Das entspricht einer Abwandlung des D-/E-Kriterium auf alle EW.
		EW = eig(invFM);
		I = sum(EW);