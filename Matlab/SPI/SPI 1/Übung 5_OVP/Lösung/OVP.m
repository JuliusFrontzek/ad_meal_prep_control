% Funktion zur Optimalen Versuchsplanung

function ERGEBNIS = OVP(VP,FM_old,invC)

VP.CONS.xmin_idx = find(isfinite(VP.CONS.xmin));
VP.CONS.xmax_idx = find(isfinite(VP.CONS.xmax));
VP.CONS.ymin_idx = find(isfinite(VP.CONS.ymin));
VP.CONS.ymax_idx = find(isfinite(VP.CONS.ymax));

u0 = VP.u(2,:);
uOpt = fmincon(@guete_ovp,u0,[],[],[],[],VP.CONS.umin,VP.CONS.umax,@cons_fcn,[]);

%% Gütefunktion für PI

	function I = guete_ovp(uVP)
		
		VP.u = [VP.u(1,:); uVP; VP.u(3,:)];
		
% -- Fishermatrix

        FM = biomodell_FIM(VP,VP.p,invC) + FM_old;
        % Optional: Normierung/Skalierung
		FM = diag(VP.p) * FM * diag(VP.p);
        
% -- A-Kriterium
		B=eye(length(VP.p));
        invFM=FM\B;
		EW = eig(invFM);
        % Summe der Eigenwerte
		I = sum(EW);
%         % Alternativ: Summe der Diagonalelemente
%         I = sum(diag(invFM));

        % mod E
%         I = max(EW)/min(EW);

        % D-Kriterium
%         I = det(invFM);
		
	end % function guete_ovp

%% Funktion zur Berechnung der nichtlinearen Beschränkungen

	function [c,ceq] = cons_fcn(uVP)
		
		VP.u = [VP.u(1,:); uVP; VP.u(3,:)];
		
% -- Fishermatrix und Simulation
		
		[FM,ERGEBNIS] = biomodell_FIM(VP,VP.p,invC);
		
		c = reshape([  ...
			VP.CONS.xmin(VP.CONS.xmin_idx) * ones(1,length(VP.t)) - ERGEBNIS.x(VP.CONS.xmin_idx,:) ; ...
			ERGEBNIS.x(VP.CONS.xmax_idx,:) - VP.CONS.xmax(VP.CONS.xmax_idx) * ones(1,length(VP.t)) ; ...
			VP.CONS.ymin(VP.CONS.ymin_idx) * ones(1,length(VP.t)) - ERGEBNIS.y(VP.CONS.ymin_idx,:) ; ...
			ERGEBNIS.y(VP.CONS.ymax_idx,:) - VP.CONS.ymax(VP.CONS.ymax_idx) * ones(1,length(VP.t)) ] ...
			,(length(VP.CONS.xmin_idx)+length(VP.CONS.xmax_idx)+length(VP.CONS.ymin_idx)+length(VP.CONS.ymax_idx))*length(VP.t),1);
		ceq = [];
		
	end % function cons_fcn

end % function OVP
