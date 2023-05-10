% Gütefunktion für PI

function I = guete_pi_WLS(p,MESS,invC)

I = 0;
for idxMess = 1:length(MESS)	% Laufindex über Experimente
	
%% Simulation

options = odeset('RelTol',1e-3,'AbsTol',1e-6);

	x0 = MESS(idxMess).x0;
	y_sim = biomodell_mgl(x0);
	for idxT = 2:length(MESS(idxMess).t)
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
		[T,X] = ode15s(@biomodell_zdgl,tspan,x0,options,MESS(idxMess).u,p);
		x0 = X(end,:)';
		y_sim(:,idxT) = biomodell_mgl(x0);
	end % idxT
	
%% Differenz zwischen Simulation und Messung
	
	diff = MESS(idxMess).y - y_sim; % ergibt eine MAtrix 2x12 für Messpunkte x Zeitintervalle
	
%% Gütefunktion
	
% multiplizieren mit 2x2 der Covarianz der Messung und aufsummieren über N
    for i = 1:length(diff)
        I = I + diff(:,i)' * invC * diff(:,i);
    end

end % for idxMess
