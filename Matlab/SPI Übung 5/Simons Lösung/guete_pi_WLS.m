% G�tefunktion f�r Parameteridentifikation

function I = guete_pi_WLS(p,MESS,invC)

I = 0;  % Initialisierung
for idxMess = 1:length(MESS)	% Laufindex �ber Experimente
	
%% Simulation

options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    
    % simuliere das Systemverhalten �ber alle Messzeitpunkte (je in den
    % Zeitintervallen zwischen den Messungen):
	x0 = MESS(idxMess).x0;      % Startwert
	y_sim = biomodell_mgl(x0);
    
    % integriere die System-DGLs �ber alle Mess-Intervalle des aktuellen 
    % Experiments und berechne zu den Messzeitpunkten den Systemausgang: 
    for idxT = 2:length(MESS(idxMess).t)
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
        % hier wird ode15s in der alten Methodik aufgerufen, n�mlich gem��
        % [t,y] = ode15s(odefun,tspan,y0,options,[feste Variable 1 nach t,x],[feste Variable 2 nach t,x])
		[T,X] = ode15s(@biomodell_zdgl,tspan,x0,options,MESS(idxMess).u,p);
		x0 = X(end,:)';
		y_sim(:,idxT) = biomodell_mgl(x0);
	end % idxT
	
%% Differenz zwischen Simulation und Messung
	
	diff = MESS(idxMess).y - y_sim; % ergibt eine MAtrix 2x12 f�r Messpunkte x Zeitintervalle
	
%% G�tefunktion
	
% multiplizieren mit 2x2 der Covarianz der Messung und aufsummieren �ber N
    for i = 1:length(diff)
        % addiere iterativ die Kostenfunktion auf f�r alle Messintervalle:
        I = I + diff(:,i)' * invC * diff(:,i);  
    end

end % for idxMess
