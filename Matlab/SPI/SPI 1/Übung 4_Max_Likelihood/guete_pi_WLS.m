% G�tefunktion f�r Parameter-Identifikation
% p: Parameter
% MESS: Struct, in dem u(t), tk (Messzeitpunkte) und y(tk) (Messungen)
% gespeichert sind
% invC: Inverse der Matrix C (Kovarianzmatrix des Messrauschens)

function I = guete_pi_WLS(p,MESS,invC)

I = 0;  % Anfangswert des G�tefunktionals
for idxMess = 1:length(MESS)	% Laufindex �ber Experimente
	
%% Simulation
	
	x0 = MESS(idxMess).x0;  % initialer Anfangszustand
	y_sim = biomodell_mgl(x0);  % simulierter Ausgang am Anfangszustand
    N = length(MESS(idxMess).t); % Anzahl Messpunkte
	for idxT = 2:N    % Iteration �ber alle Messzeitpunkte
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)]; % Delta t des Intervalls
		[~,X] = ode15s(@biomodell_zdgl,tspan,x0,[],MESS(idxMess).u,p); % integriere System-DGLs f�r einen Zeitschritt hoch vorher wurde auch noch T abgespeichert neben X
		x0 = X(end,:)';     % neuer Anfangszustand f�r n�chste Iteration
		y_sim(:,idxT) = biomodell_mgl(x0);  % simulierter Ausgang am Ende des Intervalls
	end % idxT
	
%% Differenz zwischen Simulation und Messung
    
    yMeas = MESS(idxMess).y;    % Messwerte
	diff = yMeas - y_sim;       % Differenz zwischen Mess- und simuliertem Ausgang, elem R_q,N
	
%% G�tefunktional
    for k = 1:N
        I = I + diff(:,k)' * invC * diff(:,k);  
    end % for
    I = 0.5*I; 

end % for idxMess
