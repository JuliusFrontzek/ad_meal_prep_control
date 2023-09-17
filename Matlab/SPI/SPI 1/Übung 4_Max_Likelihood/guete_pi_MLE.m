% Gütefunktion für PI mit Maximum Likelihood Estimation:

function [I,C] = guete_pi_MLE(p,MESS)

N = 0;

% Iteration über alle Experimente:
for idxMess = 1:length(MESS)
	
%% Simulation
	
	x0 = MESS(idxMess).x0;
% 	x0_ = x0;
	y_sim = biomodell_mgl(x0);
    
    % Integriere (bei festen Parametern p) über alle Messintervalle die 
    % System-DGLn auf und berechne zu den Messzeitpunkten die 
    % Ausgangsgleichungen:
	for idxT = 2:length(MESS(idxMess).t)   
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];    % Delta t zwischen zwei Messwerten
        [T,X] = ode45(@biogasmodell_zdgl,tspan,x0,[],MESS(idxMess).u,p);
		x0 = X(end,:)';
		y_sim(:,idxT) = biogasmodell_mgl(x0);  % simulierter Systemausgang zu den Messzeitpunkten (g)
	end % idxT
	
%% Differenz zwischen Simulation und Messung
	
    yMeas = MESS(idxMess).y;    % Messwerte y
	yDiff = yMeas - y_sim;      % y - g
	
%% Schätzung der Kovarianzmatrix

	N = length(yMeas);  % # Messwerte
    m = 3;              % # Parameter
    CBiased = 0;        % initial value for estimated noise covariance matrix
    
    for k = 1:N
        % bilde iterativ die Summe für eine Abschätzung der Kovarianzmatrix
        % des Messrauschens:
        CBiased = CBiased + yDiff(:,k) * yDiff(:,k)'; 
        a=1;    % only used for debugging
    end % for
    
    % normiere mit der Anzahl der Freiheitsgrade:
    C = 1/(N-m) * CBiased;
		
end % for idxMess

%% Gütefunktional MLE für den Fall einer geschätzten Kovarianzmatrix des Messrauschens C:

I = N/2*log(det(C)); 

