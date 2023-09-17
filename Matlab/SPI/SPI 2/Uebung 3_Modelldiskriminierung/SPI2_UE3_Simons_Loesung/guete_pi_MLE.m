% Gütefunktion für PI
% max. Likelihood-Schätzung ohne bekanntes Messrauschen C
function [I,C] = guete_pi_MLE(p,MESS,Modell)

C = 0;
N = 0;
for idxMess = 1:length(MESS)
	
%% Simulation
	
	x0 = MESS(idxMess).x0;
	x0_ = x0;
	y_sim = biomodell_mgl(x0);
	for idxT = 2:length(MESS(idxMess).t)
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
        
        if Modell == 1
            [T,X] = ode45(@biomodell_1_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        elseif Modell == 2
            [T,X] = ode45(@biomodell_2_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        elseif Modell == 3
            [T,X] = ode45(@biomodell_3_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        end
        
        x0 = X(end,:)';
		y_sim(:,idxT) = biomodell_mgl(x0);
	end % idxT
	
%% Differenz zwischen Simulation und Messung
	
	diff = MESS(idxMess).y - y_sim;
	
%% Schätzung der Kovarianzmatrix
% aufaddieren über die Verschiedenen experimente
% Schätzung der Covarianzmatrix  2x2 Matrix der Messstellen
	C = C + diff * diff'; 
	N = N + size(diff,2); % Anzahl der Zeitpunkte (Anzahl der Messungen in diff)
		
end % for idxMess

% C = C / N; % Teilen durch N
m = length(p);  % number of parameters
C = C / (N-m); % Teilen durch effektives N

%% Gütefunktion

% log(det(C)) ist negativ, die Det(C) soll jedoch minimiert werden
I = log(det(C)); 

