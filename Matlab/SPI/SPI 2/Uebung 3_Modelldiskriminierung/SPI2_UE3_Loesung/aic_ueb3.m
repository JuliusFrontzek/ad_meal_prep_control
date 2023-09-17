% Gütefunktion für PI

function [I,C,AIC] = aic_ueb3(p,MESS,Modell)

C = 0;
N = 0;
m = length(p);
q = 2; % Anzahl der Messgrößen
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
C = C / (N-m); % Teilen durch N

%% Gütefunktion

% log(det(C)) ist negativ, die Det(C) soll jedoch minimiert werden
I = log(det(C)); 
% I = -log(det(C)); 

%% Aufgabe: 1a
AIC = -2*N*I  + 2*length(p);
% AIC = N*I  + 2*length(p);

%% Aufgabe: 1b
% für N/m < 40
if N/m <40
    AIC = AIC + (2*m*(m+1))/(N-m-1);
end