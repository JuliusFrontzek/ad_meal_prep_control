% Gütefunktion, Kov.-Matrix des Messrauschens und Akaike-Maß für gegebene
% Parameter, Messwerte und Modell

function [I,C,AIC] = aic_ueb3(p,MESS,Modell)

C = 0;  % Initialwert für Messfehler-Kovarianzmatrix
N = 0;  % Initialwert # Messpunkte
m = length(p);  % Anzahl Parameter
q = 2;  % Anzahl Messgrößen

% Iteration über alle Experimente:
for idxMess = 1:length(MESS)
%% Simulation
	
	x0 = MESS(idxMess).x0;
	x0_ = x0;
    
    nSample = length(MESS(idxMess).t); 
	y_sim = zeros(q,nSample); 
    y_sim(:,1) = biomodell_mgl(x0);

    % Iteration über alle Messzeitpunkte:
	for idxT = 2:nSample

		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
        
        % Modellauswahl:
        if Modell == 1
            [~,X] = ode45(@biomodell_1_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        elseif Modell == 2
            [~,X] = ode45(@biomodell_2_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        elseif Modell == 3
            [~,X] = ode45(@biomodell_3_zdgl,tspan,x0,[],MESS(idxMess).u,p);
        end
        
        x0 = X(end,:)'; % neuer Startwert für nächste Iteration
		y_sim(:,idxT) = biomodell_mgl(x0);

	end % idxT
	
%% Differenz zwischen Simulation und Messung
	y_meas = MESS(idxMess).y;
	diffY = y_meas - y_sim;
	
%% Schätzung der Kovarianzmatrix
% aufaddieren über die Verschiedenen experimente
% Schätzung der Covarianzmatrix  2x2 Matrix der Messstellen
	C = C + diffY * diffY'; 
	N = N + size(diffY,2); % Anzahl der Zeitpunkte (Anzahl der Messungen in diff)
		
end % for idxMess

% C = C / N; % Teilen durch N
C = C / (N-m); % Teilen durch effektives N

%% Gütefunktion

% log(det(C)) ist negativ, die Det(C) soll jedoch minimiert werden
I = log(det(C)); 
% I = -log(det(C)); 

%% Aufgabe: 1a
AIC = N*I + 2*m;  

%% Aufgabe: 1b

    if N/m < 40
        % use modified AIC because too few measurements: 
        AIC_c = AIC + 2*m*(m+1)/(N-m-1); 
        AIC = AIC_c; 
    end
end