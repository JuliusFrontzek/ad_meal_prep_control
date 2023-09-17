% Berechne die Ausgangs-Trajektorien der Modelle mit ihren optimalen
% Parametern (könnte auch direkt als weiterer Ausgabewert in aic_ueb3.m
% integriert werden.

function y_sim = computeOutputTrajectories(p,MESS,Modell)

q = 2;  % Anzahl Messgrößen

% Iteration über alle Experimente:
for idxMess = 1:length(MESS)
%% Simulation
	
	x0 = MESS(idxMess).x0;
    
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
	
end