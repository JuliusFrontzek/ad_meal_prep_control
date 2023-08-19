%% Version
% (R2022b) Update 2
% Erstelldatum: Ende 2022
% Autor: Simon Hellmann

function [ERGEBNIS] = biogasmodell_lsg_zdglOnOffline(MESS,p,pFix)
% integriere das Systemverhalten auf bei gegebenen Parametern und
% Fütterungsplan
% Achtung: aus MESS brauchen wir nicht die Messdaten, sondern nur die
% Fütterungen

nExp = length(MESS); 
% Iteration über alle Experimente (später in der OVP wird in MESS mehr als
% ein einziges Experiment enthalten sein; dann erst wird diese Iteration
% relevant):
for idxExp = 1:nExp

    [nStates,nEval] = size(MESS(idxExp).xSim); 
    
    % placeholders: 
    xSim = zeros(nEval,nStates);
    tSim = zeros(nEval,1);   
    
	% Anfangsbedingungen für DGLn
	x0 = MESS(idxExp).x0;       % Anfangszustand
%     xSimFine = [];   % XY: trage dessen Berechnung noch nach, um die
%     Simulation besser plotten zu können! Du könntest dann eine zweite
%     Integration mit von ode ausgewählten Auswerte-Zeitpunkten hinzufügen
%     x = x0;
    
    % time vectors:
    tEvents = MESS(idxExp).u(:,1);  % time instances, when feeding is turned on/off
    tGridOn = MESS(idxExp).tOn;     % time grid of online measurements
    tGridOff = MESS(idxExp).tOff;   % time grid of offline measurements
    tOverall = unique([tGridOn; tGridOff; tEvents]); % Join and sort timestamps
   
    nIntervals = length(tEvents); 
    
    inputMat = MESS(idxExp).u(:,2:end);    % [vol flow,inlet concentrations]
    
    % iterate over all periods of constant feeding (on or off)
    for cI = 1:nIntervals
        if cI == nIntervals   % letztes Intervall geht von letztem Fütterungsimpuls (aus) bis Simulationsende ...
            tCurrent   = tEvents(end);
            tNext      = tGridOn(end);    % end of simulation time
        else    % ... alle anderen Fütterungsintervalle:
            tCurrent   = tEvents(cI);
            tNext      = tEvents(cI + 1);
        end

        % Get current feeding volume flow and inlet concentrations:
        inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
        modelOdeFun = @(t, x) biogasmodell_zdgl(t, x, p, inputVector, pFix);

        % Construct time vector for ODE (t_ode) by filtering of tOverall:
        idxTimeInterval  = (tOverall >= tCurrent & tOverall <= tNext);
        t_ode       = tOverall(idxTimeInterval); 
        if length(t_ode) == 2   % the solver would interpret t_ode as a time 
        % span otherwise in this case. but we want precisely 3 evaluations:
            t_ode   = linspace(t_ode(1), t_ode(end), 3);    % set 3 equally-spaced integration points
            [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0);
            xSim(idxTimeInterval,:) = solVec([1 end], :);  
            tSim(idxTimeInterval) = tVec([1 end]);
        else    % t_ode has more than 2 elements. These are the evaluation points for ode15s to integrate between 
            [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0); % Returns >= 3 values
            xSim(idxTimeInterval, :) = solVec;
            tSim(idxTimeInterval) = tVec;
        end % if
    
        x0 = solVec(end, :); % update initial value for next interval
    end % for
    
    % Evaluate xSim either in the online or offline time grid, discard the rest
    idxGridOn = ismember(tOverall, tGridOn); 
    idxGridOff = ismember(tOverall, tGridOff); 
    xSolOn = xSim(idxGridOn, :);
    xSolOff = xSim(idxGridOff, :);
    
    % evaluate system output for online- and offline measurements
    yOnTemp = biogasmodell_mgl(xSolOn,pFix);
    yOffTemp = biogasmodell_mgl(xSolOff,pFix);
    yOn = yOnTemp(:,1:3); 
    yOff = yOffTemp(:,4:6); 
     
   	% Ausgabe des Simulationsstructs
	ERGEBNIS(idxExp).tOn = MESS(idxExp).tOn;
    ERGEBNIS(idxExp).tOff = MESS(idxExp).tOff;
    ERGEBNIS(idxExp).tSim = MESS(idxExp).tSim;
	ERGEBNIS(idxExp).x0 = MESS(idxExp).x0;  % save last state initial value for next simulation
% 	ERGEBNIS(idxExp).x = x;
	ERGEBNIS(idxExp).xSim = xSim;
	ERGEBNIS(idxExp).u = MESS(idxExp).u;
	ERGEBNIS(idxExp).yOn = yOn;
    ERGEBNIS(idxExp).yOff = yOff;
	
end % for idxExp

