%% Version
% (R2022b) Update 5
% Erstelldatum: 28.02.23
% Autor: Simon Hellmann

% Gütefunktion für PI mit Maximum Likelihood Estimation:
function [I,COn,COff] = myGuete_pi_MLEOnOffline(pEst,MESS,pFix)

nExp = length(MESS); 

% Iteration über alle Experimente:
for idxExp = 1:nExp
    [nStates,nEval] = size(MESS(idxExp).xSim); 

    % placeholders: 
    xSim = zeros(nEval,nStates);
    tSim = zeros(nEval,1);   
    
	% initial conditions for ODEs
	x0 = MESS(idxExp).x0;       % Anfangszustand
%     x_sim = [];   % XY: trage dessen Berechnung noch nach, um die
%     Simulation besser plotten zu können! Du könntest dann eine zweite
%     Integration mit von ode ausgewählten Auswerte-Zeitpunkten hinzufügen
%     x = x0;
    
    % time vectors:
    tEvents     = MESS(idxExp).u(:,1);  % time instances, when feeding is turned on/off
    tGridOn     = MESS(idxExp).tOn;     % time grid of online measurements
    tGridOff    = MESS(idxExp).tOff;    % time grid of offline measurements
    tOverall    = unique([tGridOn; tGridOff; tEvents]); % combine time grids
    
    nIntervals  = length(tEvents);
    inputMat    = MESS(idxExp).u(:,2:end);    % [feed vol flow, inlet concentrations]
    
    for cI = 1:nIntervals
        if cI == nIntervals   % last interval goes from last feeding impuls (off) till end of simulation ...
            tCurrent   = tEvents(end);
            tNext      = tGridOn(end);    % end of simulation time
        else    % ... all other feeding intervals:
            tCurrent   = tEvents(cI);
            tNext      = tEvents(cI + 1);
        end

        % Get current feeding volume flow and inlet concentrations:
        inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
        modelOdeFun = @(t, x) biogasmodell_zdgl(t, x, pEst, inputVector, pFix);

        % derive interval time instances t_ode from overall time grid, such: 
        % that they only include periods of constant feeding (constantly on/off):
        idxTimeInterval  = (tOverall >= tCurrent & tOverall <= tNext);
        t_ode = tOverall(idxTimeInterval); 
        
        % make sure there is always >2 time instances, so t_ode is NOT
        % interpreted as a time span:
        if length(t_ode) == 2   
            t_ode = linspace(t_ode(1), t_ode(end), 3);    % set 3 equally-spaced integration points
            [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0);
            xSim(idxTimeInterval,:) = solVec([1 end], :);  
            tSim(idxTimeInterval) = tVec([1 end]);
        else    % t_ode includes more than two time instances already. Use them as integration time point: 
            [tVec, solVec] = ode15s(modelOdeFun, t_ode, x0); % Returns >= 3 values
            xSim(idxTimeInterval, :) = solVec;
            tSim(idxTimeInterval) = tVec;
        end % if
    
        x0 = solVec(end, :); % update initial value for next interval
    end % for
    
    % Evaluate xSol either in the online or offline time grid, discard the rest
    idxGridOn   = ismember(tOverall, tGridOn); 
    idxGridOff  = ismember(tOverall, tGridOff); 
    xSolOn      = xSim(idxGridOn, :);
    xSolOff     = xSim(idxGridOff, :);
    
    % evaluate system output for online- and offline measurements
%     ySimOn = biogasmodell_mgl_on(xSolOn',pFix);
%     ySimOff = biogasmodell_mgl_off(xSolOff',pFix);
    ySimOnTemp  = biogasmodell_mgl( xSolOn,pFix); 
    ySimOffTemp = biogasmodell_mgl(xSolOff,pFix); 
    ySimOn      = ySimOnTemp( :,1:3); 
    ySimOff     = ySimOffTemp(:,4:6);

    %% Differenz zwischen Simulation und Messung
    yMeasOn     = MESS(idxExp).yMeasOn;    % online Messwerte y
    yMeasOff    = MESS(idxExp).yMeasOff;   % offline Messwerte y
%     yDiffOn = yMeasOn - ySimOn;         % nicht normalisierte Differenz online
%     yDiffOff = yMeasOff - ySimOff;      % nicht normalisierte Differenz offline
    
    % normalize measurements & simulations with max measurement values: 
    maxYMeasOn  = max(yMeasOn);     % max online
    maxYMeasOff = max(yMeasOff);   % max offline
%     minYMeasOn = min(yMeasOn,[],2);     % min online
%     minYMeasOff = min(yMeasOff,[],2);   % min offline
    
    % normalized simulations:
    ySimOnNorm  = ySimOn ./ maxYMeasOn;   % online
    ySimOffNorm = ySimOff ./ maxYMeasOff; % offline
    % normalized measurements:
    yMeasOnNorm = yMeasOn ./ maxYMeasOn;   % online
    yMeasOffNorm = yMeasOff./ maxYMeasOff; % offline
    
	yDiffOnNorm = yMeasOnNorm - ySimOnNorm;   % yOnline - gOnline
    yDiffOffNorm = yMeasOffNorm - ySimOffNorm;   % yOffline - gOffnline
    
%     % überschreibe die Einträge für SIN durch Normalisierung auf [0,1]:
%     % Achtung: liefert schlechtere Fits!
%     SINSim01 = (ySimOff(1,:) - minYMeasOff(1))./(maxYMeasOff(1) - minYMeasOff(1)); 
%     SINMeas01 = (yMeasOff(1,:) - minYMeasOff(1))./(maxYMeasOff(1) - minYMeasOff(1));
%     yDiffOffNorm(1,:) = SINMeas01 - SINSim01; 
    
    %% Schätzung der Kovarianzmatrizen (eine für online, eine für offline Messungen)
	NOn  = length(yMeasOn);     % # online Messwerte
    NOff = length(yMeasOff);    % # offline Messwerte
    N    = NOn + NOff;          % # of all samples
    m    = length(pEst);        % # of Parameters
    
    COn  = 1/(NOn-m) * (yDiffOnNorm'*yDiffOnNorm); 
    COff = 1/(NOff-m) * (yDiffOffNorm'*yDiffOffNorm); 
    
end % for idxMess

%% Gütefunktional MLE für den Fall einer geschätzten Kovarianzmatrix des Messrauschens C:
% qOn = length(MESS(1).yMeasOn(1,:));      % # of online measurement signals 
% qOff = length(MESS(1).yMeasOff(1,:));    % # of offline measurement signals 
IOn     = N/2*log(det(COn));  % ignore additional summand NOn*qOn/2*log(2*pi)
IOff    = N/2*log(det(COff)); % ignore additional summand NOff*qOff/2*log(2*pi)
I       = IOn + IOff; 

% show evolution of cost function summands throughout the iterations:
disp(['IOn: ',num2str(IOn), ' IOff: ',num2str(IOff)])
end % function
