function xSol = computeXTrajectories(x0,f,inputMat,tEvents,tOverall,tGrid,tEnd,thNum,cNum,aNum,nStates)
% computes state trajectories xSol in absolute coordinates starting at x0

% xSol -        complete state trajectories in abs coordinates resulting 
% from ODE f and inputMat over time tOverall starting in x0
% x0 -          initial state value
% f -           numerical function handle of ODE
% inputMat -    array of feed volume flow u and inlet concentrations at feedings xInMat [u,xInMat] 
% tEvents -     time vector of feeding events (on or off)
% tOverall -    unified time vector of time grid tGrid and feeding events
% tGrid -       time grid
% tEnd -        end time of simulation
% thNum,cNum,aNum - arrays of time-variant, time-invariant and
% stoichiometric parameters
% nStates -     # states

% Solve ODE via iterative solution of constant feeding regimes (on or off)
xSim = zeros(length(tOverall), nStates);% allocate memory
tSim = zeros(length(tOverall),1);       % allocate memory

nIntervals = length(tEvents); 
% integrate ODEs for each interval (=time when feeding constantly =on or
% =off):
tic
for cI = 1:nIntervals
    if cI == nIntervals   % letztes Intervall geht von letztem Fütterungsimpuls (aus) bis Simulationsende ...
        tCurrent   = tEvents(end);
        tNext      = tEnd;
    else    % ... alle anderen Fütterungsintervalle:
        tCurrent   = tEvents(cI);
        tNext      = tEvents(cI + 1);
    end
    
    % Get current feeding volume flow and inlet concentrations:
    inputVector = interp1(tEvents, inputMat, tCurrent, 'nearest', 0); 
    % split into current values of feedVolFlow and xIn: 
    feedVolFlowCurr = inputVector(1); 
    xInCurr = inputVector(2:end)'; 
    odeFun = @(t,x) f(x,feedVolFlowCurr,xInCurr,thNum,cNum,aNum); 
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode           = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % in this case, the solver would interpret 
        % t_ode as a time span and choose integration time points on his own
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % request to evaluate at exactly 3 time points
        [tVec, solVec] = ode15s(odeFun, t_ode, x0);
        % of the 3 time points evaluated, only save 2 (first and last):
        xSim(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode has more than two time instants, which are the times 
        % when the integration is evaluated:
        [tVec, solVec] = ode15s(odeFun, t_ode, x0); % Returns >= 3 values
        xSim(idxTimeInterval,:) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0 = solVec(end, :);    % update initial value for next interval
end
toc

% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSol = xSim(idxGrid,:);

end