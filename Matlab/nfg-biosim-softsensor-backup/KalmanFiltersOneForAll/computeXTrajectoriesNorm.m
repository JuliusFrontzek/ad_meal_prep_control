function xSolNorm = computeXTrajectoriesNorm(x0Norm,fNorm,inputMat,tEvents,tOverall,tGrid,tEnd,thNum,cNum,aNum,nStates,TxNum,TuNum)
% computes state trajectories xSolNorm in absolute coordinates starting at 
% x0Norm

% xSolNorm -        complete state trajectories in normalized coordinates
% resulting from ODE fNorm and inputMat over time tOverall starting in x0Norm
% x0Norm -      normalized initial state value
% f -           numerical function handle of ODE
% inputMat -    array of feed volume flow u and inlet concentrations at feedings xInMat [u,xInMat] 
% tEvents -     time vector of feeding events (on or off)
% tOverall -    unified time vector of time grid tGrid and feeding events
% tGrid -       time grid
% tEnd -        end time of simulation
% thNum,cNum,aNum - arrays of time-variant, time-invariant and
% stoichiometric parameters
% nStates -     # states
% TxNum -       normalization vector of states
% TuNum -       normalization vector of inputs

% Solve ODE via iterative solution of constant feeding regimes (on or off)
xSimNorm = zeros(length(tOverall), nStates);% allocate memory
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

    % apply normalization:
    uCurrNorm = feedVolFlowCurr./TuNum; 
    xInCurrNorm = xInCurr./TxNum;

    % normalized function handle:
    odeFunNorm = @(t,xNorm) fNorm(xNorm,uCurrNorm,xInCurrNorm,thNum,cNum,aNum,TxNum,TuNum); 
    
    % Construct time vector for ODE (t_ode) by filtering of tOverall:
    idxTimeInterval = (tOverall >= tCurrent & tOverall <= tNext);
    t_ode           = tOverall(idxTimeInterval); 
    if length(t_ode) == 2   % in this case, the solver would interpret 
        % t_ode as a time span and choose integration time points on his own
        t_ode   = linspace(t_ode(1), t_ode(end), 3);    % request to evaluate at exactly 3 time points
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm);
        % of the 3 time points evaluated, only save 2 (first and last):
        xSimNorm(idxTimeInterval,:) = solVec([1 end], :);  
        tSim(idxTimeInterval) = tVec([1 end]);
    else    % t_ode has more than two time instants, which are the times 
        % when the integration is evaluated:
        [tVec, solVec] = ode15s(odeFunNorm, t_ode, x0Norm); % Returns >= 3 values
        xSimNorm(idxTimeInterval,:) = solVec;
        tSim(idxTimeInterval) = tVec;
    end
    
    x0Norm = solVec(end, :);    % update initial value for next interval
end
toc

% Evaluate xSol only in tGrid, discard the rest
idxGrid = ismember(tOverall, tGrid); 
xSolNorm = xSimNorm(idxGrid,:);

end