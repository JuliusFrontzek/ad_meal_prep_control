function xProp = stateTransitionFunction(xOld,tSpan,fullInputMat,params)
% compute the propagated state xPorp during measurement update based on the
% old estimate xOld
% tSpan - current measurement interval
% fullInputMat - [tFeed, feedVolFlow, xIn']
  
    % separate struct params
    modParams = params.modParams; 
    physParams = params.phyParams;
    
    % find feeding times that matter during the current measurement
    % interval:
    tEventsFromStart = fullInputMat(:,1);    % feeding time points (on/off) plus start
    feedVolFlows = fullInputMat(:,2);
%     tFeed = tEventsFromStart(2); 
    % find relevant feed changes (without zero feed in beginning):
    idxFristFeedChange = find(feedVolFlows > 0, 1); 
    tFeedChanges = tEventsFromStart(idxFristFeedChange:end);

    % see which feed changes happen during tSpan: 
    idxRelFeedChanges = tSpan(1) <= tFeedChanges & tFeedChanges <= tSpan(2); 
    tRelFeedChanges = tFeedChanges(idxRelFeedChanges);

    % alt: 
%     idxRelFeedChanges = find(tSpan(1) <= tFeedChanges & tFeedChanges <= tSpan(2)); 
%     tRelFeedChanges = tFeedChanges(idxRelFeedChanges);
    
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort(round([tSpan; tRelFeedChanges],4)));

    % we can only perform integration when feeding is constant!
    % XY: überprüfe, wie du hier interp1 einbauen kannst! 
    % Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
    if length(tOverall) == 2% isempty(tRelFeedChanges)
        % verwende interp1, um das aktuell anliegende Fütterungsregime [volFlow, inletConcentrations]
        % zu finden: 
        feedRegime = interp1(tEventsFromStart,fullInputMat(:,2:end),tOverall(1),'previous');
        % XY: überprüfe, ob das immer korrekte Ergebnisse liefert!
        
        % verwende die aktuell wirksame Fütterung. es sind aber nur die
        % Mengen und Eingangskonzentrationen relevant:
%         feedVector = fullInputMat(:,2:end); 
        % only evaluate integration at three defined points:
        tEval = linspace(tSpan(1),tSpan(2),3);  
%         odeFun = @(t,X) f(t,X,feedVector,AC);
        odeFun = @(t,x) BMR4_AB_ode(t,x,physParams,feedRegime,modParams);
        % compute state propagation: 
        [~,xSol] = ode15s(odeFun,tEval,xOld);
        xProp = xSol(end,:); % we're only interested in the result of the last integration step
    
    % Fall b: veränderliche Fütterung während Messintervalls:
    else 
%         % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
%         tOverall = unique(sort([tSpan; tRelFeedChanges]));
        nIntervals = length(tOverall) - 1; 
        x0 = xOld;    % Startwerte für erstes Intervall

        % Fallunterscheidung: fallen Fütterungszeitpunkte und
        % Messzeitpunkte zusammen oder sind sie unterschiedlich?
        
        % relevant feeding regime during intervals is determined by feeding 
        % at previous time point:
        feedRegime = interp1(tEventsFromStart,fullInputMat(:,2:end),tOverall(1:end-1),'previous');
        % XY: überprüfe, ob das immer korrekte Ergebnisse liefert!

        % integriere Intervall-weise, sodass während der Intervalle konstante
        % Fütterungen herrschen:
        for k = 1:nIntervals
%             feedVector = fullInputMat(jj,2:end);
            % evaluate integration only at three points: 
            tEval = linspace(tOverall(k), tOverall(k+1), 3);
            % BMR4_AB_ode(t,x,physParams,input,modParams)
%             odeFun = @(t,x) f(t,x,feedVector,AC);
            odeFun = @(t,x) BMR4_AB_ode(t,x,physParams,feedRegime(k,:),modParams);
            % compute state propagation: 
            [~,xSol] = ode15s(odeFun,tEval,x0);
            % overwrite initial value for next interval:
            x0 = xSol(end,:); % we're only interested in the result of the last integration step
        end
        xProp = x0;
    end
end