% integrate system dynamics from time k to time k+1

% XY: argumente und Beschreibung anpassen! 

function xPlusNorm = stateTransitionFcnADM1R4_norm(xMinusNorm,tSpan,feedInfoNorm,params,fNorm,TxNum,TuNum)
% Julius' alte argumente der Funktion: (xMinus, u, tSpan, p)

% extract constant parameters out of struct: 
th = params.th; 
c = params.c; 
a = params.a;

tEvents = feedInfoNorm(:,1);    % feeding time points (an/aus)
% find relevant feeding events (during current measurement interval):
filterRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2); 
tRelEvents = tEvents(filterRelEvents); 

% integrate across intervals with constant feeding:
% Case a: constant feeding during measurement interval:
if isempty(tRelEvents)
    feedVolFlowNorm = feedInfoNorm(2); 
    xInCurrNorm = feedInfoNorm(3:end)';  % current inlet concentrations
    tEval = tSpan; 
    odeFunNorm = @(t,xNorm) fNorm(xNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
    [~,xSolNorm] = ode15s(odeFunNorm,tEval,xMinusNorm);
    xPlusNorm = xSolNorm(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    x0Norm = xMinusNorm;    % initial value for first interval
    % integrate each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlowNorm = feedInfoNorm(jj,2); 
        xInCurrNorm = feedInfoNorm(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFunNorm = @(t,xNorm) fNorm(xNorm,feedVolFlowNorm,xInCurrNorm,th,c,a,TxNum,TuNum); 
        [~,xSolNorm] = ode15s(odeFunNorm,tEval,x0Norm);
        % update initial value for next interval:
        x0Norm = xSolNorm(end,:)';
    end % for
    xPlusNorm = x0Norm;
end % if

end % function