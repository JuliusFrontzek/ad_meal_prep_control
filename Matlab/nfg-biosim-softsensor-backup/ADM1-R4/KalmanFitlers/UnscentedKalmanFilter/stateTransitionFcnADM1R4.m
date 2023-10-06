% integrate system dynamics from time k to time k+1

% XY: argumente und Beschreibung anpassen! 

function xPlus = stateTransitionFcnADM1R4(xMinus,tSpan,feedInfo,params,f)
% Julius' alte argumente der Funktion: (xMinus, u, tSpan, p)

% extract constant parameters out of struct: 
th = params.th; 
c = params.c; 
a = params.a;

tEvents = feedInfo(:,1);    % feeding time points (an/aus)
% find relevant feeding events (during current measurement interval):
filterRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2); 
tRelEvents = tEvents(filterRelEvents); 

% integrate across intervals with constant feeding:
% Case a: constant feeding during measurement interval:
if isempty(tRelEvents)
    feedVolFlow = feedInfo(2); 
    xInCurr = feedInfo(3:end)';  % current inlet concentrations
    tEval = tSpan;
    % XY: hier muss dx/dt zurückgegeben werden! 
    odeFun = @(t,x) f(x,feedVolFlow,xInCurr,th,c,a); 
%         dfP_dt(xP,feedVolFlow,xInCurr,params,Q,f,dfdx);    % XY Abhängigkeiten korrigieren!
    [~,xSol] = ode15s(odeFun,tEval,xMinus);
    xPlus = xSol(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    x0 = xMinus;    % initial value for first interval
    % integrate each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlow = feedInfo(jj,2); 
        xInCurr = feedInfo(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFun = @(t,x) f(x,feedVolFlow,xInCurr,th,c,a); 
        %             odeFun = @(t,x) dfP_dt(xP,feedVolFlow,xInCurr,params,Q,f,dfdx);    % XY Abhängigkeiten korrigieren!
        [~,xSol] = ode15s(odeFun,tEval,x0);
        % update initial value for next interval:
        x0 = xSol(end,:)';
    end % for
    xPlus = x0;
end % if

end % function