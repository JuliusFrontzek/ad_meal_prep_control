function xProp = stateTransitionFunAbklinger(xOld,tSpan,inputVector,params)
% compute the propagated state xPorp during measurement update based on the
% old estimate xOld
% tSpan - current measurement interval
% inputVector - [tFeed, feedVolFlow, xIn']

% global counterX 

    % separate struct params
    modParams = params.modParams; 
    physParams = params.phyParams;
    
    % XY: hier immer Fütterung mit dem gleichem regime einbauen:
    tEval = linspace(tSpan(1),tSpan(2),3);  
    odeFun = @(t,x) BMR4_AB_ess_ode(t,x,physParams,inputVector,modParams);
    % compute state propagation: 
    [~,xSol] = ode15s(odeFun,tEval,xOld);
    xProp = xSol(end,:); % we're only interested in the result of the last integration step
%     counterX = counterX + 1
end