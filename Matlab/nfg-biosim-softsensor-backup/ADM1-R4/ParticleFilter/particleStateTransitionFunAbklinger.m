function particlesProp = particleStateTransitionFunAbklinger(particlesOld,tSpan,inputVector,params,sigma)
% compute the propagated particles during measurement update based on the
% old estimate xOld
% tSpan - current measurement interval
% inputVector - [tFeed, feedVolFlow, xIn']

% global counterX 

    % separate struct params
    modParams = params.modParams; 
    physParams = params.phyParams;
    
    [nParticles,nStates] = size(particlesOld); 
    particlesProp = zeros(size(particlesOld));  % allocate memory
    
    tEval = linspace(tSpan(1),tSpan(2),3); 
    
    % propagate every particle by one time step: 
    odeFun = @(t,x) BMR4_AB_ess_ode(t,x,physParams,inputVector,modParams);
    
    for k = 1:nParticles
        x0 = particlesOld(k,:); 
        [~,xSol] = ode15s(odeFun,tEval,x0);
        xProp = xSol(end,:); % we're only interested in the result of the last integration step     
        processNoise = mvnrnd(zeros(1,nStates),sigma,1);    % create a random number from multivariate normal distribution with mean and sigma (matrix) 
        kQ = 0.01;  % XY: sicherlich viel zu geringer Wert!
        particlesProp(k,:) = xProp + kQ*processNoise; 
    end

%     counterX = counterX + 1
end