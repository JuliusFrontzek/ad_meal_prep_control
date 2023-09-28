%% Version
% (R2022b) Update 5
% Erstelldatum: 25.04.2023
% last modified: 24.09.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = extendedKalmanFilterMultiRate(xOld,POld,feedInfo,...
                        yMeas,params,Q,R,f,g,dfdx,dhdx, ...
                        tSpan,nStates,qOn,qOff,flagAugmented,flagDelayPeriod,flagMajor)

% compute time & measurement update with Joseph-version of EKF for
% multirate measurements

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% xOld - old state estimate
% POld - old state error covariance matrix
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector (minor or major instance)
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise
% f - function handle of ODEs of system equations
% g - function handle of output equations 
% dfdx - function handle of partial derivatives of df/dx 
% dhdx - function handle of partial derivatives of df/dx 
% tSpan - time interval between old and new measurement (& state estimate)
% nStates - # states (without augmentation)
% qOn - # online measurement signals
% qOff - # offline measurement signals
% flagAugmented - 0: non-augmented, 1: augmented
% flagDelayPeriod - 0: NOT waiting for sample from lab to return. 1: waiting
% flagMajor - 0: minor instance, 1: major instance

global counter

% make print-statement if xOld contains any negative concentrations: 
if any(xOld<0)
    % print statement: 
    disp(['neg. Konzentration bei ',num2str(tSpan(1)),' Tagen'])
    xOld(xOld<0)
    % apply clipping: 
    xOld(xOld < 0) = 0; 
    counter = counter + 1
end 

%% Time Update
xPOld = [xOld;reshape(POld,[],1)];  % combined vector of xHat and P

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
    odeFun = @(t,xP) dfP_dtAug(xP,feedVolFlow,xInCurr,params,Q,f,dfdx,nStates,flagAugmented); 
    [~,xPSol] = ode15s(odeFun,tEval,xPOld);
    xPMinus = xPSol(end,:)';
% Case b: feeding changes during measurement interval:
else 
    % create fine time grid of relevant feeding events and measurements:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    xP0 = xPOld;    % initial value for first interval
    % integrate each interval sequentially to ensure constant feedings:
    for jj = 1:nIntervals
        feedVolFlow = feedInfo(jj,2); 
        xInCurr = feedInfo(jj,3:end)';   % current inlet concentrations
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFun = @(t,xP) dfP_dtAug(xP,feedVolFlow,xInCurr,params,Q,f,dfdx,nStates,flagAugmented); 
        [~,xPSol] = ode15s(odeFun,tEval,xP0);
        % update initial value for next interval:
        xP0 = xPSol(end,:)';
    end
    xPMinus = xP0;
end

% separate states x and covariance matrix P:
% differentiate between augmented and non-augmented case:
if flagAugmented == 1
    xAugMinus = xPMinus(1:2*nStates);
    PAugMinus = reshape(xPMinus(2*nStates+1:end),[2*nStates,2*nStates]);
    xMinus = xAugMinus(1:2*nStates); 
    PMinus = PAugMinus; 
else % non-augmented case:
    xMinus = xPMinus(1:nStates);
    PMinus = reshape(xPMinus(nStates+1:end),[nStates,nStates]);
end

%% Measurement Update
xMinusCore = xMinus(1:nStates); % without augmentation
h = g(xMinusCore,params.c);     % predicted model output
if flagAugmented == 1
    Ha = zeros(qOn+qOff,2*nStates);     % allocate memory
    Ha(:,1:nStates) = dhdx(xMinusCore,params.c);  % the remaining nStates cols stay 0
    H = Ha; 
else % non-augmented case:
    H = dhdx(xMinusCore,params.c);  % partial derivative of output, evaluated at xMinus
end

% Slicing: during minor instance...
if flagMajor == 0
    % ... reduce H and R to # of measurement signals during minor instance:
    H = H(1:qOn,:); 
    R = R(1:qOn,1:qOn); 
    % ... extract reduced # of measurement signals during minor instance:
    h = h(1:qOn); 
    yMeas = yMeas(1:qOn); 
end

% ...whereas during major instance, use all of the measurements available
% (merged online and offline measurements together)
S = H*PMinus*H' + R;    % auxiliary matrix
SInv = S\eye(size(R));  % efficient least squares (numerically more efficient alternative for inverse)
K = PMinus*H'*SInv;     % Kalman Gain matrix
Kv = K*(yMeas' - h);    % effective correction of Kalman Gain on state estimate (n,1); 
% Kv = zeros(nStates,1);  % XY Rania

% Achtung: only during delay period suppress measurement update of sample
% state through indicator matrix:
if flagDelayPeriod == 1
    M = blkdiag(eye(nStates),zeros(nStates)); % indicator matrix
    xPlus = xAugMinus + M*Kv;    % updated state estimation
    % Joseph-Algorithm for measurement update of P: 
    leftMat = eye(2*nStates) - M*K*H; 
    rightMat = leftMat'; 
    PPlus = leftMat*PMinus*rightMat + K*R*K'; % updated covariance of estimation error
else % outside delay period, i.e. when not waiting for lab measurement to 
    % return, update all of x and P, regardless whether augmented or not:
    xPlus = xMinus + Kv;    % updated state estimation
    % Joseph-Algorithm for measurement update of P: 
    leftMat = eye(size(xMinus)) - K*H; 
    rightMat = leftMat'; 
    PPlus = leftMat*PMinus*rightMat + K*R*K'; % updated covariance of estimation error
end

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 

end
