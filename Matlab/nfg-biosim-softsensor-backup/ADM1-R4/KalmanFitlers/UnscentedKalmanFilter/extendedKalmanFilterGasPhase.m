%% Version
% (R2022b) Update 2
% Erstelldatum: 25.04.2023
% Autor: Simon Hellmann

function [xPlus,PPlus,Kv] = extendedKalmanFilterGasPhase(xOld,POld,tSpan,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx)
% compute time and measurement update acc. to Joseph-version of the EKF
% only use pCH4 and gasVolFlow for gas phase estimation (no pCO2 anymore)

% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% params - struct with stoichiometric coefficients a, aggregated constants
% c and time-variant parameters th
% Q - power spectral density matrix of process noise
% R - covariance matrix of measurement noise
% f - function handle of ODEs of system equations
% g - function handle of output equations 
% dfdx - function handle of partial derivatives of df/dx 
% dhdx - function handle of partial derivatives of df/dx 

global counterX

% make print-statement if xOld contains any negative concentrations: 
if any(xOld<0)
    % print statement: 
    disp(['neg. Konzentration bei ',num2str(tSpan(1)),' Tagen'])
    xOld(xOld<0)
    % apply clipping: 
    xOld(xOld < 0) = 0; 
    counterX = counterX + 1
end 

nStates = length(xOld); 

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
    odeFun = @(t,xP) dfP_dt(xP,nStates,feedVolFlow,xInCurr,params,Q,f,dfdx);    % XY Abhängigkeiten korrigieren!
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
        odeFun = @(t,xP) dfP_dt(xP,nStates,feedVolFlow,xInCurr,params,Q,f,dfdx);    % XY Abhängigkeiten korrigieren!
        [~,xPSol] = ode15s(odeFun,tEval,xP0);
        % update initial value for next interval:
        xP0 = xPSol(end,:)';
    end
    xPMinus = xP0;
end

% separate states and covariance matrix:
xMinus = xPMinus(1:nStates);
PMinus = reshape(xPMinus(nStates+1:end),[nStates,nStates]);

%% Measurement Update
h = g(xMinus,params.c); % predicted model output
H = dhdx(xMinus,params.c);  % partial derivative of output, evaluated at xMinus
% exclude CO2 measurement to avoid to overdefine the gas production: 
h = h([1,2,4:end]); 
H = H([1,2,4:end],:); 
R = R([1,2,4:end],[1,2,4:end]); 
yMeas = yMeas([1,2,4:end]); 

S = H*PMinus*H' + R;    % auxiliary matrix
SInv = S\eye(size(R));  % numerically more efficient way to compute inverse
K = PMinus*H'*SInv;     % Kalman Gain matrix

Kv = K*(yMeas' - h);    % effective correction of Kalman Gain on state estimate (n,1); 
% Kv = zeros(nStates,1);  % XY Rania

xPlus = xMinus + Kv;    % updated state estimation

% apply Joseph-Algorithm: 
leftMat = eye(nStates) - K*H; 
rightMat = leftMat'; 
PPlus = leftMat*PMinus*rightMat + K*R*K'; % updated covariance of estimation error

% old: without Joseph Algorithm:
% PPlus = PMinus - K*H*PMinus; % updated covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 
end
