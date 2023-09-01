%% Version
% (R2022b) Update 5
% Erstelldatum: 30.08.2023
% Autor: Simon Hellmann

function [xPlus,PPlus] = extendedKalmanFilterMultiRateNoDelay(xOld,POld,feedInfo,yMeas,params,Q,R,f,g,dfdx,dhdx,tSpan,nStates)

% compute time & measurement update with Joseph-version of EKF for
% multirate measurements without delay (Intensiv-Beprobung offline)

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

global counter

% make print-statement if xOld contains any negative concentrations: 
if nStates > 13 && any(xOld<0) % beim ADM1-R3 kann S_ion negativ werden
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
    odeFun = @(t,xP) dfP_dt(xP,feedVolFlow,xInCurr,params,Q,f,dfdx); 
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
        odeFun = @(t,xP) dfP_dt(xP,feedVolFlow,xInCurr,params,Q,f,dfdx); 
        [~,xPSol] = ode15s(odeFun,tEval,xP0);
        % update initial value for next interval:
        xP0 = xPSol(end,:)';
    end
    xPMinus = xP0;
end

% separate states x and covariance matrix P:
xMinus = xPMinus(1:nStates);
PMinus = reshape(xPMinus(nStates+1:end),[nStates,nStates]);

%% Measurement Update
h = g(xMinus,params.c);     % full predicted model output
H = dhdx(xMinus,params.c);  % partial derivative of output, evaluated at xMinus

% figure out which measurements are available and reduce h, yMeas, H and R
% accordingly: 
idxMeas = ~isnan(yMeas);    % 1 where measurement, 0 where NaN
h = h(idxMeas);         % still col vector
yMeas = yMeas(idxMeas); % still row vector
H = H(idxMeas,:); 
R = R(idxMeas,idxMeas); 

S = H*PMinus*H' + R;    % auxiliary matrix
SInv = S\eye(size(R));  % efficient least squares (numerically more efficient alternative for inverse)
K = PMinus*H'*SInv;     % Kalman Gain matrix
Kv = K*(yMeas' - h);    % effective correction of Kalman Gain on state estimate (n,1); 
% Kv = zeros(nStates,1);  % XY Rania

xPlus = xMinus + Kv;    % updated state estimation
% Joseph-Algorithm for measurement update of P: 
leftMat = eye(nStates) - K*H; 
rightMat = leftMat'; 
PPlus = leftMat*PMinus*rightMat + K*R*K'; % updated covariance of estimation error

% clipping of negative state estimations: 
% xPlus(xPlus<0) = 0; 

end
