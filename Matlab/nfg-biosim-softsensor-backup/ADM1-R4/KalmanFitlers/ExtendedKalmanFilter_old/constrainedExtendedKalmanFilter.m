%% Version
% (R2022b) Update 2
% Erstelldatum: Januar 2023
% Autor: Simon Hellmann

function [xPlus,PPlus,Kv] = constrainedExtendedKalmanFilter(xOld,POld,tSpan,feedInfo,yMeas,AC,R)
% time and measurement update acc. to the central difference Kalman Filter
% xPlus - new state estimate
% PPlus - new state error covariance matrix
% Kv - effective Kalman Gains (same dimension/units as states)
% xOld - old state estimate
% POld - old state error covariance matrix
% tSpan - time interval between old and new measurement (& state estimate)
% feedInfo - combination of feeding information [tEvents; feedVolFlow; inlet concentrations]
% yMeas - latest measurement vector
% AC - struct with stoichiometric coefficients and aggregated constant
% model parameters
% R - covariance matrix of measurement noise

%% Tuning
nStates = length(xOld); 
nOutputs = length(yMeas);
% Q = diag([0.0644, 0.6848, 1.6561, 348.1543, 2.0443, 0.4208, 0.6206, 0.5350, 1.5051, 0.5701, 0.0912]);  
Q = diag([0.016, 0.555, 0.563, 958.4, 1.263, 2.816, 2.654, 0.972, 2.894, 10, 0.374, 0.948]);
% Diagonalmatrix aus initialem Startwert x0SS (ad-hoc Ansatz)
% XY: noch anzupassen!

% fester Parametersatz: (Hydrolysekonstanten)
th = [3, 0.25, 0.2, 0.1, 0.02, 0.75]; % [kchF, kchS, kpr, kli, kdec, fracChFast] 
% füge den zum struct AC hinzu: 
AC.th = th; 

%% Time Update
xPOld = [xOld;reshape(POld,[],1)];  % combined vector of xHat and P

tEvents = feedInfo(1,:);    % Fütterungszeitpunkte (an/aus)
% alt: 
% idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
% tRelEvents = tEvents(idxRelEvents);   % Auswertung anhand Index

% neu: 
filterRelEvents = tEvents >= tSpan(1) & tEvents <= tSpan(2); 
tRelEvents = tEvents(filterRelEvents); % Auswertung anhand bool-Wertes

% integriere über Intervalle konstanter Fütterung: 
% Fall a: konstante Fütterung während ges. Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
    odeFun = @(t,xP) dfdt_P(t,xP,feedVector,AC,Q);
    [~,xPSol] = ode15s(odeFun,tEval,xPOld);
    xPMinus = xPSol(end,:)';
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus relevanten Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    xP0 = xPOld;    % Startwert für erstes Intervall
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen
    for jj = 1:nIntervals
        feedVector = feedInfo(:,jj);
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFun = @(t,xP) dfdt_P(t,xP,feedVector,AC,Q);
        [~,xPSol] = ode15s(odeFun,tEval,xP0);
        % Update Startwert für nächstes Intervall:
        xP0 = xPSol(end,:)';
    end
    xPMinus = xP0;
end

% separiere Zustände und Kovarianzmatrix:
xMinus = xPMinus(1:nStates);
PMinus = reshape(xPMinus(nStates+1:end),[nStates,nStates]);

%% Measurement Update
H = dhdx(xMinus,AC.c);  % partial derivative of output, evaluated at xMinus
S = H*PMinus*H' + R;    % auxiliary matrix
SInv = S\eye(nOutputs); % efficient least squares solution of matrix inverse
K = PMinus*H'*SInv;     % Kalman Gain matrix

h = BMR4_AB_frac_mgl(xMinus,AC.c);  % predicted model output
Kv = K*(yMeas - h);     % effective correction of Kalman Gain on state estimate (n,1); 

% alte Version aus dem EKF:
% xPlus = xMinus + Kv;    % updated state estimation
% PPlus = PMinus - K*H*PMinus; % updated covariance of estimation error

% Neu: Berücksichtige für das measurement update, dass die Zustände positiv
% bleiben müssen. Löse also ein beschränktes Optimierungsproblem: 
costFun = @(x) evaluateCEKFCostFun(x,yMeas,R,xMinus,PMinus,AC.c); 
A = -eye(nStates); 
b = zeros(nStates,1);  
xPlus = fmincon(costFun,xMinus,A,b); 

% Joseph-Form of update of state error covariance matrix: 
PPlus = (eye(nStates) - K*H)*PMinus*(eye(nStates) - K*H)' + K*R*K';

% hard correction of negative state estimations: 
% xPlus(xPlus<0) = 0; 
end




