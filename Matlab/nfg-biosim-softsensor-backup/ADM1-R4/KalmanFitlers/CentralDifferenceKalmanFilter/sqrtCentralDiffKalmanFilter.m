function [xPlus,PPlus,Kv] = sqrtCentralDiffKalmanFilter(xOld,POld,tSpan,feedInfo,yMeas,AC)
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

nStates = length(xOld); 
% nMeas = length(yMeas); 

%% Tuning
% Initialisierung
sigmaV = 0.2*24;    % Trommelgaszähler [L/h], umgerechnet in [L/d]
sigmaCh4 = 0.002;   % [%]
sigmaCo2 = 0.002;   % [%]
sigmaSIN = 0.12;    % [g/L]
sigmaTS = 0.01;     % Trocknungswaage
sigmaVS = 0.0035; 
sigma = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
% R = diag(sigma.^2); % Kovarianzmatrix des Messrauschens (aus Sensordatenblättern)
SR = diag(sigma);   % cholesky factor of measurement noise covariance

Q = diag([0.0644, 0.6848, 1.6561, 348.1543, 2.0443, 0.4208, 0.6206, 0.5350, 1.5051, 0.5701, 0.0912]); % ad-hoc Ansatz: 
% Diagonalmatrix aus initialem Startwert xMinus
% XY: Tuning noch anzupassen!
SQ = chol(Q, "lower");       % cholesky factor of process noise covariance

% fester Parametersatz: (Hydrolysekonstanten)
th = [0.25, 0.2, 0.1, 0.02];   % [kch, kpr, kli, kdec] 
% füge den zum struct AC hinzu: 
AC.th = th; 

%% Time Update (TU) in 3 steps

%% 1. Choose Sigma Points X of Time Update (TU)
h = sqrt(3);  % optimal value for second-order approximation (King Skript)

try 
    SxOld = chol(POld,'lower'); % van der Merwe requires lower triangular matrix
catch 
    % make sure POld is symmetric pos. definite: 
    eps = 1e-5;     % all negative eigenvalues of POld are increased up to eps
    POldspd = makePSymmPosDef(POld, eps); 
    SxOld = chol(POldspd,'lower');
end

XTU = [xOld, repmat(xOld,[1,nStates]) + h.*SxOld, ...
             repmat(xOld,[1,nStates]) - h.*SxOld]; 

%% 2. Propagate Sigma Points
% placeholder for all propagated Sigma Points (only x-component needed, not w)
XTUMinus = zeros(size(XTU)); 
% mind index shift: v.d.M. starts at 0, Matlab at 1

tEvents = feedInfo(1,:);    % feeding time points (on/off)
idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
    odeFun = @(t,X) f(t,X,feedVector,AC);
    % auf die for-Schleife kannst du nicht verzichten, denn die Matlab
    % ode-Solver kennen nur Vektoren als Startwerde, keine Matrizen!
    for ll = 1:2*nStates + 1
        [~,XTUSol] = ode15s(odeFun,tEval,XTU(:,ll));
        XTUMinus(:,ll) = XTUSol(end,:)';
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    XAtIntEnd = XTU;% Startwerte für erstes Intervall (Achtung Matrix mit allen 2n+1 Sigmapunkten)
    % 'XAtIntEnd' = X at the end of the interval
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for jj = 1:nIntervals
        feedVector = feedInfo(:,jj);
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFun = @(t,X) f(t,X,feedVector,AC);
        for lll = 1:2*nStates + 1
            [~,XTUSol] = ode15s(odeFun,tEval,XAtIntEnd(:,lll));
            XAtIntEnd(:,lll) = XTUSol(end,:)';
        end
    end
    XTUMinus = XAtIntEnd; % übernimm Werte am Ende des letzten Intervalls
end

%% 3. Aggregate Sigma Points to Priors for x and P
% create weights vector for state prior: 
wm0 = (h^2-nStates)/h^2; % XY: Achtung: negativer Wert!
wmi = 1/(2*h^2); 
wm = [wm0, repmat(wmi,1,2*nStates)];

% aggregate state prior:
xMinus = sum(wm.*XTUMinus,2);  
% xMinusTest = XTUMinus*wm';    % alternative Darstellung nach Arrese

% weights for prior of state error cov. matrix P:
wc1 = 1/(4*h^2);
wc2 = (h^2 - 1)/(4*h^2);

% aggregate this to cholesky-Factor of PMinus:
M1 = sqrt(wc1)*(XTUMinus(:,2:nStates+1) - XTUMinus(:,nStates + 2:end)); 
M2 = sqrt(wc2)*(XTUMinus(:,2:nStates+1) + XTUMinus(:,nStates + 2:end) - ...
                2*repmat(XTUMinus(:,1),[1,nStates]));
M = [M1, M2, SQ]; 
SxMinus = qr(M',"econ")'; % liefert nur den relevanten Teil von R ohne Nullzeilen.
% Transposition von M und nach qr, damit man lower triangular Matrix der Dim (n,n) erhält 

%% Measurement Update (MU) in 3 Steps (XY: wie viele sind es wirklich?): 

%% 1. Choose Sigma Points of Measurement Update

XMU = [xMinus,  repmat(xMinus,[1,nStates]) + h.*SxMinus, ...
                repmat(xMinus,[1,nStates]) - h.*SxMinus]; 

%% 2. Derive Sigma-Measurements and aggregate them
% XY: dieser Schritt ist gefährlich, denn XMU kann auch negative
% Sigmapunkte enthalten. Dann kommt bullshit in der Messgleichung raus.
Y = BMR4_AB_mgl_mat(XMU,AC.c); 

% aggregate expected measurement: 
ySim = sum(wm.*Y,2);
% hTest = Y*wm'; % alternative Darstellung nach Arrese

% aggregate this to cholesky-Factor of Pyy:
My1 = sqrt(wc1)*(Y(:,2:nStates+1) - Y(:,nStates + 2:end)); 
My2 = sqrt(wc2)*(Y(:,2:nStates+1) + Y(:,nStates + 2:end) - ...
                2*repmat(Y(:,1),[1,nStates]));
My = [My1, My2, SR]; 
Sy = qr(My',"econ")'; % liefert nur den relevanten Teil von R ohne Nullzeilen.
% Transposition von M und nach qr, damit man lower triangular Matrix der Dim (n,n) erhält 

Pxy = sqrt(wc1)*SxMinus*[Y(:,2:nStates+1) - Y(:,nStates + 2:end)]';

%% 4. Kalman Gain and actual MU:

K = (Pxy/Sy')/Sy; 
Kv = K*(yMeas - ySim);
xPlus = xMinus + Kv;

U = K*Sy;
[~,nCols] = size(U); 

% create a copy to overwrite in every iteration:
STemp = SxMinus';   % for cholupdate we need upper triangular matrices!  
for col = 1:nCols
    STempUpdate = cholupdate(STemp,U(:,col),"-")';
    STemp = STempUpdate; % overwrite STemp for next iteration
    PTempUpdate = STempUpdate'*STempUpdate; 
    lambda = eig(PTempUpdate)

    % predict the covariance matrix and chol. factor of the next iteration:
    PTempNext = PTempUpdate - U(:,col+1)*U(:,col+1)'; 
    % and compute its eigenvalues to see if it will lose its pos. def'ness:
    lambdaNext = eig(PTempNext)
    all(lambdaNext > 0)
    
    try 
        chol(PTempNext) 
    catch 
        % make sure POld is symmetric pos. definite: 
        eps = 1e-5;     % all negative eigenvalues of POld are increased up to eps
        PTempNextspd = makePSymmPosDef(PTempNext, eps); 
        chol(PTempNextspd,'lower');
    end

end

SxPlus = STempUpdate';     % turn into lower Cholesky Factor again 

% transform to PPlus: 
PPlus = SxPlus*SxPlus'; 

end
