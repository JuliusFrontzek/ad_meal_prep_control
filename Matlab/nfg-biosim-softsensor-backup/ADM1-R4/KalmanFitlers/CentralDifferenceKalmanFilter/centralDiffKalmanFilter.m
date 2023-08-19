function [xPlus,PPlus,Kv] = centralDiffKalmanFilter(xOld,POld,tSpan,feedInfo,yMeas,AC)
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

disp(['t: ', num2str(tSpan(1))]);

nStates = length(xOld); 
nMeas = length(yMeas); 

%% Tuning
% Initialisierung
sigmaV = 0.2*24;    % Trommelgaszähler [L/h], umgerechnet in [L/d]
sigmaCh4 = 0.002;   % [%]
sigmaCo2 = 0.002;   % [%]
sigmaSIN = 0.12;    % [g/L]
sigmaTS = 0.01;     % Trocknungswaage
sigmaVS = 0.0035; 
sigma = [sigmaV, sigmaCh4, sigmaCo2, sigmaSIN, sigmaTS, sigmaVS]; 
R = diag(sigma.^2); % Kovarianzmatrix des Messrauschens (aus Sensordatenblättern)

Q = diag([0.0644, 0.6848, 1.6561, 348.1543, 2.0443, 0.4208, 0.6206, 0.5350, 1.5051, 0.5701, 0.0912]); % ad-hoc Ansatz: 
% Diagonalmatrix aus initialem Startwert xMinus
% XY: noch anzupassen!

% fester Parametersatz: (Hydrolysekonstanten)
th = [0.25, 0.2, 0.1, 0.02];   % [kch, kpr, kli, kdec] 
% füge den zum struct AC hinzu: 
AC.th = th; 

%% Time Update (TU) in 4 steps

%% 1. State Augmentation (with process noise w)
xOldW = [xOld;zeros(nStates,1)];    % augmented state vector
POldW = zeros(2*nStates);           % placeholder 
POldW(1:nStates,1:nStates) = POld; 
POldW(nStates+1:end,nStates+1:end) = Q; % augmented state error cov. matrix

%% 2. Choose Sigma Points X of Time Update
% XY: durch Cholesky-Zerlegung hier potentielle Fehlerquelle!
h = sqrt(3);  % optimal value for second-order approximation (King Skript)
XTU0 = xOldW; 
% insert checking whether POldW is positive definite:  
lambdaTU = eig(POldW);
flagPosDefTU = all(lambdaTU > 0); 

try 
    cholPOldW = chol(POldW);    % Cholesky Factorization
catch 
    % make sure POldW is symmetric pos. definite
    eps = 1e-3;     % all negative eigenvalues of POld are increased up to eps
    POldWspd = makePSymmPosDefWhile(POldW, eps); 
    cholPOldW = chol(POldWspd);  % cholesky factorization
end

% placeholder for all Sigma Points of Time Update
XTU = zeros(2*nStates, 4*nStates + 1); 
% mind index shift: King starts at 0, Matlab at 1:
XTU(:,1) = XTU0; 

% fill remaining columns in two iterations
% 1:2n
for k = 2 : 2*nStates + 1
    XTU(:,k) = xOldW + h*cholPOldW(:,k-1);
end 
% 2n+1:4n
for kk = 2*nStates+2 : 4*nStates + 1
    XTU(:,kk) = xOldW - h*cholPOldW(:,kk-2*nStates-1);
end

% partition Sigma Points XTU = [XTUx; Xw]:
XTUx = XTU(1:nStates,:); 
Xw = XTU(nStates+1:end,:); 

%% 3. Propagate Sigma Points
% placeholder for all propagated Sigma Points (only x-component needed, not w)
XTUxMinus = zeros(nStates, 4*nStates + 1); % mind index shift: King starts at 0, Matlab at 1
% XTU(:,1) = XTU0;

tEvents = feedInfo(1,:);    % feeding time points (on/off)
idxRelEvents = find(tEvents >= tSpan(1) & tEvents <= tSpan(2));
tRelEvents = tEvents(idxRelEvents);

% we can only perform integration when feeding is constant!
% Fall a: konst. Fütterung während gesamtem Messintervalls (keine Änderung)
if isempty(tRelEvents)
    feedVector = feedInfo; % verwende die aktuell wirksame Fütterung
    tEval = tSpan;
    odeFun = @(t,X) f(t,X,feedVector,AC);
    for ll = 1:4*nStates + 1
        [~,XTUxSol] = ode15s(odeFun,tEval,XTUx(:,ll));
        XTUxMinus(:,ll) = XTUxSol(end,:)' + Xw(:,ll);
    end 
% Fall b: veränderliche Fütterung während Messintervalls:
else 
    % erstelle 'feines' Zeitraster aus Fütterungs- und Messzeiten:
    tOverall = unique(sort([tSpan, tRelEvents]));
    nIntervals = length(tOverall) - 1; 
    X0 = XTUx;    % Startwerte für erstes Intervall (Achtung Matrix mit allen 4n+1 Sigmapunkten)
    % integriere Intervall-weise, sodass während der Intervalle konstante
    % Fütterungen herrschen:
    for jj = 1:nIntervals
        feedVector = feedInfo(:,jj);
        tEval = [tOverall(jj), tOverall(jj+1)];
        odeFun = @(t,X) f(t,X,feedVector,AC);
        for lll = 1:4*nStates + 1
            [~,XTUxSol] = ode15s(odeFun,tEval,X0(:,lll));
            X0(:,lll) = XTUxSol(end,:)' + Xw(:,lll);
        end
    end
    XTUxMinus = X0;
end

%% 4. Aggregate Sigma Points to Priors for x and P
% create weights vector for state prior: 
wTU0 = (h^2-nStates)/h^2; 
wTUi = 1/(2*h^2); 
wmTU = [wTU0, repmat(wTUi,1,4*nStates)];

% weights for prior of state error cov. matrix P:
wc1 = 1/(4*h^2);
wc2 = (h^2 - 1)/(4*h^2);

% aggregate state prior:
xMinus = sum(wmTU.*XTUxMinus,2);   

% aggregate prior of state error cov. matrix P:
PMinus = zeros(nStates);    % placeholder
% Achtung Index-Shift! (Matlab startet bei 1, King bei 0)
for ii = 1:2*nStates
    PMinus = PMinus + ...
    wc1.*(XTUxMinus(:,ii+1) - XTUxMinus(:,2*nStates+ii+1))*(XTUxMinus(:,ii+1) - XTUxMinus(:,2*nStates+ii+1))' + ...
    wc2.*(XTUxMinus(:,ii+1) + XTUxMinus(:,2*nStates+ii+1) - 2*XTUxMinus(:,1))*...
         (XTUxMinus(:,ii+1) + XTUxMinus(:,2*nStates+ii+1) - 2*XTUxMinus(:,1))';
end
% XY: geht das vielleicht auch direkt mit Vektor-Operationen statt mit der
% for-Schleife? 

%% Measurement Update (MU) in 4 Steps: 

%% 1. State Augmentation (with measurement noise v)
xMinusV = [xMinus;zeros(nMeas,1)];  % augmented state vector
PMinusV = zeros(nStates+nMeas);     % placeholder 
PMinusV(1:nStates,1:nStates) = PMinus; 
PMinusV(nStates+1:end,nStates+1:end) = R; % augmented state error cov. matrix

% hard correction of negative entries in P's diagonal:
% XY: checke VL-Inhalte zum Joseph-Algorithmus (abseits des harten Resets)!
% diagPMinus = diag(PMinusVTemp); 
% diagPMinus(diagPMinus < 0) = 1e-3;    % we want positive definite matrices, so 0 is unacceptable!
% % remove diag from PMinusVTemp...
% PMinusVNoDiag = PMinusVTemp - diag(diag(PMinusVTemp)); 
% % ... and add the corrected diag instead: 
% PMinusV = PMinusVNoDiag + diag(diagPMinus);

%% 2. Choose Sigma Points of Measurement Update
XMU0 = xMinusV; 
% insert checking whether PMinusV is positive definite: 
lambdaMU = eig(PMinusV); 
flagPosDefMU = all(lambdaMU > 0); 

try 
    cholPMinusV = chol(PMinusV);    % Cholesky Factorization
catch 
    % make sure POldW is symmetric pos. definite: 
    PMinusVspd = makePSymmPosDef(PMinusV, eps); 
    cholPMinusV = chol(PMinusVspd); 
end

% placeholder for all Sigma Points of Measurement Update:
XMU = zeros(nStates+nMeas, 2*(nStates+nMeas)+1); 
% mind index shift: King starts at 0, Matlab at 1:
XMU(:,1) = XMU0; 

% fill remaining columns in two iterations
% 1:n+q
for k = 2 : nStates+nMeas+1
    XMU(:,k) = xMinusV + h*cholPMinusV(:,k-1);
end 
% n+q+1 : 2(n+q)
for kk = nStates+nMeas+2 : 2*(nStates+nMeas)+1
    XMU(:,kk) = xMinusV - h*cholPMinusV(:,kk-(nStates+nMeas)-1);
end

% partition Sigma Points XMU = [XMUx; Xv]:
XMUx = XMU(1:nStates,:); 
Xv = XMU(nStates+1:end,:); 

%% 3. Derive Sigma-Measurements and aggregate them

Y = BMR4_AB_mgl_mat(XMUx,AC.c) + Xv; 

% we don't need all 4n+1 entries of wmTU for the measurement update anymore:
wmMU = wmTU(1:2*(nStates+nMeas)+1); 
% aggregate expected measurement: 
h = sum(wmMU.*Y,2);

% (auto-)covariance of measurement: 
Pyy = zeros(nMeas);    % placeholder
% beachte Index-Shift Matlab<>King:
for ii = 1:nStates+nMeas
    Pyy = Pyy + ...
    wc1.*(Y(:,ii+1) - Y(:,nStates+nMeas+ii+1))*(Y(:,ii+1) - Y(:,nStates+nMeas+ii+1))' + ...
    wc2.*(Y(:,ii+1) + Y(:,nStates+nMeas+ii+1) - 2*Y(:,1))*...
         (Y(:,ii+1) + Y(:,nStates+nMeas+ii+1) - 2*Y(:,1))';
end
% XY: geht das vielleicht auch ohne for-Schleife und dafür mit schlauer
% Vektor-Matrix-Multiplikation?

% cross covariance matrix states/measurements:
% insert checking whether PMinus is positive definite: 
lambdaPMinus = eig(PMinus); 
flagPosDefPMinus = all(lambdaPMinus > 0);

try 
    cholPMinus = chol(PMinus);    % Cholesky Factorization
catch 
    % make sure POldW is symmetric pos. definite: 
    PMinusSpd = makePSymmPosDef(PMinus, eps); 
    cholPMinus = chol(PMinusSpd); 
end

% beachte Index Shift Matlab<>King:
Pxy = sqrt(wc1)*cholPMinus*(Y(:,2:nStates+1) - ...
                            Y(:,nStates+nMeas+2:2*nStates+nMeas+1))';

%% 4. Kalman Gain and actual MU:
K = Pxy*inv(Pyy); 
Kv = K*(yMeas - h);

xPlusTemp = xMinus + Kv; 
PPlus = PMinus - K*Pyy*K'; 

% hard correction of negative state estimations: 
xPlus = xPlusTemp; 
xPlus(xPlus<0) = 0; 

% hard correction of negative entries in P's diagonal:
% XY: checke VL-Inhalte zum Joseph-Algorithmus (abseits des harten Resets)!
% diagP = diag(PPlusTemp); 
% diagP(diagP < 0) = 1e-2;    % we want positive definite matrices, so 0 is unacceptable!
% % remove diagP from PPlus...
% PPlusNoDiag = PPlusTemp - diag(diag(PPlusTemp)); 
% % ... and add the corrected diag instead: 
% PPlus = PPlusNoDiag + diag(diagP); 

%% check if PPlus is pos. def! 
% insert checking whether PMinusV is positive definite: 
lambdaPPlus = eig(PPlus); 
flagPosDefPPlus = all(lambdaPPlus > 0); 

end
