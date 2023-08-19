%% Version
% (R2022b) Update 2
% Erstelldatum: November 2022
% Autor: Simon Hellmann

function dxPdt = dfdt_P(t,xP,feedVector,AC,Q)
% dxPdt - right-hand side of ODE (of states and dynamics of state error
% covariance matrix)
% t - time vector
% xP - states and covariance matrix (reshaped as vector), stacked above
% each other as column vector
% feedVector - combination of [feedVolFlow; inlet concentrations]
% AC - struct with stoichiometric coefficients and aggregated constant
% model parameters
% Q - power spectral density matrix of process noise

nStates = 12; 

% extract constant parameters out of struct: 
th = AC.th; 
c = AC.c; 
a = AC.a;

dxPdt = zeros(size(xP));    % placeholder

% separiere Zustände und Kovanrianzmatrix:
x = xP(1:nStates); 
P = reshape(xP(nStates+1:end),[nStates,nStates]);

% hard reset negativer Konzentrationen:
% x(x < 0) = 0; 

% Stellgröße und Eingangskonzentrationen:
u = feedVector(1);  % feedVolumeFlow
xi = feedVector(2:end); % Note: für x11, x12 (Gasphase) gibt es prinzipiell 
% keine Eingangskonzentrationen, daher hat xi auch nur 9 Einträge

%% DGLs für Zustände:

f = BMR4_AB_frac_ode(t,x,feedVector,AC);

dxPdt(1:nStates) = f; 

%% DGLs für Kovarianzmatrix (dynamics of state error covariance matrix):

% partielle Ableitungen der rechten Seite nach den Zuständen, ausgewertet 
% am aktuellen Schätzwert für x:
F = dfdx(x,u,th,c,a);

dPdt = F*P + P*F' + Q;      % dynamics of state error covariance matrix

dxPdt(nStates+1:end) = reshape(dPdt,[],1);
end