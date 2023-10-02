% Originalmodell

function MESS = biomodell_messdaten(MESS)
% erzeugt neue Messwerte (MESS), basierend auf einigen vorgegebenen Feldern
% im Struct MESS (t, x0 & u), der als Eingang übergeben wird

t = MESS.t;
x0 = MESS.x0;
u = MESS.u; 

% Initialisierung
x = x0;     % in x werden nur die Zustände zu den Messzeitpunkte gespeichert
x_sim = []; % Platzhalter für Simulationsergebnisse

% Itegriere über alle Messintervalle, d.h. simuliere das "echte" (=gestörte) Systemverhalten:
for idxT = 2:length(t)
	tspan = [t(idxT-1) t(idxT)];
	[T,X] = ode45(@biomodell_orig_zdgl,tspan,x0,[],u);  % Originalmodell aufintegrieren während tspan
	x_sim = [x_sim [T';X']];    % erweitere x_sim mit jedem Messintervall um 
    % eine Blockmatrix (erste Zeile: Simulationszeitpunkte, 2.-4. Zeile:
    % Zustände zu Simulationszeitpunkten
    
	% neuer Anfangswert = alter Endwert:
    x0 = X(end,:)'; 
	x(:,idxT) = x0;
end % idxT

% Fülle die restlichen Felder des Structs MESS mit den neu ermittelten
% Daten aus obiger Iteration:
MESS.x = x;
MESS.x_sim = x_sim;
MESS.y = biomodell_mgl(MESS.x);
% Messwerte werden mit normalverteiltem Rauschen überlagert:
MESS.y = MESS.y + sqrt(0.05) .* randn(size(MESS.y,1),length(t));

%% DGL (als nested function)

function dxdt = biomodell_orig_zdgl(t,x,u)

%% -- Grenzen der Zustände

% keine negativen Konzentrationen erlauben:
for idxX = 1:length(x)
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% -- Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);

%% -- Definition der Parameter

mumax = 0.1;
KS = 0.2;
YXS = 0.6;
muMmax = 0.4; % maximale Wachstumsrate
KM = 30;
%% -- Definition der Stellgrößen

tu = u(1,:);    % Fütterungszeitpunkte
q = u(2,find(t >= tu,1,'last'));
cSF = u(3,find(t >= tu,1,'last'));

%% -- Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS); % Wachstum (Biomasse aus Substrat)
muM = muMmax * cS / (cS + KM); % Erhaltungsterm als Monod-Kinetik (M: Maintenance)
% muM = muMmax * cS ; % Erhaltungsterm als lineare Kinetik 
%% -- DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q 	- muM * mX;
dxdt(3,1) = q;
