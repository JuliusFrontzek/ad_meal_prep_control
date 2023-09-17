% DGL zu biologischen Modell: 
% das einfache (3 Parameter) ist bereits (unvollständig) dargestellt, das 
% erweiterte (mit 4 bzw. 5 Parametern) ist anfänglich noch auskommentiert

function dxdt = biomodell_zdgl(t,x,u,p)
% liefert die rechte Seite des DGL-Systems bei fester Fütterung u und
% Parametern p

%% Grenzen der Zustände
for idxX = 1:length(x)  % negative Konzentrationen und Volumina sind nicht physikalisch. 
    % Offenbar kann der Algorithmus auch solche Ergebnisse generieren, die
    % dann natürlich Unsinn ergeben. Also harte Korrektur.
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zustände
mX = x(1);
mS = x(2);
V = x(3);

%% Definition der Parameter
muMax = p(1);
KS = p(2);
YXS = p(3);
% muMmax = p(4); % einfaches Maintenance-Modell
% KM = p(5); % zusätzlicher Term - komplexes Maintenancemodell

%% Definition der Stellgrößen (Achtung: in u steht mehr drin als nur die Stellgrößen (Fütterungsmengen), 
% sondern auch noch die Fütterungszeitpunkte und die Input-Konzentration!)

tu = u(1,:);                        % Fütterungszeitpunkte: 1. Zeile
q = u(2,find(t >= tu,1,'last'));    % Fütterungsmengen: 2. Zeile
cSF = u(3,find(t >= tu,1,'last'));  % Konzentration des Feeds: 3. Zeile

%% Zwischengrößen

cS = mS/V;                      % Substrat-Konzentration
mu = muMax * cS / (cS + KS);    % Wachstums-Kinetik Substratabbau zu Biomasse
% muM = muMmax * cS ;             % Kinetik einfaches Maintenance-Modell
% muM = muMmax * cS / (cS + KM);  % Kinetik komplexes Maintenance-Modell

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q; % ohne Maintenance-Modell
% dxdt(2,1) = -1/YXS * mu * mX + cSF * q - muM * mX; % mit Maintenance-Modell
dxdt(3,1) = q;
