% DGL zu biologischen Modell: 
% das einfache (3 Parameter) ist bereits (unvollst�ndig) dargestellt, das 
% erweiterte (mit 4 bzw. 5 Parametern) ist anf�nglich noch auskommentiert

function dxdt = biomodell_zdgl(t,x,u,p)
% liefert die rechte Seite des DGL-Systems bei fester F�tterung u und
% Parametern p

%% Grenzen der Zust�nde
for idxX = 1:length(x)  % negative Konzentrationen und Volumina sind nicht physikalisch. 
    % Offenbar kann der Algorithmus auch solche Ergebnisse generieren, die
    % dann nat�rlich Unsinn ergeben. Also harte Korrektur.
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zust�nde
mX = x(1);
mS = x(2);
V = x(3);

%% Definition der Parameter
muMax = p(1);
KS = p(2);
YXS = p(3);
% muMmax = p(4); % einfaches Maintenance-Modell
% KM = p(5); % zus�tzlicher Term - komplexes Maintenancemodell

%% Definition der Stellgr��en (Achtung: in u steht mehr drin als nur die Stellgr��en (F�tterungsmengen), 
% sondern auch noch die F�tterungszeitpunkte und die Input-Konzentration!)

tu = u(1,:);                        % F�tterungszeitpunkte: 1. Zeile
q = u(2,find(t >= tu,1,'last'));    % F�tterungsmengen: 2. Zeile
cSF = u(3,find(t >= tu,1,'last'));  % Konzentration des Feeds: 3. Zeile

%% Zwischengr��en

cS = mS/V;                      % Substrat-Konzentration
mu = muMax * cS / (cS + KS);    % Wachstums-Kinetik Substratabbau zu Biomasse
% muM = muMmax * cS ;             % Kinetik einfaches Maintenance-Modell
% muM = muMmax * cS / (cS + KM);  % Kinetik komplexes Maintenance-Modell

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q; % ohne Maintenance-Modell
% dxdt(2,1) = -1/YXS * mu * mX + cSF * q - muM * mX; % mit Maintenance-Modell
dxdt(3,1) = q;
