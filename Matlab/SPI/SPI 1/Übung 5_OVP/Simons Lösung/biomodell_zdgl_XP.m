% Erweiterte DGL zu einfachem biologischen Modell
% liefert rechte Seiten f�r x und X, die alle in dxdt zusammengefasst werden
function dxdt = biomodell_zdgl_XP(t,x,u,p)

m = length(p);  % # Parameter 
n = m;          % # Zust�nde (hier Zufall)
%% Grenzen der Zust�nde

% harte Korrektur bei negativen Konzentrationen (unphysikalisch):
for idxX = 1:n 
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zust�nde

mX = x(1);
mS = x(2);
V = x(3);
XP = reshape(x(n+1:end),n,m); % die Eintr�ge hinter x sind X, diese werden in ...
% die korrekte Form einer (n,m)-Matrix gebracht

%% Definition der Parameter

muMax = p(1);
KS = p(2);
YXS = p(3);

%% Definition der Stellgr��en
% das ist eine entscheidende Stelle im Skript! Hier wird geschaut, wann das
% System tats�chlich gef�ttert wird, (wenn ja) wie viel und mit welcher
% Eingangskonzentration cSF
tu = u(1,:);                        % F�tterungszeiten
q = u(2,find(t >= tu,1,'last'));    % F�tterungsmengen
cSF = u(3,find(t >= tu,1,'last'));  % Eingangs-Konzentrationen

%% Zwischengr��en

cS = mS/V;                          % Substrat-Konzentration
mu = muMax * cS / (cS + KS);        % Monod-Kinetik

%% DGL

% eigentlicher Systemzustand (f) in den ersten drei Eintr�gen:
dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q;
dxdt(3,1) = q;

%% Aufgabe 1a: Vervollst�ndigen Sie die Formel f�r Xp
% partielle Ableitungen...
dfdx = biomodell_dzdgldx(x,p);
dfdp = biomodell_dzdgldp(x,p);

% ...f�r die DGL der Sensitivit�ten:
dXPdt = dfdx * XP + dfdp;

% Als weitere Eintr�ge von dxdt werden die Werte der Sensitivit�t des ...
% Zustands Xp dazu geschrieben. Das ist eigentlich eine n*m Matrix, deren 
% Eintr�ge nun untereinander geschrieben werden:


% Achtung: x ist hier nicht der reine Zustandsvektor, sondern der ...
% um X erweiterte Zustandsvektor, also length(x) = n + n*m = 3 + 3*3 = 12.
dxdt(3+[1:n*m],1) = dXPdt(1:n*m)';
