% DGL zu Einfachem biologischen Modell

function dxdt = biomodell_3_zdgl(t,x,u,p)

%% Grenzen der Zustände

for idxX = 1:length(x)
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);

%% Definition der Parameter

mumax = p(1);
KS = p(2);
YXS = p(3);
muMmax = p(4); % einfaches Maintenance-Modell
KM = p(5); % zusätzlicher Term - komplexes Maintenancemodell

%% Definition der Stellgrößen

tu = u(1,:);
q = u(2,find(t >= tu,1,'last'));
cSF = u(3,find(t >= tu,1,'last'));

%% Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS);
% muM = muMmax * cS ; % einfaches Maintenance-Modell
muM = muMmax * cS / (cS + KM); % Komplexes Maintenance-Modell

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q; % ohne Maintenance-Modell
dxdt(2,1) = -1/YXS * mu * mX + cSF * q - muM * mX; % mit Maintenance-Modell
dxdt(3,1) = q;
