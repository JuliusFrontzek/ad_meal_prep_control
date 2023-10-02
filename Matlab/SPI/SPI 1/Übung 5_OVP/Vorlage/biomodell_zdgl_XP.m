% Erweiterte DGL zu einfachem biologischen Modell

function dxdt = biomodell_zdgl_XP(t,x,u,p)

%% Grenzen der Zustände

for idxX = 1:3
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);
XP = reshape(x(4:end),3,3);

%% Definition der Parameter

mumax = p(1);
KS = p(2);
YXS = p(3);

%% Definition der Stellgrößen

tu = u(1,:);
q = u(2,find(t >= tu,1,'last'));
cSF = u(3,find(t >= tu,1,'last'));

%% Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS);

dfdx = biomodell_dzdgldx(x,p);
dfdp = biomodell_dzdgldp(x,p);

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q;
dxdt(3,1) = q;

%% Aufgabe 1a: Vervollständigen Sie die Formel für Xp
dXPdt = ;
dxdt(3+[1:length(p)^2],1) = dXPdt(1:length(p)^2)';
