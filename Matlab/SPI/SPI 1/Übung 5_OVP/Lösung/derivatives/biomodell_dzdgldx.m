% Jacobimatrix des einfachen biologischen Modells

function dfdx = biomodell_dzdgldx(x,p)

%% Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);

%% Definition der Parameter

mumax = p(1);
KS = p(2);
YXS = p(3);

%% Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS);

dmudcS = mumax * KS / (cS+KS)^2;

%% Jacobimatrix

dfdx = [
	       mu        dmudcS/V*mX      -dmudcS*mS/V^2*mX
	-1/YXS*mu -1/YXS*dmudcS/V*mX 1/YXS*dmudcS*mS/V^2*mX
	        0                  0                      0
	];
