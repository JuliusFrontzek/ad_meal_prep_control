% Ableitung der DGL nach den Parametern

function dfdp = biomodell_dzdgldp(x,p)

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

dmudKS = - mumax * cS / (cS+KS)^2;

%% 

dfdp = [cS/(cS+KS)*mX           dmudKS*mX           0
        -1/YXS*cS/(cS+KS)*mX    -1/YXS*dmudKS*mX    1/YXS^2*mu*mX
	    0                       0                   0];
