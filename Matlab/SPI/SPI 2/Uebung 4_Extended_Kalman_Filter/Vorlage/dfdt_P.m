function dxdt = dfdt_P(t,x,u,p_KF)

%% Parameter

p = p_KF;

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

%% Definition der Stellgrößen

q = u;
cSF = 100;

%% Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS);

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q;
dxdt(3,1) = q;

end