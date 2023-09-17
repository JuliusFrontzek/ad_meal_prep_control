function dfdx = dfdx(x,p_KF)


%% Zustände
mX = x(1);
mS = x(2);
V = x(3);

%% Parameter
p = p_KF;
mumax = p(1);
KS = p(2);
YXS = p(3);

%% Hilfsgrößen
cS = mS/V;
mu = mumax * cS / (cS + KS);

dmudcS = mumax * KS / (cS+KS)^2;


dfdx = [mu        dmudcS/V*mX           -dmudcS*mS/V^2*mX
        -1/YXS*mu -1/YXS*dmudcS/V*mX    1/YXS*dmudcS*mS/V^2*mX
	    0         0                     0];
end