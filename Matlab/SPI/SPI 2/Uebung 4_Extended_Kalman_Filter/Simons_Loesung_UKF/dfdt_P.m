function dxPdt = dfdt_P(t,xP,u,p,Q)

dxPdt = zeros(size(xP));

%% Parameter
mumax = p(1);
KS = p(2);
YXS = p(3);

%% separiere Zust�nde von Kovanrianzmatrix:
x = xP(1:3); 
P = reshape(xP(4:end),[3,3]);

%% Grenzen der Zust�nde
% hard reset nedativer Konzentrationen:
for idxX = 1:length(x)
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% Definition der Zust�nde

mX = x(1);
mS = x(2);
V = x(3);

%% Definition der Stellgr��en
q = u;
cSF = 100;

%% Zwischengr��en

cS = mS/V;
mu = mumax * cS / (cS + KS);

%% DGLs f�r Zust�nde:

dxPdt(1) = mu * mX;
dxPdt(2) = -1/YXS * mu * mX + cSF * q;
dxPdt(3) = q;

%% DGLs f�r Kovarianzmatrix:

% partielle Ableitungen der rechten Seite nach den Zust�nden, 
% ausgewertet am aktuellen Sch�tzwert f�r x:
F = dfdx(x,p);

dPdt = F*P + P*F' + Q; 

dxPdt(4:end) = reshape(dPdt,[],1);
end