% DGL mit der M�glichkeit f�r den JacGen Zust�nde, Stellgr��en und Parameter
% symbolisch definiert zu werden

function dxdt = Bsp_biomodell(t,x,u,p,varargin)

%% Grenzen der Zust�nde



%% Symbolische Definition der Zust�nde und Parameter f�r den JacGen

if isempty(t)
    syms('mumax','KS','YXS','mX','mS','V','q','cSF');
else
    
for idxX = 1:3
if x(idxX) < 0
    x(idxX) = 0;
end % if
end % for idxX

%% Definition der Zust�nde & Parameter mit Werten falls aufgerufen durch ode-solver
mX = x(1);
mS = x(2);
V = x(3);

% Definition der Parameter

mumax = p(1);
KS = p(2);
YXS = p(3);

%% Definition der Stellgr��en

tu = u(1,:);
q = u(2,find(t >= tu,1,'last'));
cSF = u(3,find(t >= tu,1,'last'));
end
%% Zwischengr��en

cS = mS/V;
mu = mumax * cS / (cS + KS);

%% DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q;
dxdt(3,1) = q;