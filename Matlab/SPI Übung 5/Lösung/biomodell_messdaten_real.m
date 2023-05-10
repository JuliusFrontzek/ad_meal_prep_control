% Originalmodell

function MESS = biomodell_messdaten(MESS)

t = MESS.t;
x0 = MESS.x0;
x = x0;
x_sim = [];
u = MESS.u;

for idxT = 2:length(t)
	tspan = [t(idxT-1) t(idxT)];
	[T,X] = ode45(@biomodell_orig_zdgl,tspan,x0,[],u);
	x_sim = [x_sim [T';X']];
	x0 = X(end,:)';
	x(:,idxT) = x0;
end % idxT

MESS.x = x;
MESS.x_sim = x_sim;
MESS.y = biomodell_mgl(MESS.x);
MESS.y = MESS.y + sqrt(0.05) .* randn(size(MESS.y,1),length(t));

%% DGL

function dxdt = biomodell_orig_zdgl(t,x,u)

%% -- Grenzen der Zustände

for idxX = 1:length(x)
	if x(idxX) < 0
		x(idxX) = 0;
	end % if
end % for idxX

%% -- Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);

%% -- Definition der Parameter

mumax = 0.1;
KS = 0.2;
YXS = 0.6;
muMmax = 0.4; %maximaler erhaltungsterm
KM=30;
%% -- Definition der Stellgrößen

tu = u(1,:);
q = u(2,find(t >= tu,1,'last'));
cSF = u(3,find(t >= tu,1,'last'));

%% -- Zwischengrößen

cS = mS/V;
mu = mumax * cS / (cS + KS); % Wachstum
muM = muMmax * cS / (cS + KM); % Erhaltung
% muM = muMmax * cS ; % Erhaltung
%% -- DGL

dxdt(1,1) = mu * mX;
dxdt(2,1) = -1/YXS * mu * mX + cSF * q 	- muM * mX;
dxdt(3,1) = q;
