function dxdt = dfdt_P(t,x,u,Q_noise,p_KF)

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

% Alternative Modellunsicherheit
% dxdt = 0.9*dxdt;

P0 = reshape(x(4:end),[3,3]);

% % Dynamische Bestimmung des Zustandsrauschen
% dfdp = dfdp(x,p_KF);
% Q_noise = dfdp * 1e-1*eye(3) * dfdp';
% Q_noise(3,3) = 0.5;

% F = dfdx(x,p_KF);
% dPdt = F * P0 + P0 * F' + Q_noise;

% Matrix Riccati Gleichung
F = dfdx(x,p_KF); 
dPdt =  F * P0 + P0 * F' + Q_noise;
dxdt(4:12,1) = dPdt(1:9)';
end