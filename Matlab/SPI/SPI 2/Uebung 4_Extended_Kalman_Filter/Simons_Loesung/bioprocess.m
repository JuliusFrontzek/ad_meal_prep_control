function [x_real, y_noise] = bioprocess(tspan,x0,u)
        
        % add noise to initial values (to account for model uncertainty)
%         nStates = length(x0); 
        addedNoise = 0.03.*x0.*randn(size(x0)); % 2% of current value
        x0 = x0 + addedNoise;
        
        [t,x] = ode45(@(t,x) bioprocess_dxdt(t,x,u),tspan,x0);
        x_real = x(end,:)';
        
        % Messung
%         y_noise(1,1) = x_real(1)./x_real(3) + 1.15^2 * randn(1);
        y_noise(1,1) = x_real(1)./x_real(3) + min([1.15^2,0.05*x_real(1)./x_real(3)])*randn(1);
        y_noise(1,2) = x_real(3) + 0.25^2*randn(1);
        y_noise = y_noise';
	
end


%% Systemgleichung

function dxdt = bioprocess_dxdt(t,x,u)
%% Parameter

p = [0.1,0.2,0.6];

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
