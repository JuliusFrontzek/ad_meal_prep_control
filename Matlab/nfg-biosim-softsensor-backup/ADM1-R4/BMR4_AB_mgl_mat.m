% Messgleichung für BMR4_AB für den ganzen Zeitschrieb, also 

function y = BMR4_AB_mgl_mat(x,c)
% computes the measurement outputs y of dim [nSample,nMeas] 
% y - model output for all sampling times
% x - matrix of state trajectories as row vectors -> dim [nSample,nStates]
% c - fixed, aggregated system paramters

[nSample,~] = size(x); 

% Messgleichungen
volFlow = (c(6).*x(:,10).^2 + c(7).*x(:,10).*x(:,11) + c(8).*x(:,11).^2 + c(9).*x(:,10) + c(10).*x(:,11) + ones(nSample,1)*c(11))./24;   % [L/d] -> [L/h]
pch4 = c(12).*x(:,10);                                          % [bar]
pco2 = c(13).*x(:,11);                                          % [bar]
SIN = x(:,3);                                                   % [g/L]
TS = ones(nSample,1) - ones(nSample,1).*x(:,4)./c(14);          % [-]
VS = ones(nSample,1) - ones(nSample,1).*x(:,9)./(c(14)-x(:,4)); % [-]

y = [volFlow,pch4,pco2,SIN,TS,VS];

end
