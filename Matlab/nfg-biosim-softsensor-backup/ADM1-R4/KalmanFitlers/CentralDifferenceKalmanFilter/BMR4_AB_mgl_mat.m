% Messgleichung für BMR4_AB für den ganzen Zeitschrieb, also 
% x - Matrizen von dim = (n,nSamples)

function y = BMR4_AB_mgl_mat(x,c)
% computes the measurement outputs y of dim [ny,nSample] 
% y - model output for all sampling times
% x - matrix of state trajectories as row vectors -> dim [nStates,nSample]
% c - fixed, aggregated system paramters

[~,nSample] = size(x); 

% Messgleichungen
volFlow = (c(9).*x(10,:).^2 + c(10).*x(10,:).*x(11,:) + c(11).*x(11,:).^2 + c(12).*x(10,:) + c(13).*x(11,:) + c(14))./24;   % [L/d] -> [L/s]
pch4 = c(15).*x(10,:);                                          % [bar]
pco2 = c(16).*x(11,:);                                          % [bar]
SIN = x(3,:);                                                   % [g/L]
TS = ones(1,nSample) - ones(1,nSample).*x(4,:)./c(17);          % [-]
VS = ones(1,nSample) - ones(1,nSample).*x(9,:)./(c(17)-x(4,:)); % [-]

y = [volFlow;
     pch4; 
     pco2; 
     SIN; 
     TS; 
     VS];

end
