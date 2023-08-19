% Messgleichung für BMR4_AB_frac für den ganzen Zeitschrieb, also 
% x - Matrizen von dim = (n,nSamples)

function y = BMR4_AB_frac_mgl_mat(x,c)
% computes the measurement outputs y of dim [ny,nSample] 
% y - model output for all sampling times
% x - matrix of state trajectories as row vectors -> dim [nStates,nSample]
% c - fixed, aggregated system paramters

[nSample,~] = size(x); 

% Messgleichungen
volFlow = (c(9).*x(:,11).^2 + c(10).*x(:,11).*x(:,12) + c(11).*x(:,12).^2 + c(12).*x(:,11) + c(13).*x(:,12) + c(14))./24;  % L/d --> L/h
pch4 = c(15).*x(:,11); 
pco2 = c(16).*x(:,12); 
SIN = x(:,3); 
TS = ones(nSample,1) - ones(nSample,1).*x(:,4)./c(17);
VS = ones(nSample,1) - ones(nSample,1).*x(:,10)./(c(17)-x(:,4));

y = [volFlow,pch4,pco2,SIN,TS,VS];

end
