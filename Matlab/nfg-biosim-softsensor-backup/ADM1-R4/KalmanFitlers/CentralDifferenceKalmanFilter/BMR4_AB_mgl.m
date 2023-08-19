% Messgleichung für einzelne Datenpunkte x mit dim = (n,1)

function y = BMR4_AB_mgl(x,c)
% compute the measurement outputs y of dim [ny,1] from: 
% x - state vector at time k of dim [nStates,1] (row vector)
% c - fixed, aggregated system paramters

% [~,nSample] = size(x); 

% Messgleichungen
y = [(c(9)*x(10)^2 + c(10)*x(10)*x(11) + c(11)*x(11)^2 + c(12)*x(10) + c(13)*x(11) + c(14))/24;     % Konversion von L/d -> L/h
     c(15)*x(10); 
     c(16)*x(11); 
     x(3); 
     1 - x(4)/c(17); 
     1 - x(9)/(c(17)-x(4))];
end
