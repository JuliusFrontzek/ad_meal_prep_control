% Messgleichung für einzelne Datenpunkte x mit dim = (n,1)

function y = BMR4_AB_frac_mgl(x,c)
% compute the measurement outputs y of dim [ny,nSample] from: 
% x - state vector of dim (nStates,1) -> row vector
% c - fixed, aggregated system paramters

% [nSample,~] = size(x); 

% Messgleichungen
y = [(c(6)*x(11)^2 + c(7)*x(11)*x(12) + c(8)*x(12)^2 + c(9)*x(11) + c(10)*x(12) + c(11))/24;    % [L/h]
     c(12)*x(11);
     c(13)*x(12);   
     x(3);
     1 - x(4)/c(14);
     1 - x(10)/(c(14)-x(4))];  % water and ash concentrations in [kg/l]

end
