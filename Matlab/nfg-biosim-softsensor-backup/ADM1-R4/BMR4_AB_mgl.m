% Messgleichung für einzelne Datenpunkte x mit dim = (1,n)

function y = BMR4_AB_mgl(x,c)
% compute the measurement outputs y of dim [1,ny] from: 
% x - state vector at time k of dim [1,nStates] (row vector)
% c - fixed, aggregated system paramters

    % [nSample,~] = size(x); 
    
    % enforce positive concentrations:
    x (x<0) = 0; 

    % Messgleichungen
    volFlow = (c(6)*x(10)^2 + c(7)*x(10)*x(11) + c(8)*x(11)^2 + c(9)*x(10) + c(10)*x(11) + c(11))/24;   % [L/d] -> [L/h]
    pch4 = c(12)*x(10);             % [bar]
    pco2 = c(13)*x(11);             % [bar]
    SIN = x(3);                     % [g/L]
    TS = 1 - x(4)/c(14);            % [-]
    VS = 1 - x(9)/(c(14)-x(4));     % [-]
    
    y = [volFlow,pch4,pco2,SIN,TS,VS];

end
