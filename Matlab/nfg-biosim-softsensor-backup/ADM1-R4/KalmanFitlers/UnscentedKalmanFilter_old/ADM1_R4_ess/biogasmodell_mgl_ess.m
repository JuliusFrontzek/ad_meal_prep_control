% Messgleichung

function y = biogasmodell_mgl(x,pFix)
% compute the measurement outputs y of dim ny from: 
% x - vector of state trajectories of dim nStates
% pFix - fixed system paramters

%% extract fixed parameters for better understanding: 
% Kch4 = pFix(1); 
% Kco2 = pFix(2); 
R = pFix(3);
T = pFix(4);
% kla = pFix(5);
kp = pFix(6);
ph2o = pFix(7);
% Vl = pFix(8);
% Vg = pFix(9);
patm = pFix(10);
% rho = pFix(11);

% global counterY

%% compute 3 online measurements:
nStates = length(x); 

% allow only positive concentrations: 
x(x<0) = 0; 

pCh4       = R*T.*x(nStates-1)./16;   % [bar]
pCo2       = R*T.*x(nStates)./44;     % [bar]
pGas       = pCh4 + pCo2 + ph2o;        % [bar]
qGas       = kp.*(pGas - patm).*pGas./patm./3600; % [L/d] -> [L/s]

%% Messgleichungen
y = zeros(1,3);     % place holder

y(1:3) = [qGas,pCh4,pCo2];    % online measurements

% counterY = counterY + 1

end
