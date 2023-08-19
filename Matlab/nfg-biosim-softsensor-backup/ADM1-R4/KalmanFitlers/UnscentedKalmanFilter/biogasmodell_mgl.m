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
rho = pFix(11);

% global counterY

%% compute 3 online measurements:
nStates = length(x); 

% allow only positive concentrations: 
x(x<0) = 0; 

pCh4       = R*T.*x(nStates-1)./16;   % [bar]
pCo2       = R*T.*x(nStates)./44;     % [bar]
pGas       = pCh4 + pCo2 + ph2o;        % [bar]
qGas       = kp.*(pGas - patm).*pGas./patm./24; % [L/d] -> [L/h]

%% compute 3 offline measurements (assumed to be available at the same frequency): 
S_IN = x(3);                              % [g/l]
TS = 1 - 1/rho .* x(4);               % [-]
VS = 1 - 1/(rho - x(4)) .* x(9);   % [-]

%% Messgleichungen
y = zeros(1,6);     % place holder

y(1:3) = [qGas,pCh4,pCo2];    % online measurements
y(4:6) = [S_IN,TS,VS];        % offline measurements (incorrectly assumed available at the same time instances as the online ones)

% counterY = counterY + 1

end
