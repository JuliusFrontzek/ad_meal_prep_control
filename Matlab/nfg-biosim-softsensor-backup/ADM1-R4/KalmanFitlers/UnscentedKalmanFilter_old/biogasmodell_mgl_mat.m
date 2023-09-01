% Messgleichung

function y = biogasmodell_mgl_mat(x,pFix)
% compute the measurement outputs y of dim [ny,nSample] from: 
% x - matrix of state trajectories of dim [nSample,nStates]
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

%% compute 3 online measurements:
[nSample,nStates] = size(x); 

% allow only positive concentrations: 
x(x<0) = 0; 

pCh4       = R*T.*x(:,nStates-1)./16;   % [bar]
pCo2       = R*T.*x(:,nStates)./44;     % [bar]
pGas       = pCh4 + pCo2 + ph2o;        % [bar]
qGas       = kp.*(pGas - patm).*pGas./patm./24; % [L/d] -> [L/h]

%% compute 3 offline measurements (assumed to be available at the same frequency): 
S_IN = x(:,3);                              % [g/l]
TS = ones(nSample,1) - ones(nSample,1)/rho .* x(:,4);               % [-]
VS = ones(nSample,1) - ones(nSample,1)./(rho - x(:,4)) .* x(:,9);   % [-]

%% Messgleichungen
y = zeros(nSample,6);     % place holder

y(:,1:3) = [qGas,pCh4,pCo2];    % online measurements
y(:,4:6) = [S_IN,TS,VS];        % offline measurements (incorrectly assumed available at the same time instances as the online ones)

end
