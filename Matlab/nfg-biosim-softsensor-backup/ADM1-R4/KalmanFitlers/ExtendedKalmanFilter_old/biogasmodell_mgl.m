% Messgleichung für's ADM1-R4-frac

function y = biogasmodell_mgl(x,pFix)
% compute the measurement outputs y of dim [ny,nSample] from: 
% x - matrix of state trajectories of dim [nStates,nSample] -> row vectors
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
[n,nSample] = size(x); 

pCh4       = R*T.*x(n-1,:)./16;
pCo2       = R*T.*x(n,:)./44;
pGas       = pCh4 + pCo2 + ph2o;
qGas       = kp.*(pGas - patm).*pGas./patm./24; % [L/d] --> [L/h]

%% compute 3 offline measurements (assumed to be available at the same frequency): 
% XY: das ist sicher noch nicht richtig
S_IN = x(3,:); 
TS = ones(1,nSample) - ones(1,nSample)/rho .* x(4,:);
VS = ones(1,nSample) - ones(1,nSample)./(rho - x(4,:)) .* x(10,:);

%% Messgleichungen
y = zeros(6,nSample);     % place holder

y(1:3,:) = [qGas;pCh4;pCo2];    % online measurements
y(4:6,:) = [S_IN;TS;VS];        % offline measurements (incorrectly assumed available at the same time instances as the online ones)

end
