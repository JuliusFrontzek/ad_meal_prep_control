% Messgleichung

function y = biogasmodell_mgl_off(x,pFix)
% compute the offline measurement outputs y of dim [ny,nSample] from: 
% x - matrix of state trajectories of dim [nStates,nSample] -> row vectors
% pFix - fixed system paramters

%% extract fixed parameters for better understanding: 
% Kch4 = pFix(1); 
% Kco2 = pFix(2); 
% R = pFix(3);
% T = pFix(4);
% kla = pFix(5);
% kp = pFix(6);
% ph2o = pFix(7);
% Vl = pFix(8);
% Vg = pFix(9);
% patm = pFix(10);
rho = pFix(11);

%% compute 3 online measurements:
[~,nSample] = size(x); 

% pCh4       = R*T.*x(n-1,:)./16;
% pCo2       = R*T.*x(n,:)./44;
% pGas       = pCh4 + pCo2 + ph2o;
% qGas       = kp.*(pGas - patm).*pGas./patm;

%% compute 3 offline measurements (assumed to be available at the same frequency): 
S_IN = x(3,:); 
TS = ones(1,nSample) - ones(1,nSample)/rho .* x(4,:);
VS = ones(1,nSample) - ones(1,nSample)./(rho - x(4,:)) .* x(9,:);

%% Messgleichungen
y = zeros(3,nSample);     % place holder

% y(1:3,:) = [pCh4;pCo2;qGas];    % online measurements
y(1:3,:) = [S_IN;TS;VS];        % offline measurements (incorrectly assumed available at the same time instances as the online ones)

end
