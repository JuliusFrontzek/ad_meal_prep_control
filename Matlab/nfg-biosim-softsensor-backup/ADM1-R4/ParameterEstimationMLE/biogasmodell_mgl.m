%% Version
% (R2022b) Update 5
% Erstelldatum: 28.02.23
% Autor: Simon Hellmann

% Messgleichung
function y = biogasmodell_mgl(x,pFix)
% compute the online measurement outputs y of dim [ny,nSample] from: 
% x - matrix of state trajectories of dim [nSample,nStates] -> column vectors
% pFix - fixed system paramters

%% extract fixed parameters for better understanding: 
% K_ch4 = pFix(1); 
% K_co2 = pFix(2); 
R       = pFix(3);
T_op    = pFix(4);
% kla = pFix(5);
k_p     = pFix(6);
p_h2o   = pFix(7);
% V_liq = pFix(8);
% V_gas = pFix(9);
p_atm   = pFix(10);
rho_liq = pFix(11);

% compute 3 online measurements:
[nSample,n] = size(x); 

p_ch4   = R*T_op.*x(:,n-1)./16;
p_co2   = R*T_op.*x(:,n)./44;
p_gas   = p_ch4 + p_co2 + p_h2o;
q_qas   = k_p.*(p_gas - p_atm).*p_gas./p_atm;

% compute 3 offline measurements (assumed to be available at the same frequency): 
S_IN    = x(:,3); 
TS      = ones(nSample,1) - ones(nSample,1)/rho_liq .* x(:,4);
VS      = ones(nSample,1) - ones(nSample,1)./(rho_liq - x(:,4)) .* x(:,9);

% summarize in matrix:
y       = [q_qas,p_ch4,p_co2,S_IN,TS,VS];    % online and offline measurements
end
