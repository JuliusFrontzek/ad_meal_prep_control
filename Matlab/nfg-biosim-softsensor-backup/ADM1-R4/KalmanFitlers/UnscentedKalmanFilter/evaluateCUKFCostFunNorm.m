%% Version
% (R2022b) Update 5
% Erstelldatum: 03.10.2023
% last modified: 03.10.2023
% Autor: Simon Hellmann

function J = evaluateCUKFCostFunNorm(sigmaXNorm,sigmaXMinusNorm, ...
                yMeasNorm,RNorm,PMinusNorm,gNorm,c,TxNum,TyNum)
% compute the scalar value of the nonlinear cost function acc. to 
% Kolas et al. (2009), Eq. (32), in normalized coordinates

% J - scalar value of cost function
% sigmaXNorm - updated (normalized) sigma point to compute posterior (design variable of optimization)
% sigmaXMinusNorm - propagated (normalized) sigma points to compute prior
% yMeasNorm - normalized measurement vector
% RNorm - normalized measurement noise cov. matrix
% PMinusNorm - normalized a priori state error cov. matrix

% call output equation:
yNorm = gNorm(sigmaXNorm,c,TxNum,TyNum);

q = numel(yMeasNorm);               % # outputs
nStates = length(sigmaXMinusNorm);  % # states

% use efficient least squares for matrix inversion (for numerical stability): 
RNormInv = RNorm\eye(q);
PMinusNormInv = PMinusNorm\eye(nStates); 

% compute cost function:
J = (yMeasNorm - yNorm)'*RNormInv*(yMeasNorm - yNorm) + ...
    (sigmaXNorm - sigmaXMinusNorm)'*PMinusNormInv*(sigmaXNorm - sigmaXMinusNorm); 

end