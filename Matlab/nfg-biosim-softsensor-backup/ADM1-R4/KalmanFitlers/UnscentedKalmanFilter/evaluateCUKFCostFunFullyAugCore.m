%% Version
% (R2022b) Update 5
% Erstelldatum: 21.10.2023
% Autor: Simon Hellmann

function J = evaluateCUKFCostFunFullyAugCore(sigmaX,sigmaXMinus, ...
                addNoiseOnSigmapointsY,yMeas,R,PMinus,g)
% compute the scalar value of the nonlinear cost function acc. to 
% Kolas et al. (2009), Eq. (32), in normalized coordinates

% J - scalar value of cost function
% sigmaX - updated sigma point to compute posterior (design variable of optimization)
% sigmaXMinus - propagated sigma points to compute prior
% yMeas - measurement vector
% R - measurement noise cov. matrix
% PMinus - a priori state error cov. matrix

% call output equation:
yNom = g(sigmaX);
y = yNom + addNoiseOnSigmapointsY; 

q = numel(yMeas);               % # outputs
nStates = numel(sigmaX);        % # states

% use efficient least squares for matrix inversion (for numerical stability): 
RInv = R\eye(q);
PMinusInv = PMinus\eye(nStates); 

% compute cost function:
J = (yMeas - y)'*RInv*(yMeas - y) + ...
    (sigmaX - sigmaXMinus)'*PMinusInv*(sigmaX - sigmaXMinus); 

% XY for Terrance: hier bitte Gradienten berechnen!

end