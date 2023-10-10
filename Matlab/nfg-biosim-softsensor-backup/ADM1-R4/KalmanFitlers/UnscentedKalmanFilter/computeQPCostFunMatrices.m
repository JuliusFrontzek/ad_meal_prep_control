%% Version
% (R2022b) Update 5
% Erstelldatum: 09.10.2023
% Autor: Simon Hellmann

function [HMat,fTranspose] = computeQPCostFunMatrices(sigmaXMinus,yMeas,R,PMinus)
% compute matrix H and vector f of QP optimization problem acc. to 
% Kolas et al. (2009), Eq. (33)

% HMat - quadratic matric H
% fTranspose - transposed vector f
% yMeas - measurement vector
% R - measurement noise cov. matrix
% PMinus - a priori state error cov. matrix

q = numel(yMeas);   % # outputs
nStates = length(sigmaXMinus);  % # states

% use efficient least squares for matrix inversion (for numerical stability): 
RInv = R\eye(q);
PMinusInv = PMinus\eye(nStates); 

% linear output matrix:
D = [1,0,0,0,0,0;   % Sch4
     0,1,0,0,0,0;   % SIC
     0,0,0,0,0,1];  % Xbac

% H and f acc. to Kolas: 
HKolas = D'*RInv*D + PMinusInv; 
fTransposeKolas = yMeas*RInv*D + sigmaXMinus'*PMinusInv; 

% H and f in shape required by Matlab quadprog (multiplied by 1/2):
HMat = 2*HKolas; 
HMat = 0.5*(HMat + HMat');  % ensure symmetry 
fTranspose = -2*fTransposeKolas; 

end