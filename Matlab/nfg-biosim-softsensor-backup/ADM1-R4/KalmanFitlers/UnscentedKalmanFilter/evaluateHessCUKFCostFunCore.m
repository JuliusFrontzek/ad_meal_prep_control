%% Version
% (R2022b) Update 5
% Erstelldatum: 03.10.2023
% last modified: 06.10.2023
% Autor: Simon Hellmann

function [hessJ] = evaluateHessCUKFCostFunCore(sigmaX,lambda,R,PMinus)

% hessJ - Hessian of cost function J
% sigmaX - current value of sigma point (design variable)
% lambda - struct with Lagrange multiplier vectors associated with
% nonlinear constraints (see Matlab docu at
% https://www.mathworks.com/help/optim/ug/writing-scalar-objective-functions.html#bu2w6a9-1)
% R - measurement noise cov. matric
% PMinus - a priori state error cov. matrix

Jac_h_x = [1 0 0 0 0 0
           0 1 0 0 0 0
           0 0 0 0 0 1];

q = 3;
nStates = 6;

RInv = R\eye(q);
PMinusInv = PMinus\eye(nStates); 

%%%%%%%%%%%%%%%%%%%%%
%% calcualtion of the hessian %
%%%%%%%%%%%%%%%%%%%%%

hessJ = 2*PMinusInv' + 2 * (RInv*Jac_h_x)' * Jac_h_x;

end