%% Version
% (R2022b) Update 5
% Erstelldatum: 03.10.2023
% last modified: 06.10.2023
% Autor: Simon Hellmann

% function [hessJ] = evaluateHessCUKFCostFunCore(sigmaX,sigmaXMinus, ...
%                 yMeas,R,PMinus,g)

% the Hessian function has to follow the syntax below, therefore the 
% remaining variables are defined as global: 
function [hessJ] = evaluateHessCUKFCostFunCore(sigmaX,lambda)

% hessJ - Hessian of cost function J
% sigmaX - current value of sigma point (design variable)
% lambda - struct with Lagrange multiplier vectors associated with
% nonlinear constraints (see Matlab docu at
% https://www.mathworks.com/help/optim/ug/writing-scalar-objective-functions.html#bu2w6a9-1)

% XY for Simon

global R
global PMinus 

Jac_h_x = [1 0 0 0 0 0
           0 1 0 0 0 0
           0 0 0 0 0 1];

%% Dummies zum testen
% RInv = eye(q);
% PMinusInv = eye(nStates); 

q = 3;
nStates = 6;

RInv = R\eye(q);
PMinusInv = PMinus\eye(nStates); 

%%%%%%%%%%%%%%%%%%%%%
%% calcualtion of the hessian %
%%%%%%%%%%%%%%%%%%%%%

hessJ = 2*PMinusInv' + 2 * (RInv*Jac_h_x)' * Jac_h_x;

end