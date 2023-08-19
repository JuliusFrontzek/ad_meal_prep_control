function V = evaluateCUKFCostFun(x,yMeas,R,sigmaXProp,PMinus,c)

% sigmaXProp - specific sigma point (nStates,1)

% rufe Messgleichung auf:
yModel = BMR4_AB_mgl_h2o(x,c); 

nOutputs = length(yMeas);
nStates = length(sigmaXProp); 
% use efficient least squares solution of matrix inverses: 
RInv = R\eye(nOutputs);
PInv = PMinus\eye(nStates); 

V = 0.5*(yMeas - yModel)*RInv*(yMeas - yModel)' + ...
    0.5*(x - sigmaXProp)'*PInv*(x - sigmaXProp); 

end