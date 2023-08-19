function V = evaluateCEKFCostFun(x,yMeas,R,xMinus,PMinus,c)

% rufe Messgleichung auf:
yModel = BMR4_AB_frac_mgl(x,c); 

nOutputs = length(yMeas);
nStates = length(xMinus); 
% use efficient least squares solution of matrix inverses: 
RInv = R\eye(nOutputs);
PInv = PMinus\eye(nStates); 

V = 0.5*(yMeas - yModel)'*RInv*(yMeas - yModel) + ...
    0.5*(x - xMinus)'*PInv*(x - xMinus); 

end