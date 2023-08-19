function sumRightHandSide = callFunctionHandle(funHandle,x,u,a)

rightHandSide = funHandle(x,u,a); 
sumRightHandSide = sum(rightHandSide)