% compute output y of current state x
function yNorm = measurementFcnADM1R4_norm(xNorm,gNorm,c,TxNum,TyNum)

yNorm = gNorm(xNorm,c,TxNum,TyNum); % predicted model output

end