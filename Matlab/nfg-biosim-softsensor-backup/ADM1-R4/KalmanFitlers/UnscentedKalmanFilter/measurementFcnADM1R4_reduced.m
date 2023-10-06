% compute output y of current state x but not gasVolFlow
function y = measurementFcnADM1R4_reduced(x,g,c)

fullOutput = g(x,c);    % all model outputs
y = fullOutput([1,2,4:end]);  % exclude pCO2 from measurements

end