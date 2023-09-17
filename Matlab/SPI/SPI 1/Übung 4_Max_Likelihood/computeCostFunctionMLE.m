function I = computeCostFunctionMLE(p,MESS)
% gibt nur das Gütemaß zurück und lässt die Messfehler-Kovarianz-Matrix C
% weg
    [I,~] = guete_pi_MLE(p,MESS); 
end