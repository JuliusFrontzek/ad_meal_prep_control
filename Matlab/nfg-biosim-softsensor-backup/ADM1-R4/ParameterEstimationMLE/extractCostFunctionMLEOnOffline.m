function I = extractCostFunctionMLEOnOffline(p,MESS,pFix)
% gibt nur das Gütemaß zurück, lässt die Messfehler-Kovarianz-Matrix C weg
    
    [I,~,~] = myGuete_pi_MLEOnOffline(p,MESS,pFix); 
end