% Startwerte der Sensitivit‰t
% Im Normalfall nur Nullen (auﬂer wenn eine AB unbekannt und als Funktion
% der Parameter angeschrieben ist):
function dx0dp = biomodell_dzdgl0dp(x,p)
dx0dp = zeros(length(x),length(p));     % Nullmatrix mit Dim. (n x m)
