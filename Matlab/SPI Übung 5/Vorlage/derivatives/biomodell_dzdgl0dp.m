% Startwerte der Sensitivit�t

function dx0dp = biomodell_dzdgl0dp(x,p)

dx0dp = zeros(length(x),length(p));
