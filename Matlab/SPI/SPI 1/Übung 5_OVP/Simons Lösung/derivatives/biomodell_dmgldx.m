% Ableitung der Messgleichung g nach den Zust�nden x

function dgdx = biomodell_dmgldx(x)

%% Definition der Zust�nde

mX = x(1);
mS = x(2);
V = x(3);

%%

dgdx = [
	1/V   0 -mX/V^2
	  0 1/V -mS/V^2
	];
