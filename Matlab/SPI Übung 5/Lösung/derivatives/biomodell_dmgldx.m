% Ableitung der Messgleichung nach den Zuständen

function dhdx = biomodell_dmgldx(x)

%% Definition der Zustände

mX = x(1);
mS = x(2);
V = x(3);

%%

dhdx = [
	1/V   0 -mX/V^2
	  0 1/V -mS/V^2
	];
