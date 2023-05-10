% Fisher'sche Informationsmatrix f�r einfaches biologisches Modell

function [ERGEBNIS] = biomodell_lsg_zdgl(MESS,p)

% Anzahl der Zust�nde
n = length(MESS(1).x0);

np=length(p);
% hier noch allgemeine routine schreiben zum finden von y
ny=2; 

for idxMess = 1:length(MESS)
	
	% Anfangsbedingungen f�r DGLn mit Sensitivit�ten
	x0 = [MESS(idxMess).x0];	
    x_sim = [];
	y = biomodell_mgl(x0);
    x=x0;
	for idxT = 2:length(MESS(idxMess).t)
		
%% Simulation des Systems ohne Sensitivit�ten
		
		% Simulation zwischen zwei Messpunkten
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
		
		[T,X] = ode45(@biomodell_zdgl,tspan,x0,[],MESS(idxMess).u,p);
		x0 = X(end,:)';
		x(:,idxT) = x0;
		x_sim = [x_sim [T'; X(:,1:n)']];
		y(:,idxT) = biomodell_mgl(x0);
		
	end % for idxT
	
	% Ausgabe des Simulationsstructs
	ERGEBNIS(idxMess).t = MESS(idxMess).t;
	ERGEBNIS(idxMess).x0 = MESS(idxMess).x0;
	ERGEBNIS(idxMess).x = x;
	ERGEBNIS(idxMess).x_sim = x_sim;
	ERGEBNIS(idxMess).u = MESS(idxMess).u;
	ERGEBNIS(idxMess).y = y;
	
end % for idxMess

