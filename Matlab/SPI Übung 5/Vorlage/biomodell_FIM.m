% Fisher'sche Informationsmatrix f�r einfaches biologisches Modell

function [FM,ERGEBNIS] = biomodell_FIM(MESS,p,invC)

% Anzahl der Zust�nde
n = length(MESS(1).x0);

% Anfangsbedingungen der Sensitivit�ten
XP0 = biomodell_dzdgl0dp(MESS(1).x0,p);

% Initialisierung der FIM
FM = zeros(length(p),length(p));

for idxMess = 1:length(MESS)
	
	% Anfangsbedingungen f�r DGLn mit Sensitivit�ten
	x0_ = [MESS(idxMess).x0;XP0(1:length(p)^2)'];	% Spalten der Matrix werden untereinander gepackt
	
	x_ = x0_;
	x = MESS(idxMess).x0;
	x_sim = [];
	y = biomodell_mgl(x0_);
	for idxT = 2:length(MESS(idxMess).t)
		
%% Simulation des Systems mit Sensitivit�ten
		
		% Simulation zwischen zwei Messpunkten
		tspan = [MESS(idxMess).t(idxT-1) MESS(idxMess).t(idxT)];
		
		[T,X] = ode45(@biomodell_zdgl_XP,tspan,x0_,[],MESS(idxMess).u,p);
		x0_ = X(end,:)';
		x_(:,idxT) = x0_;
		x(:,idxT) = x0_(1:n);
		x_sim = [x_sim [T'; X(:,1:n)']];
		y(:,idxT) = biomodell_mgl(x0_);
		
%%Berechnung der Fisher'schen Informationsmatrix
		
		% Berechnung von dhdx �ber Funktion
		dhdx = biomodell_dmgldx(x(:,idxT));
		
		% Sensitivit�ten aus Simulation		
		XP = reshape(x_(n+1:end,idxT),n,length(p));
		
		% dhdp
		dhdp = zeros(size(dhdx,1),size(XP,2));
		
        
        
        
		%% Aufgabe 1b: Berechnung von dhdtheta und FIM
		dhdtheta = ;
		FM = FM + ;
		
        
        
        
	end % for idxT
	
	% Ausgabe des Simulationsstructs
	ERGEBNIS(idxMess).t = MESS(idxMess).t;
	ERGEBNIS(idxMess).x0 = MESS(idxMess).x0;
	ERGEBNIS(idxMess).x = x;
	ERGEBNIS(idxMess).x_sim = x_sim;
	ERGEBNIS(idxMess).u = MESS(idxMess).u;
	ERGEBNIS(idxMess).y = y;
	
end % for idxMess
