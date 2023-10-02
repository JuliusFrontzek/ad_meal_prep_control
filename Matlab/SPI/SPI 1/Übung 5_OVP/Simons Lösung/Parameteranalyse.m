function [Corr,StdDev,relStdDev,lambda,EV,CN] = Parameteranalyse(CV,p)
%[Corr,StdDev,relStddev,lambda,EV,CN] = Parameteranalye(CV)
% Korrelationsmatrix, Standardabweichung, relative Standardabweichung,
% Eigenwerte, Eigenvektoren, Konditionszahl
% CV - inverse Fishermatrix (untere Grenze für Kov.Matrix des Schätzfehlers --> Cramer-Rao Bound)
% p - Schätzung der Parameter

m = length(CV);  
n = m;  % Achtung: hier Zufall!
CorrPre = zeros(m);
% beachte: die Korrelationsmatrix ist symmetrisch und enthält auf der
% Hauptdiagonalen nur Einsen. Daher genügt es, die rechte obere
% Dreiecksmatrix ab der ersten Nebendiagonalen zu berechnen (CorrPre). Die 
% gesamte Korrelationsmatrix wird dann erhalten, indem die Transponierte 
% von CorrPre zu CorrPre dazuaddiert und schließlich die Hauptdiagonale mit 
% lauter Einsen dazuaddiert wird. Bilde zunächst CorrPre:
for i=1:m-1 % rechte obere Diagonale von CV
    for j=i+1:m
        	CorrPre(i,j)=CV(i,j)/sqrt(CV(i,i)*CV(j,j)); % ausrechnen der Korrelation
    end    
end

%% 1c
% nun wird Mat1 mit ihrer Transponierten und einer Diagonalmatrix aus
% Einsen (Einheitsmatrix) addiert:
Corr = CorrPre + CorrPre' + eye(n);

% Eigenwerte lambda und Eigenvektoren EV von CV (der Schätzfehler-Kov.-Matrix):
[EV,lambda] = eig(CV);

% Konditionszahl (Maß für Steifigkeit)  
CN = cond(CV); % Verhältnis aus größtem und kleinstem Singulärwert 

% Parameterstandardabweichung
StdDev = sqrt(diag(CV));    % Wurzel der (Auto-)Varianzen, die auf der Diagonalen der Kov.-Matrix der Schätzfehler stehen

% Relative Parameterstandardabweichung (Std.-Abw./Parameter-Schätzwert):
relStdDev = StdDev./p;