function [Corr,StdDev,relStdDev,EW,EV,CN]=Parameteranalyse(CV,p)
%[Corr,StdDev,relStddev,EW,EV,CN]=Parameteranalye(CV)
% Korrelationsmatrix, Standardabweichung, relative Standardabweichung,
% Eigenwerte, Eigenvektoren, Konditionszahl

[n,m]=size(CV);
Corr=zeros(n,m);
for i=1:n-1 %rechte obere Diagonale von CV
    for j=i+1:m
        	Corr(i,j)=CV(i,j)/sqrt(CV(i,i)*CV(j,j)); %ausrechnen der Korrelation
    end    
end

%% 1c

Corr=Corr+Corr'+eye(n);

% Eigenwerte EW und Eigenvektoren EV
[EV,EW]=;
EW=;

% Konditionszahl
CN=; % oder Cond

% Parameterstandardabweichung
StdDev=;

% Relative Parameterstandardabweichung
relStdDev=;