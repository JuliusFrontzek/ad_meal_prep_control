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
Corr=Corr+Corr'+eye(n);

% Eigenwerte EW und Eigenvektoren EV
[EV,EW]=eig(CV);
EW=diag(EW);

% Konditionszahl
CN=max(EW)/min(EW); % oder Cond

% Parameterstandardabweichung
StdDev=sqrt(diag(CV))';

% Relative Parameterstandardabweichung
relStdDev=(StdDev./p');