
function Cxy = kreuzkovar(x,y,Nmax)
% asymetrische Covarianzfunktion
% Cxy = kreuzkovar(sig4,sig5,Nmax);
x=x(:)';
y=y(:)';

mx = mean(x(Nmax+1:end-Nmax));
my = mean(y(Nmax+1:end-Nmax));
% -- Variante 1 mit Skalarprodukt
Cxy=zeros(1,2*Nmax+1);
for tau = 0:Nmax	% Laufindex über verschiedene zeitliche Verschiebungen
	Cxy(1+Nmax-tau) = (y([Nmax+1:end-Nmax]-tau)-my) * (x(Nmax+1:end-Nmax)-mx)' / (length(x) - 2*Nmax-1);
end % for tau
for tau = 1:Nmax	% Laufindex über verschiedene zeitliche Verschiebungen
	Cxy(1+Nmax+tau) = (y([Nmax+1:end-Nmax]+tau)-my) * (x(Nmax+1:end-Nmax)-mx)' / (length(x) - 2*Nmax-1);
end % for tau
end 