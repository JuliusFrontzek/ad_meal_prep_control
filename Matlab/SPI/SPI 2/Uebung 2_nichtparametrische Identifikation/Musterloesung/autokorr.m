function Rxx = autokorr(x,Nmax)
% asymetrische Autokorrelationsfunktion
x=x(:)';
% %Initialisieren  von Rxx
Rxx=zeros(1,2*Nmax+1);
% -- Variante 1 mit Skalarprodukt
for tau = 0:Nmax	% Laufindex über verschiedene zeitliche Verschiebungen (verschiebung um -tau)
	Rxx(1+Nmax-tau) = x([Nmax+1:end-Nmax]-tau) * x(Nmax+1:end-Nmax)' / (length(x) - 2*Nmax);
end % for tau
for tau = 1:Nmax	% Laufindex über verschiedene zeitliche Verschiebungen (verschiebung um +tau)
	Rxx(1+Nmax+tau) = x([Nmax+1:end-Nmax]+tau) * x(Nmax+1:end-Nmax)' / (length(x) - 2*Nmax);
end % for tau
end % function autokorr
