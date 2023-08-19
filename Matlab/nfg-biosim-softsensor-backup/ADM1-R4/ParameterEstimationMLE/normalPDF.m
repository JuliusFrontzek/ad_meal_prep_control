%% Version
% (R2022b) Update 2
% Erstelldatum: 24.02.23
% Autor: Simon Hellmann


function pdf = normalPDF(x,mu,sigma)

vorFaktor = 1/(sqrt(2*pi)*sigma); 
z = (x - mu).^2./(2*sigma^2); 

% pdf = vorFaktor * exp(-z); 
pdf = exp(-z);

end