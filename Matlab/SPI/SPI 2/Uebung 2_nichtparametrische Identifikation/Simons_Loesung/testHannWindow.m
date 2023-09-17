%% wie sieht das Hann-Window im Zeitbereich aus? 

clear
close all
clc

Ts = 0.1;   % sampling time [s]

% different window sizes: 
M1 = 20; 
M2 = 40; 
M3 = 80; 

M = [M1,M2,M3]; 

k1 = -(M1-1):M1-1; 
k2 = -(M2-1):M2-1; 
k3 = -(M3-1):M3-1;
w1 = 0.5*(1 - cos(2*pi*k1/(M1 - 1)));
w2 = 0.5*(1 - cos(2*pi*k2/(M2 - 1)));
w3 = 0.5*(1 - cos(2*pi*k3/(M3 - 1)));

% plot results: 
figure()
plot(k1,w1)
hold on
plot(k2,w2,k3,w3)
xlabel('samples')
ylabel('amplitude')
title('Hann Window (time domain)')
legend(['M1=',num2str(M1)], ['M2=',num2str(M2)], ['M3=',num2str(M3)])

%% compute fourier transforms of them
% for the Fourier Transforms, we only want positive k's:
k1 = 0:M1-1; 
k2 = 0:M2-1; 
k3 = 0:M3-1;
w1 = 0.5*(1 - cos(2*pi*k1/(M1 - 1)));
w2 = 0.5*(1 - cos(2*pi*k2/(M2 - 1)));
w3 = 0.5*(1 - cos(2*pi*k3/(M3 - 1)));

% plot with only positive sample sizes: 
figure()
plot(k1,w1)
hold on
plot(k2,w2,k3,w3)
xlabel('samples')
ylabel('amplitude')
title('Hann Window (time domain)')
legend(['M1=',num2str(M1)], ['M2=',num2str(M2)], ['M3=',num2str(M3)])

N1 = length(w1);
N2 = length(w2); 
N3 = length(w3);

% make even numbers out of these:
N1 = 2*floor(N1/2);   
N2 = 2*floor(N2/2);   
N3 = 2*floor(N3/2);   

% the FFT delivers a complex vector. The abs. values spread half-and-half
% on the pair of complex-conjugate frequencies. As we only want to plot real
% frequencies, we take the double of all the vector except the first entry
% (this one is only real and represents the constant component of the signal)
W1 = fft(w1)/N1; 
W2 = fft(w2)/N2;
W3 = fft(w3)/N3;

W1(2:end) = 2*W1(2:end); 
W2(2:end) = 2*W2(2:end); 
W3(2:end) = 2*W3(2:end); 

% set up frequency vector: 
fMax = 1/Ts; 
% full frequency vector
fVec1 = [0:N1-1]/N1*fMax;
fVec2 = [0:N2-1]/N2*fMax;
fVec3 = [0:N3-1]/N3*fMax;
% kürze Frequenzvektor bis zur Nyquistfrequenz
fVec1Rel = fVec1(1:N1/2);
fVec2Rel = fVec2(1:N2/2);
fVec3Rel = fVec3(1:N3/2);
% kürze FFT auf gleiche Länge wie Frequenzvektor und nimm Betrag:
W1Rel = abs(W1(1:N1/2));
W2Rel = abs(W2(1:N2/2));
W3Rel = abs(W3(1:N3/2));

% plot Fourier transforms: 
figure()
subplot(3,1,1)
stem(fVec1Rel,W1Rel);
title('Fourier Transforms of Hann Windows (diff. window sizes)')
ylabel('|W_N(f)|')
legend('window size = 20')
subplot(3,1,2)
stem(fVec2Rel,W2Rel);
ylabel('|W_N(f)|')
legend('window size = 40')
subplot(3,1,3)
stem(fVec3Rel,W3Rel);
xlabel('Frequency [Hz]')
ylabel('|W_N(f)|')
legend('window size = 80')
