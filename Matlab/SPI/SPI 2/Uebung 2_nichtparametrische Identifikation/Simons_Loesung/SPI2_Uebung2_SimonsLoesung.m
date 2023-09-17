%% Übung 2: Nichtparametrische Ident im geschl. Regelkreis

close all
clear 
clc

%% Vorbereitung: Eingangssignal aus Simulink-Modell erzeugen 
% (analog zu SPI 1 - UEB 3)

prbs_in = 0.1 * idinput(511,'prbs');    % erzeugen eines PRBS-Signals mit 512 Einträgen zw. -0.1 und 0.1
t_prbs = 0.1;                           % solange wird ein Wert des PRBS-Signals gehalten
prbs_in = [zeros(10,1);prbs_in];        % Zufügen von 0 am Anfgang, damit Verhalten bei Ruhe erkennbar

matlab_prbs.time = [];                  % Damit der Eingang in Simulink verwendet werden kann
matlab_prbs.signals.values = prbs_in;

Ts = 0.005;     % Abtastschrittweite in Simulink (Messung & Regler)
tend = 50;      % gesammte Simulationsdauer

systemname = 'Identversuch_stabil';		% Simulink-Modell
sim(systemname) % Daten erzeugen aus Simulink-file

%% Aufgabe 2a: Messdaten erzeugen für die Auto- sowie Kreuzvalidierung
% Daten zur Identifikation in Vektoren packen:
yu = [y,u];	
yw = [y,w]; 
uw = [u,w]; 

% create data objects:
data_uy = iddata(y,u,Ts);   % (output,input,sampling time)
data_wy = iddata(y,w,Ts);
data_wu = iddata(u,w,Ts); 

% Daten zur Kreuzvalidierung in Vektoren packen:
yuk = [yk,uk]; 
ywk = [yk,wk]; 
uwk = [uk,wk];

% Ident Daten plotten
figure
subplot(2,1,1)
stairs(t,y)
hold on
stairs(t,w,'k')
title('Messdaten zur Identifikation und Autovalidierung (oben) und Stellgröße (unten)')
legend('Ausgang y','Eingang w')
subplot(2,1,2)
stairs(t,u)
legend('Stellgöße')
xlabel('t in s')			
hold off

% Kreuzvalidierungsdaten Plotten
figure
subplot(2,1,1)
stairs(t,yk)
hold on
stairs(t,wk,'k')
title('Messdaten zur Validierung (oben) und Stellgröße (unten)')
legend('Ausgang y','Eingang w')
subplot(2,1,2)
stairs(t,uk)
legend('Stellgöße')
xlabel('t in s')			
hold off

%% Aufgabe 2b: Berechnung der Kreuzleistungsdichte und Nichtparametrische Identifikation

% Abschätzung der Totzeit der Strecke: 
nk = delayest(data_uy); 

% Kreuzkorrelation zwischen w und y: 
[R,lags] = xcorr(w,y); 
figure()
subplot(2,1,1)
plot(lags,R); 
xlabel('\tau'); 
ylabel('R_{wy}')
title('Kreuzkorrelation und -Kovarianz zwischen w und y')

% Kreuzkovarianz zwischen w und y: 
Nmax = 1000; % max. shift 
Cwy = kreuzkovar(w,y,Nmax); 
subplot(2,1,2)
plot(-Nmax:Nmax,Cwy)
xlabel('\tau')
ylabel('C_{wy}')

% etfe, spa, welch-methoden ausprobieren
winSize = 500;          % Fensterbreite, legt auch Anzahl Freq.-Punkte fest
noverlap = round(winSize/2);  % Anzahl Ueberlappungen (50%) für Welch-Methode

% Übertragungsfunktionen von w nach y (geschl. RK):
% Kreisfrequenz-Vektor festlegen: 
N = length(w);  % number of samples
nEvalEtfe = floor(N/2);% number of frequencies at which PSD is evaluated. 
% See Video SPI2, VL-Woche 1, 02 ETFE Details. Liegt an der FFT und zu 
% welchen Frequenzen diese ausgewertet wird
fMax = 1/Ts;        % sampling frequency [Hz]
wmax = 2*pi*fMax;   % max. Abfastfrequenz
wNyq = 0.5*wmax;    % Nyquist frequency [rad/s]
wEval = [1:nEvalEtfe]/nEvalEtfe*wNyq; % angular frequencies, at which TF is evaluated 

% über Korrelogramm (spa)
G_wy_spa = spa(data_wy,winSize,wEval);
G_wu_spa = spa(data_wu,winSize,wEval);
G_uy_spa = spa(data_uy,winSize,wEval);

% über Periodogramm (etfe):
G_wy_etfe = etfe(data_wy,winSize,nEvalEtfe);
G_wu_etfe = etfe(data_wu,winSize,nEvalEtfe);
G_uy_etfe = etfe(data_uy,winSize,nEvalEtfe);

% turn idfrd- into iddata-objects... 
G_wy_spa_iddata = iddata(G_wy_spa);     
G_uy_spa_iddata = iddata(G_uy_spa);
G_wy_etfe_iddata = iddata(G_wy_etfe);     
G_uy_etfe_iddata = iddata(G_uy_etfe);

% ... to extract transfer functions and frequencies: 
% G_wy (closed-loop)
% frequencies in Hz:
f_wy_spa = G_wy_spa_iddata.Frequency./(2*pi);
f_wy_etfe = G_wy_etfe_iddata.Frequency./(2*pi);
% complex transfer functions (MANually extracted): 
G_wy_spa_man = G_wy_spa_iddata.OutputData; 
G_wy_etfe_man = G_wy_etfe_iddata.OutputData; 

% G_uy (open-loop)
% frequencies in Hz:
f_uy_spa = G_uy_spa_iddata.Frequency./(2*pi);
f_uy_etfe = G_uy_etfe_iddata.Frequency./(2*pi);
% complex transfer functions (MANually extracted): 
G_uy_spa_man = G_uy_spa_iddata.OutputData; 
G_uy_etfe_man = G_uy_etfe_iddata.OutputData; 

% plot Bode plots of both TransferFunctions: 
% G_wy
figure()
bode(G_wy_spa,G_wy_etfe)
legend('G_{wy,spa}','G_{wy,etfe}','Location','southwest');

% G_wu
figure()
bode(G_wu_spa,G_wu_etfe)
legend('G_{wu,spa}','G_{wu,etfe}','Location','southwest')

% beachte: die Befehle etfe und spa berücksichtigen automatisch, nur bis
% zur Nyquistfrequenz zu gehen; denn ihnen wurden iddata-Objekte übergeben,
% bei deren Definition wiederum die Abtastzeit Ts mit angegeben wurde

%% weiter 2b) berechne ÜF G_wy_man manuell über (Kreuz-)Leistungsdichten aus Periodogrammen: 
% Leistungsdichte w: 
W = fft(w); 
absW = abs(W); 
PhiW = Ts/N*absW.^2; 

% Kreuzleistungsdichte wy: 
Y = fft(y); 
absY = abs(Y); 
PhiWY = Ts/N*absY.*absW; 
PhiWYcc = Ts/N*Y.*conj(W);    % complex-conjugate

% ÜF von w nach y: 
G_wy_man = PhiWYcc./PhiW; 
G_wy_man_approx = PhiWY./PhiW; 
% kürze die ÜF bis auf die Nyquistfrequenz ein (RELevanter Frequenzbereich), 
% denn danach kommen nur noch komplex-konjugierte Amplituden: 
NRel = floor(N/2);
% relevant frequency domain (till nyquist frequency): 
G_wy_man_rel = G_wy_man(1:NRel);  
G_wy_man_approx_rel = G_wy_man_approx(1:NRel);  

% dazu passender Frequenzvektor (NRel/2 Auswertestellen): 
% wVec = linspace(0,1,NRel)*wNyq; 
wVec = [1:NRel]./NRel*wNyq;

% plotte G_wy_man als Bode-Diagramm (nur Amplitude): 
figure()
semilogx(wVec,20*log10(abs(G_wy_man_rel)),'LineWidth',2);
hold on
semilogx(wVec,20*log10(abs(G_wy_man_approx_rel)));
[~,~,nSamples] = size(G_wy_etfe.ResponseData); 
myG_wy_etfe = reshape(G_wy_etfe.ResponseData,[nSamples,1]); % turn multi-dim array into vector
semilogx(wEval,20*log10(abs(myG_wy_etfe)))
xlabel('angular frequency')
ylabel('Amplitude in dB')
legend('manuell correct', 'manuell approx', 'etfe')
title('Approx. & exakte Berechnung der Kreuzleistungsdichten sind gleich')

%% weiter 2b) Berechne ÜF mit Welch-Verfahren: 
% Methode 1: aus (Auto-)Leistungsdichten (nur bei minimalphasigen Systemen
% sinnvoll, sonst falsche Ergebnisse wegen in |G(jw)| fehlende Information 
% über Phasengang)
window = 2000; 
nOverlap = window/2;    %  50% overlap
fs = 1/Ts;              % sampling frequency
[Puu,fWelch] = pwelch(u,window,nOverlap,[],fs); 
[Pyy,~] = pwelch(y,window,nOverlap,[],fs);
[Pww,~] = pwelch(w,window,nOverlap,[],fs);
G_uy_welch1 = sqrt(Pyy./Puu); 
G_wy_welch1 = sqrt(Pyy./Pww); 

% Methode 2: aus Kreuz- und Autoleistungsdichten (bei nicht-minimalphasigen
% Systemen zu bevorzugen)
[Puy,~] = cpsd(u,y,window,nOverlap,[],fs);
[Pwy,~] = cpsd(w,y,window,nOverlap,[],fs);
G_uy_welch2 = Puy./Puu; % Achtung komplex!
G_wy_welch2 = Pwy./Pww; % Achtung komplex!

% Methode 3: NUR aus Kreuzleistungsdichten (print-Skript, S.186):
[Pwy,~] = cpsd(w,y,window,nOverlap,[],fs);
[Pwu,~] = cpsd(w,u,window,nOverlap,[],fs);
G_S_welch = Pwy./Puy;   % Achtung komplex!

% plotte beide Methoden als Bode-Diagramme: 
figure()
% G_uy (open loop): 
subplot(2,1,1) 
semilogx(fWelch,20*log10(G_uy_welch1));
hold on
semilogx(fWelch,20*log10(abs(G_uy_welch2)));
semilogx(fWelch,20*log10(abs(G_S_welch)));
semilogx(f_uy_spa,20*log10(abs(G_uy_spa_man)),'LineWidth',2);
semilogx(f_uy_etfe,20*log10(abs(G_uy_etfe_man)));
xlabel('frequency')
ylabel('|G_{uy}(jw)| in dB')
legend('Welch über Leistungsdichten', 'Welch über Auto-/Kreuzleistungsdichten', ...
       'Welch NUR über Kreuzleistungsdichten (Skript S.186)', ...
       'spa manuell extrahiert', 'etfe manuell extrahiert', ...
       'Location', 'SouthWest')

% G_wy (closed loop): 
subplot(2,1,2) 
semilogx(fWelch,20*log10(G_wy_welch1));
hold on
semilogx(fWelch,20*log10(abs(G_wy_welch2)));
semilogx(f_wy_spa,20*log10(abs(G_wy_spa_man)),'LineWidth',2);
semilogx(f_wy_etfe,20*log10(abs(G_wy_etfe_man)));
xlabel('frequency')
ylabel('|G_{wy}(jw)| in dB')
legend('Welch über Leistungsdichten', ...
       'Welch über Auto-/Kreuzleistungsdichten', ...
       'spa manuell extrahiert', 'etfe manuell extrahiert')

%% Aufgabe 2c: Überführen in ein parametrisches Modell
% mit myVectorFit - UEB 1 - SPI 2
% für die biasfreie Schätzung der Strecken-ÜF G_S = G_uy siehe G_S_welch!

fMin = 20; 
fMax = 80; 
order = 3;  % XY: hiermit herumspielen, bis Ergebnis zufriedenstellend

stateSpaceStrecke = myVectorFit(fMin, fMax, fWelch, G_S_welch, order); 

%% nach MuLö: 
% hier wird glaube ich auch kein richtiges Ergebnis erzielt, denn G wird
% bereits als rein reelle Größe übergeben und nicht als komplexe ÜF!

% [w_welch,f_pwelch]   = pwelch(w,window,noverlap,[],1/Ts);
% [wy_pwelch, ~] = cpsd(w,y,window,noverlap,[],1/Ts);
% [wu_pwelch, f_pwelch] = cpsd(w,u,window,noverlap,[],1/Ts);
% G_wy_pwelch = abs(wy_pwelch./w_welch);
% G_wu_pwelch = abs(wu_pwelch./w_welch);
% 
% freq    = f_pwelch;
% Gjw     = G_wy_pwelch./G_wu_pwelch;
% stateSpaceStreckeMuLo = myVectorFit(fMin, fMax, f_pwelch, Gjw, order); 

%% Aufgabe 2d: Vergleich mit dem Simulink Modell
A = stateSpaceStrecke.A;
B = stateSpaceStrecke.B;
C = stateSpaceStrecke.C;
D = stateSpaceStrecke.D;
[num, den] = ss2tf(A,B,C,D); 

% -> liefert keine vernünftige Anpassungsgüte, wenn man (wie im Skript
% angegeben) die komplexen Vectoren der Kreuzleistungsdichten übergibt. In
% der Musterlösung werden jeweils nur die Beträge der Kreuzleistungsdichten
% übergeben, was natürlich wieder einen reellen Vektor übergibt. Die
% Anpassungsgüte ist aber auch hier nicht ideal. 

%% zusätzliche Aufgaben

% 2e) ist die Kreuzleistungsdichte kommutativ? 
% untersuche das anhand der Kreuzleistungsdichte Phi_uy:

Phi_uy = cpsd(u,y,window,nOverlap,[],fs);
Phi_yu = cpsd(y,u,window,nOverlap,[],fs);
