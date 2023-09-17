%% Übung 2: Nichtparametrische Ident im geschl. Regelkreis

close all
clear 
clc

%% Vorbereitung
% Eingangssignal erzeugen
% analog zu SPI 1 - UEB 3
prbs_in = 0.5 * idinput(511,'prbs');    % erzeugen eines PRBS-Signals mit 512 Einträgen zw. -0.1 und 0.1
t_prbs = 0.1;                           % solange wird ein Wert des PRBS-Signals gehalten
prbs_in = [zeros(10,1);prbs_in];        % Zufügen von 0 am Anfgang, damit Verhalten bei Ruhe erkennbar

matlab_prbs.time = [];                  % Damit der Eingang in Simulink verwendet werden kann
matlab_prbs.signals.values = prbs_in;

Ts = 0.005;                             % Abtastschrittweite des Reglers und der Messung in Simulink
tend = 50;                              % gesammte Messzeit bzw. Simulationsdauer

%% Aufgabe 2a: Messdaten erzeugen für die Auto- sowie Kreuzvalidierung
% yu = [y,u];	% Daten zur Identifikation in einen Vektor packen
% yw = [y,w]; % Daten zur Identifikation in einen Vektor packen
% uw = [u,w]; % Daten zur Identifikation in einen Vektor packen
% 
% % data_yu = iddata(
% 
% yuk = [yk,uk];	% Daten zur Identifikation in einen Vektor packen
% ywk = [yk,wk]; % Daten zur Identifikation in einen Vektor packen
% uwk = [uk,wk]; % Daten zur Identifikation in einen Vektor packen
systemname = 'Identversuch_stabil';		% Simulink-Modell
sim(systemname)					% Daten erzeugen

yu = [y,u];	% Daten zur Identifikation in einen Vektor packen
yw = [y,w]; % Daten zur Identifikation in einen Vektor packen
uw = [u,w]; % Daten zur Identifikation in einen Vektor packen

yuk = [yk,uk];	% Daten zur Identifikation in einen Vektor packen
ywk = [yk,wk]; % Daten zur Identifikation in einen Vektor packen
uwk = [uk,wk]; % Daten zur Identifikation in einen Vektor packen

% Ident Daten plotten
figure
subplot(2,1,1)
stairs(t,y)
hold on
stairs(t,w,'k')
title('Messdaten zur Identifikation (oben) und Stellgröße (unten)')
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

%% Kreuzkorrelation, Kreuzkovarianz, Totzeit abschätzen
% Kreuzkorrelation mit xcorr
[r,lags] = xcorr(w,y);
figure
plot(lags, r)
xlabel('lag')
title('xcorr')
grid on
box on

% Kreuzkorvarianz
Nmax=1000;
Cxy = kreuzkovar(u,y,Nmax);
figure
plot(-Nmax:Nmax,Cxy)
ylabel('Crosscovariance')
xlabel('tau')
title('Kreuzkovarianz')

data = iddata(y,u);
delay = delayest(data)

%% Aufgabe 2b: Berechnung der Kreuzleistungsdichte und Nichtparametrische Identifikation
window = 2000;  % Fensterbreite, legt auch Anzahl Freq.-Punkte fest
noverlap = round(window/2);  % Anzahl Ueberlappung (50%)

% etfe, spa, welch ausprobieren

G_spawy       = spa(iddata(y,w,Ts), window);    % Auffruf spa mit Fuehrungsgroeße und Ausgang
d_freq_wy     = iddata(G_spawy);                % Auslesen der Daten aus Spa
Gjw_wy        = abs(d_freq_wy.OutputData);      % Berechnen der Kreuzleistungsdichte
freq_wy       = d_freq_wy.Frequency / (2 * pi); % Auslesen der Frequenz f, inkl. Umrechnung in Hz

G_spawu       = spa(iddata(u,w,Ts), window);    % Auffruf spa mit Fuehrungsgroeße und Stellgroeße
d_freq_wu     = iddata(G_spawu);
Gjw_wu        = abs(d_freq_wu.OutputData);
freq_wu       = d_freq_wu.Frequency / (2 * pi);

G_etfe_wu        = etfe(iddata(u,w,Ts), window, floor(numel(u)/2)); % Auffruf etfe mit Fuehrungsgroeße und Ausgang, mit Einschraenkung N/2
d_freq_etfe_wu   = iddata(G_etfe_wu);
Gjw_etfe_wu      = abs(d_freq_etfe_wu.OutputData);
freq_etfe_wu     = d_freq_etfe_wu.Frequency / (2 * pi);

G_etfe_wy        = etfe(iddata(y,w,Ts), window, floor(numel(u)/2));
d_freq_etfe_wy   = iddata(G_etfe_wy);
Gjw_etfe_wy      = abs(d_freq_etfe_wy.OutputData);
freq_etfe_wy     = d_freq_etfe_wy.Frequency / (2 * pi);

% 1. Methode mit Welch
[y_welch,fy_welch]   = pwelch(y,window,noverlap,[],1/Ts);
[w_welch,fw_welch]   = pwelch(w,window,noverlap,[],1/Ts);
[u_welch,fu_welch]   = pwelch(u,window,noverlap,[],1/Ts);
Gwy_pwelch_abs       = sqrt(y_welch./w_welch);
Gwu_pwelch_abs       = sqrt(u_welch./w_welch);

% 2. Methode mit Welch
[wy_pwelch, ~] = cpsd(w,y,window,noverlap,[],1/Ts);
G_wy_pwelch = abs(wy_pwelch./w_welch);

[wu_pwelch, f_pwelch] = cpsd(w,u,window,noverlap,[],1/Ts);
G_wu_pwelch = abs(wu_pwelch./w_welch);

% Plot
% Fuehrungsgroeße und Ausgang
figure
hold on
    plot(freq_wy, 20*log10(Gjw_wy), 'Linewidth', 2, 'Displayname', 'abs(Spa wy)')
    plot(freq_etfe_wy, 20*log10(Gjw_etfe_wy), 'Linewidth', 2, 'Displayname', 'abs(etfe wy)')
    plot(fy_welch, 20*log10(Gwy_pwelch_abs), 'r-', 'Linewidth', 2, 'Displayname', 'pwelch wy')
    plot(f_pwelch, 20*log10(G_wy_pwelch), 'Linewidth', 2, 'Displayname', 'cpsd wy')
    grid on
    box on
legend
ylabel('Amplitude')
xlabel('Hz')

% Fuehrungsgroeße und Stellgroeße
figure
hold on
    plot(freq_wu, 20*log10(Gjw_wu), 'Linewidth', 2, 'Displayname', 'abs(Spa wu)')
    plot(freq_etfe_wu, 20*log10(Gjw_etfe_wu), 'Linewidth', 2, 'Displayname', 'abs(etfe wu)')
    plot(fy_welch, 20*log10(Gwu_pwelch_abs), 'k-', 'Linewidth', 2, 'Displayname', 'pwelch wu')
    plot(f_pwelch, 20*log10(G_wu_pwelch), 'Linewidth', 2, 'Displayname', 'cpsd wu')
    grid on
    box on
legend
ylabel('Amplitude')
xlabel('Hz')

%% Aufgabe 2c: Überführen in ein parametrisches Modell
% mit myVectorFit - UEB 1 - SPI 2
% Auswahl der Methode

freq    = f_pwelch;
Gjw     = G_wy_pwelch./G_wu_pwelch;
% Simon: ich glaube, das ist falsch, denn Gjw ist eine reelle ÜF, da
% bereits G_wy_pwelch und W_wu_pwelch reell waren. Hier wird zu viel
% approxomiert. Gjw muss eigentlich komplex sein!

% Mit VectorFit
fmin    = 20;	% kleinste Frequenz, die übergeben werden soll;
fmax    = 80;	% größte Frequenz
order   = 3;   % Ordnung des gesuchten Modells
% Gjw     = data_ident.OutputData;
% freq    = data_ident.Frequency / (2 * pi);
G_vectorfit = myVectorFit(fmin, fmax,freq,Gjw,order)

% Vergleich mit ursprünglichem System
figure
bode(G_vectorfit)
legend('Gvectorfit')

%% Aufgabe 2d: Vergleich mit dem Simulink Modell
[num, den] = ss2tf(G_vectorfit.A, G_vectorfit.B, G_vectorfit.C, G_vectorfit.D);
