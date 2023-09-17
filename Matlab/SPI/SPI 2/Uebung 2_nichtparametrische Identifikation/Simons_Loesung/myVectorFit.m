function Fit_SS = myVectorFit(fmin,fmax,freq,Gjw,order)
%MYVECTORFIT passt eine Übertragungsfunktion mit VectorFit an den
%Frequenzgang an
%  Fit_SS = myVectorFit(fmin,fmax,freq,Gjw,order)

f = freq;		
G = Gjw; 

% Modell entsprechend der gewählten minimalen und maximalen Frequenz beschneiden

% alte, fehleranfällige Darstellung:
index_max = find(f>fmax);
index_max = index_max(1);
index_min = find(f>fmin);
index_min = index_min(1);

f_red = f(index_min:index_max);
G_red = G(index_min:index_max);

% Simons Alternativ-Darstellung: 
% idxMinMax = f<=fmax & f>=fmin; 
% f_red = f(idxMinMax); 
% G_red = G(idxMinMax); 

figure
subplot(2,1,1)
plot(f,20*log10(abs(G)),f_red,20*log10(abs(G_red)),'r.')
legend('G','für Vectorfit')
xlabel('f [Hz]')
title('Ausgewählte Punkte von G für eine Anpassung mit Vectorfit')
ylabel('Amplitude [dB]')
subplot(2,1,2)
plot(f,180/pi*unwrap(angle(G)),f_red,180/pi*unwrap(angle(G_red)),'r.')
xlabel('f [Hz]')
ylabel('Phase [deg]')

% reduziertes, nichtparametrisches Modell speichern
save G_red G_red f_red

%% Vectorfitting

% Skript nutzt Vectorfitting zur Anpassung von Frequenzgängen.
% Initialisierungspole werden solange durch die neuen Pole ersetzt, bis
% sich keine Verbesserung mehr ergibt. Das Ganze wird für zwei verschiedene
% Gewichte ausprobiert und das beste Ergebnis geliefert.
% Gregor Gelbert, 2006

% Für genaue Erklärungen zur Funktionsweise von Vectorfit gibt es im Ordner
% .\Vectorfit_orginal einen Userguide und ein Paper.

% clear all; % darf nicht auskommentiert werden, führt sonst zu Fehler
clearvars -except order

load G_red		% Laden der beschnittenen Daten aus reduce_f_points.m

% Unwichtig, nur damit die Bilder nicht alle übereinanderliegen
bdwidth		= 15;
topbdwidth	= 30;
fig_width	= 560;
fig_heigth	= 420;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1	= [bdwidth, scnsize(4) - (fig_heigth + topbdwidth)-250, fig_width, fig_heigth];
pos2	= [scnsize(3) - (fig_width + bdwidth), pos1(2), fig_width, fig_heigth];
pos3	= [bdwidth, scnsize(4) / 2 - (fig_heigth + topbdwidth + 25), fig_width, fig_heigth];
pos4	= [scnsize(3) - (fig_width + bdwidth), pos3(2), fig_width, fig_heigth];
figure(1), set(1,'position',pos1);
figure(2), set(2,'position',pos2);

G = G_red.';	% Komplexe Übertragungsfunktion
fp = f_red';

N = order;		% Hier die Ordnung der zu identifizierenden Übertragungsfunktion angeben

Niter1max = 20;                   
Niter2max = 20;         
          
w = 2 * pi * fp;
s = i .* w;
Ns = length(s);

%=====================================
% Rational function approximation of f(s):
%=====================================

% Complex starting poles :
bet = linspace(w(1),w(Ns),N/2);
% bet = linspace(w(1)/2,w(Ns)*2,N/2);
poles = [];
for n = 1:length(bet)
  alf	= -bet(n) * 1e-1;
  poles = [poles (alf - i * bet(n)) (alf + i * bet(n)) ]; 
end

weight = ones(1,Ns); %No weighting
% weight = 1 ./ abs(fp); %Weighting with inverse of magnitude function

VF.relax = 1;      %Use vector fitting with relaxed non-triviality constraint
VF.kill = 2;       %Enforce stable poles
%VF.kill = 0;      %unstable poles are not changed
VF.asymp = 1;      %Include both D, E in fitting    
VF.skip_pole = 0;  %Do not skip pole identification
VF.skip_res = 0;   %Do not skip identification of residues (C,D,E) 
VF.use_normal = 1; %Use Normal Equations
VF.use_sparse = 1; %Use sparse computations
VF.cmplx_ss = 0;   %Create real-only state space model

VF.spy1 = 0;       %No plotting for first stage of vector fitting
VF.spy2 = 1;       %Create magnitude plot for fitting of f(s) 
VF.logx = 0;       %Use linear abscissa axis
VF.logy = 0;       %Use logarithmic ordinate axis 
VF.errplot = 1;    %Include deviation in magnitude plot
VF.phaseplot = 1;  %Include plot of phase angle
VF.legend = 0;     %do NOT include legends in plots

%%%%%%%%%%%%%%%%Step1%%%%%%%%%%%%%
disp('vector fitting step 1...')

dummy = 1;
iter = 1;
nobetfit = 0;
rmserr_best = 0;
format long

while dummy ~= 2
	poles_old = poles;
	[SER,poles,rmserr,fit] = vectfit2(G,s,poles,weight,VF); 
	rmserr
	
	if (rmserr<rmserr_best)||(iter==1)
		rmserr_best=rmserr;    
		poles_best=poles_old;  
		nobetfit=0;
		Best_in=iter;
	else         %  was the privious fit better?
		nobetfit=nobetfit+1;
	end 
	
	if nobetfit >= Niter1max
		dummy=2;
	end
	iter=iter+1;
end

% [SER,poles,rmserr,fit] = vectfit2(G,s,poles_best,weight,VF); %%nur zur Kontrolle
disp('end step 1...')
disp(['rmserr_best: ' num2str(rmserr_best) ' in iteration ' num2str(Best_in)])
disp(' ')

weight_Step1 = weight;
rmserr_Step1 = rmserr_best;
poles_Step1 = poles_best;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Step2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VF.kill = 0;
poles = poles_best;
% poles = starting_poles;
clear dummy, dummy = input('Weiter 1; Ende 2?');
disp('vector fitting step 2...')
% weight = 1./(abs(fp));
weight = ones(1,Ns);

% for iter = 1:Niter
% 	if iter == Niter2
% 		VF.legend = 1;
% 	end %Include legend in final plot
% 	disp(['   Iter ' num2str(iter)]) 
% 	dummy = 1;
iter = 1;
nobetfit = 0;
rmserr_best = 0;

while dummy ~= 2
	poles_old=poles;
	[SER,poles,rmserr,fit]=vectfit2(G,s,poles,weight,VF); 
	rmserr
	
	if (rmserr < rmserr_best) ||(iter==1)
		rmserr_best = rmserr;    
		poles_best = poles_old;  
		nobetfit = 0;
		Best_in = iter;
	else         %  was the privious fit better?
		nobetfit = nobetfit + 1;
	end
	
	if nobetfit >= Niter2max
		dummy = 2;
	end
	iter = iter+1;
end

% [SER,poles,rmserr,fit] = vectfit2(G,s,poles_best,weight,VF); %%nur zur Kontrolle
disp('end step 2...')
disp(['rmserr_best: ' num2str(rmserr_best) ' in iteration ' num2str(Best_in)])
disp(' ')

% if fit in Step 1 was better as in step2
% Bestfit
if rmserr_Step1 < rmserr_best
	weight = weight_Step1;
	poles_best = poles_Step1;
	disp('Fit in Step1 was better ...')
else
	disp('Fit in Step2 was better ...')  
end

display('bester Fit:')

[SER,poles,rmserr,fit] = vectfit2(G,s,poles_best,weight,VF); %%nur zur Kontrolle
rmserr

A = SER.A; B = SER.B; C = SER.C; D = SER.D; E = SER.E; %Renaming variables

A = full(A); B = full(B);

Fit_SS = ss(A,B,C,D);       %Zustandsraummodell der identfizierten Matrizen

H = freqresp(Fit_SS,2*pi*fp); H = (squeeze(H)).';
hold
plot(fp,180/pi*unwrap(angle(G)-angle(H)),'g') %Fehler in der Phase mit einzeichnen
legend('org','vecfit','error')

%Bode Diagramm der identifizierten Übertragungsfunktion
figure
% set(3,'position',pos3);
bode(Fit_SS);
legend('Modell von Vectorfit')

% figure
% set(4,'position',pos4)
% pzmap(Fit_SS)
% legend('Modell von Vectorfit')

end

