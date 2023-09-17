function   SPI_JacGen_Wilms
% 1. Definition der symbolischen Variablen (Zustände, Parameter, Stellgr.)
% 2. Erzeugung der rhs (rechte Seite der DGL) mit symbolischen Variablen
% 3. Berechnung der partiellen Ableitungen df/dX aus dieser rhs
% 4. Speichern von 3. in einem ext. m-file als Jakobi-Matrix
% 5. Berechnung der partiellen Ableitungen df/dP aus dieser rhs
% 6. Speichern von 3. in einem ext. m-file als Jakobi-Matrix
%
% Terrance Wilms - SPI

%% 1. Automatic Generation of the Jacobians for Bsp_biomodell
% [Terrance Wilms]

% Symbolische Definition der Zustände, Stellgrößen und Parameter aus
% Bsp_biomodell
syms('mumax','KS','YXS','mX','mS','V','u','q','cSF');

% Definition der Zustände
Var = [mumax KS YXS];% ...

% Definition der Parameter nach denen abgeleitet wird
Par = [mX mS V];

% Definition der Stellgröße
Uvar = [q cSF];

%% 2. Erzeugung der rhs mit symbolic variables
% Einmaliges Aufrufen der Funktion mit dem Modell (Bsp_biomodell):
[Diffsys] = feval(@Bsp_biomodell,[],[],[],[],[],[],[],[],[]);

% Diffsys sollte so aussehen
% Diffsys =
%              (mS*mX*mumax)/(V*(KS + mS/V))
%  cSF*q - (mS*mX*mumax)/(V*YXS*(KS + mS/V))
%                                          q

%% 3. Berechnung der analytischen partiellen Ableitung df/dX aus Diffsys

% alter Weg vor 2019
%
% for i2 = 1:length(Diffsys)
%     for j2 = 1:length(Var)
%         JacX(i2,j2) = diff(Diffsys(i2),Var(j2));
%     end
% end
% 

% Neue Matlab Funktion
JacX = jacobian(Diffsys,Var);

% JacX sollte so aussehen
% JacX =  
% [      (mS*mX)/(V*(KS + mS/V)),    -(mS*mX*mumax)/(V*(KS + mS/V)^2),                                   0]
% [ -(mS*mX)/(V*YXS*(KS + mS/V)), (mS*mX*mumax)/(V*YXS*(KS + mS/V)^2), (mS*mX*mumax)/(V*YXS^2*(KS + mS/V))]
% [                            0,                                   0,                                   0]

%% 4. creation of the M-file
S1 = vectorize(JacX);
% S2 = S1(10:end-3);
S2 = S1;
% perform some regular expression replacements: 
S3 = regexprep(S2,' ','');
% S4 = regexprep(S3,',',' ');
S5 = regexprep(S3, ';', ';\n');

% open the file (create it if non-existent)...
fid = fopen('Bsp_JacX.m','w');

% ... write into it...
fprintf(fid,'%%created automatically with SPI_JacGen_Wilms\n');
fprintf(fid,'%%    %s\n\n',datestr(now));
fprintf(fid,'JacX = %s;',S5);

% ... and close it:
fclose(fid);

%% 5. Berechnung der analytischen partiellen Ableitung df/dp aus Diffsys

% alter Weg vor 2019
%
% for i2 = 1:length(Diffsys)
%     for j2 = 1:length(Par)
%         JacP(i2,j2) = diff(Diffsys(i2),Par(j2));
%     end
% end

% Neue Matlab Funktion
JacP = jacobian(Diffsys,Par);

%JacP sollte so aussehen
% JacP =
%  
% [      (mS*mumax)/(V*(KS + mS/V)),         (mX*mumax)/(V*(KS + mS/V)) - (mS*mX*mumax)/(V^2*(KS + mS/V)^2),         (mS^2*mX*mumax)/(V^3*(KS + mS/V)^2) - (mS*mX*mumax)/(V^2*(KS + mS/V))]
% [ -(mS*mumax)/(V*YXS*(KS + mS/V)), (mS*mX*mumax)/(V^2*YXS*(KS + mS/V)^2) - (mX*mumax)/(V*YXS*(KS + mS/V)), (mS*mX*mumax)/(V^2*YXS*(KS + mS/V)) - (mS^2*mX*mumax)/(V^3*YXS*(KS + mS/V)^2)]
% [                               0,                                                                      0,                                                                             0]
 

%% 6. creation of the M-file (same procedure as for JacX)
S1 = vectorize(JacP);
% S2 = S1(10:end-3);
S2 = S1;
S3 = regexprep(S2,' ','');
S4 = regexprep(S3,',',' ');
S5 = regexprep(S4, ';', '\n');
fid = fopen('Bsp_JacP.m','w');

fprintf(fid,'%%created automatically with SPI_JacGen_Wilms\n');
fprintf(fid,'%%    %s\n\n',datestr(now));

fprintf(fid,'JacP = %s;',S5);
fclose(fid);
end
