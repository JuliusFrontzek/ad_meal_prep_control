function   generate_dfdp_sym
% 1. Definition der symbolischen Variablen (Zustände, Parameter, Stellgr.)
% 2. Erzeugung der rhs (rechte Seite der DGL) mit symbolischen Variablen
% 3. Berechnung der partiellen Ableitungen df/dx aus dieser rhs
% 4. Speichern von 3. in einem ext. m-file als Jakobi-Matrix
%
% Code-Grundlage von Terrance Wilms - SPI - TU Berlin

%% 1. Automatic Generation of the Jacobians for Bsp_biomodell
% [Terrance Wilms]

% Symbolische Definition der Zustände, Stellgrößen und Parameter aus
% Bsp_biomodell
x = sym('x',[8,1]);
syms u
p = sym('th',[4,1]);

% Zustände, nach denen abgeleitet wird:
Var = x;

% Parameter, nach denen abgeleitet wird:
Par = p;

%% 2. Erzeugung der rhs mit symbolic variables
[Diffsys] = BMR4_AB_ess_ode_sym;

%% 3. Berechnung der analytischen partiellen Ableitung df/dX aus Diffsys
JacX = jacobian(Diffsys,Var);

%% 4. creation of the M-file
S1 = vectorize(JacX);
% S2 = S1(10:end-3);
S2 = S1;
% perform some regular expression replacements: 
S3 = regexprep(S2,' ','');
% S4 = regexprep(S3,',',' ');
S5 = regexprep(S3, ';', ';\n');

% open the file (create it if non-existent)...
fid = fopen('jacX.m','w');

% ... write into it...
fprintf(fid,'%%created automatically with SPI_JacGen_Wilms\n');
fprintf(fid,'%%    %s\n\n',datestr(now));
fprintf(fid,'JacX = %s;',S5);

% ... and close it:
fclose(fid);

%% 5. Berechnung der analytischen partiellen Ableitung df/dp aus Diffsys

JacP = jacobian(Diffsys,Par);

%% 6. creation of the M-file (same procedure as for JacX)
S1 = vectorize(JacP);
% S2 = S1(10:end-3);
S2 = S1;
S3 = regexprep(S2,' ','');
% S4 = regexprep(S3,',',' ');
S5 = regexprep(S3, ';', ';\n');

% open the file (create it if non-existent)...
fid = fopen('jacP.m','w');

% ... write into it...
fprintf(fid,'%%created automatically with SPI_JacGen_Wilms\n');
fprintf(fid,'%%    %s\n\n',datestr(now));
fprintf(fid,'JacP = %s;',S5);

% ... and close it:
fclose(fid);
end
