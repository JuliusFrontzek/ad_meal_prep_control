function   generate_dhdx_sym
% 1. Definition der symbolischen Variablen (Zustände, Parameter, Stellgr.)
% 2. Erzeugung der Messgleichung mit symbolischen Variablen
% 3. Berechnung der partiellen Ableitungen dh/dx aus dieser rhs
% 4. Speichern von 3. in einem ext. m-file als Jakobi-Matrix
%
% Code-Grundlage von Terrance Wilms - SPI - TU Berlin

%% 1. Automatic Generation of the Jacobians for ADM1-R4
% [Terrance Wilms]

% Symbolische Definition der Zustände, Stellgrößen und Parameter aus
% Bsp_biomodell
x = sym('x',[12,1]);
% syms u
% t = sym('t',[4,1]);

% Zustände, nach denen abgeleitet wird:
Var = x;

% % Parameter, nach denen abgeleitet wird:
% Par = t;

%% 2. Erzeugung der rhs mit symbolic variables
[Diffsys] = BMR4_AB_frac_mgl_sym;

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
% % replace "x3" with "x(3)", same as "c15" with "c(15)"
% pat = "x" + digitsPattern(1,2);     % detect 1 or 2 digits 
% extrPat = extract(S5,pat);    % extracted pattern
% S6 = replace(S5,pat,"(" + extrPat + ")");

% open the file (create it if non-existent)...
fid = fopen('jacMglX.m','w');

% ... write into it...
fprintf(fid,'%%created automatically with SPI_JacGen_Wilms\n');
fprintf(fid,'%%    %s\n\n',datestr(now));
fprintf(fid,'JacMglX = %s;',S5);

% ... and close it:
fclose(fid);

end
