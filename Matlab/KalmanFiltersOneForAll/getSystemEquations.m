%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  [f,g] = getSystemEquations(flagModel,flagFrac,nStates,sza,szc,szth)
% returns system equations as function handles

% f - ODE function handle
% g - measurement equation function handle
% flagModel -       3: ADM1-R3; 4: ADM1-R4
% flagFrac -        0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% nStates -         # states
% sza, szc, szth -  sizes of arrays aNum, cNum and thNum

% define symbolic ("S") variables (all vector are defined as column vectors)
xS = sym('x', [nStates 1]);     % states as col. vector
syms uS real                    % input
xiS = sym('xi', [nStates,1]);   % inlet concentrations (assumed known) 
thS = sym('th', szth);  % time-variant parameters (theta)
cS = sym('c', szc);     % known & constant time-invariant parameters 
aS = sym('a', sza);     % petersen matrix with stoichiometric constants

switch flagModel
    case 3
        if flagFrac == 0             
            % compute symbolic objects:
            dynamics = ADM1_R3_ode_sym(xS, uS, xiS, thS, cS, aS);
            outputs = ADM1_R3_mgl_sym(xS,cS); 
            % transform into numeric function handles. Note that the independet
            % variables are explicitely defined. Their order must be followed when 
            % using the function handle!
            f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
            g = matlabFunction(outputs, 'Vars', {xS, cS}); 
        end % if
    case 4
        if flagFrac == 0 % R4 in abs coordinates  
            % compute symbolic objects:
            dynamics = ADM1_R4_ode_sym(xS, uS, xiS, thS, cS, aS);
            outputs = ADM1_R4_mgl_sym(xS,cS); 
            
            % transform into numeric function handle. Note that the independet
            % variables are explicitely defined. Their order must be followed when 
            % using the function handle!
            f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
            g = matlabFunction(outputs, 'Vars', {xS, cS}); 

        else % R4-frac:
            % compute symbolic objects:
            dynamics = ADM1_R4_frac_ode_sym(xS, uS, xiS, thS, cS, aS); 
            outputs = ADM1_R4_frac_mgl_sym(xS,cS); 
            
            % transform into numeric function handle. Note that the independet
            % variables are explicitely defined. Their order must be followed when 
            % using the function handle!
            f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
            g = matlabFunction(outputs, 'Vars', {xS, cS}); 
        end
end % switch
end % fun