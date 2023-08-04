%% Version
% (R2022b) Update 5
% Erstelldatum: 12.07.2023
% Autor: Simon Hellmann

function  [f,g] = getSystemEquations(flagModel,flagFrac,nStates,sza,szc,szth)
% XY: Zweck der Funktion und Argumente beschreiben

switch flagModel
    case 3
        if flagFrac == 0 
            % define symbolic ("S") variables (all vector are defined as column vectors)
            xS = sym('x', [nStates 1]);   % states as col. vector
            syms uS real             % input
            xiS = sym('xi', [nStates,1]); % inlet concentrations (assumed known) 
            thS = sym('th', [szth]);  % time-variant parameters (theta)
            cS = sym('c', [szc]);   % known & constant time-invariant parameters 
            aS = sym('a', [sza]);% petersen matrix with stoichiometric constants
            
            dynamics = ADM1_R3_ode_sym(xS, uS, xiS, thS, cS, aS); % symbolic object
            outputs = ADM1_R3_mgl_sym(xS,cS); 
            
            % transform into numeric function handle. Note that the independet
            % variables are explicitely defined. Their order must be followed when 
            % using the function handle!
            f = matlabFunction(dynamics, 'Vars', {xS, uS, xiS, thS, cS, aS}); 
            g = matlabFunction(outputs, 'Vars', {xS, cS}); 

        end % if
end % switch
end % fun