%% Version
% (R2022b) Update 5
% Erstelldatum: 13.07.2023
% Autor: Simon Hellmann

function  [fNorm,gNorm] = getSystemEquationsNorm(flagModel,flagFrac,nStates,sza,szc,szth)
% XY: Zweck der Funktion und Argumente beschreiben

% XY: varargout in output-variablen aufnehmen, siehe hier: 
% https://www.mathworks.com/help/matlab/matlab_prog/support-variable-number-of-outputs.html

switch flagModel
    case 3
        q = 8;  % 8 measurements for R3 models
        if flagFrac == 0 
            xNormS = sym('xNorm', [nStates 1]); % normalized states as col. vector
            syms uNorm real;                    % normalized input
            xiNormS = sym('xi', [nStates,1]);   % normalized inlet concentrations 
            thS = sym('th', [szth]);  % time-variant parameters (theta)
            cS = sym('c', [szc]);   % known & constant time-invariant parameters 
            aS = sym('a', [sza]);% petersen matrix with stoichiometric constants
            TxS = sym('Tx', [nStates,1]);       % normalization matrix for states
            TyS = sym('Ty', [q,1]);             % normalization matrix for outputs
            syms Tu real                        % normalization variable for input
            
            dynamicsNorm = ADM1_R3_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
            outputsNorm = ADM1_R3_norm_mgl_sym(xNormS, cS, TxS, TyS); 
            
            % turn into numeric function handles: 
            fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
            gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 

        end % if
end % switch
end % fun