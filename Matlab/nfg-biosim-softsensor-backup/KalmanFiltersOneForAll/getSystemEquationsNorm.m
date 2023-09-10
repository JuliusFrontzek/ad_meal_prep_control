%% Version
% (R2022b) Update 5
% Erstelldatum: 13.07.2023
% Autor: Simon Hellmann

function  [fNorm,gNorm] = getSystemEquationsNorm(flagModel,flagFrac,nStates,sza,szc,szth)
% returns normalized system equations as function handles

% fNorm - normalized ODE function handle
% gNorm - normalized measurement equation function handle
% flagModel -       3: ADM1-R3; 4: ADM1-R4
% flagFrac -        0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% nStates -         # states
% sza, szc, szth -  sizes of arrays aNum, cNum and thNum

xNormS = sym('xNorm', [nStates 1]); % normalized states as col. vector
syms uNorm real;                    % normalized input
xiNormS = sym('xi', [nStates,1]);   % normalized inlet concentrations 
thS = sym('th', szth);              % time-variant parameters (theta)
cS = sym('c', szc);                 % known & constant time-invariant parameters 
aS = sym('a', sza);                 % petersen matrix with stoichiometric constants
TxS = sym('Tx', [nStates,1]);       % normalization matrix for states
syms Tu real                        % normalization variable for input
            

switch flagModel
    case 3
        q = 8;  % 8 measurements for R3 models
        TyS = sym('Ty', [q,1]);             % normalization matrix for outputs
        
        if flagFrac == 0 % R3-norm
            dynamicsNorm = ADM1_R3_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
            outputsNorm = ADM1_R3_norm_mgl_sym(xNormS, cS, TxS, TyS); 
            
        else % R3-frac-norm
            dynamicsNorm = ADM1_R3_frac_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
            outputsNorm = ADM1_R3_frac_norm_mgl_sym(xNormS, cS, TxS, TyS);             
        
        end % if

    case 4
        q = 6;  % 6 measurements for R4 models
        TyS = sym('Ty', [q,1]);             % normalization matrix for outputs
        
        if flagFrac == 0
            dynamicsNorm = ADM1_R4_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
            outputsNorm = ADM1_R4_norm_mgl_sym(xNormS,cS,TxS,TyS);
        else 
            dynamicsNorm = ADM1_R4_frac_norm_ode_sym(xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu); 
            outputsNorm = ADM1_R4_frac_norm_mgl_sym(xNormS,cS,TxS,TyS);

        end % if

end % switch

% turn into numeric function handles: 
fNorm = matlabFunction(dynamicsNorm, 'Vars', {xNormS, uNorm, xiNormS, thS, cS, aS, TxS, Tu}); 
gNorm = matlabFunction(outputsNorm, 'Vars', {xNormS, cS, TxS, TyS}); 
        
end % fun