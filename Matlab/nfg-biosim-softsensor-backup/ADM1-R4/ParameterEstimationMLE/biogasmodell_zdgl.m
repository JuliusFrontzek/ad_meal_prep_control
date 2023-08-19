%% Version
% (R2022b) Update 5
% Erstelldatum: Ende 2022
% Autor: Simon Hellmann

% DGL zum ADM1-R4: 
function dx_dt = biogasmodell_zdgl(t, x, pEst, input, pFix) % vorher:(t,x,u,p)
% liefert die rechte Seite des DGL-Systems bei fester Fütterung u und Parametern p

%% Grenzen der Zustände
% for idxX = 1:length(x)   
%     % Harte Korrektur, damit keine negative Konzentrationen erhalten werden.
% 	if x(idxX) < 0
% 		x(idxX) = 0;
% 	end % if
% end % for idxX

% delivers the right-hand side of the homogeneous ODE system of the
% BMR4+AB (incl. balancing of ash) with constant values for u and xin 
    n       = length(x);
    u       = input(1);     % volume flow of substrate feed
    xIn     = zeros(n,1);   % make sure xIn is always a column vector 
    xIn(:)  = input(2:end); % inlet concentrations

    % expand parameter vectors for better understanding
    % physical parameters first:   
    K_H_ch4 = pFix( 1); 
    K_H_co2 = pFix( 2); 
    R       = pFix( 3);
    T_op    = pFix( 4);
    k_La    = pFix( 5);
    k_p     = pFix( 6);
    p_h2o   = pFix( 7);
    V_liq   = pFix( 8);
    V_gas   = pFix( 9);
    p_atm   = pFix(10);
%     rho_liq = pFix(11); 
    % Then, unknown kinetic model parameters to be ESTimated:
    k_ch    = pEst(1); 
    k_pr    = pEst(2);
    k_li    = pEst(3);
    k_dec   = pEst(4);
    
    % 4 Algebraic equations
    p_ch4   = R*T_op * x(n-1)/16;
    p_co2   = R*T_op * x(n)/44;
    p_gas   = p_ch4 + p_co2 + p_h2o;
    q_gas   = k_p*(p_gas - p_atm)*p_gas/p_atm;
    
    % liquid and gas volume flow:
    qV      = zeros(n,1);    % ensure that qV is always a column vector
    qV(1:n-2,1) = u/V_liq;
    qV(n-1:n,1) = q_gas/V_gas;
    
    % extract states for better readability:
    xCell   = num2cell(x([1,2,5:8])); 
    [S_ch4,S_IC,X_ch,X_pr,X_li,X_bac] = xCell{:};
    
    % 6 Rate equations
    rate = [k_ch * X_ch;
            k_pr * X_pr;
            k_li * X_li;
            k_dec* X_bac;
            k_La *(S_ch4 - 16*K_H_ch4*p_ch4);
            k_La *(S_IC  - 44*K_H_co2*p_co2)];

    % Petersen-Matrix ADM1-R4 nach Weinrich 2021, https://github.com/soerenweinrich/ADM1/blob/main/Documents/Weinrich_2021_Model_structures.pdf
    % unter Berücksichtigung von X_ash (inert)
    %           S_ch4     S_IC      S_IN        S_h2o       X_ch    X_pr    X_li     X_bac      X_ash   S_ch4,g     S_co2,g  
    petersen = [ 0.24819, 0.68087, -0.02065,   -0.045576,  -1,      0,      0,       0.13716,   0,      0,          0; 
                 0.32208, 0.79543,  0.16892,   -0.45876,    0,     -1,      0,       0.17233,   0,      0,          0; 
                 0.63928, 0.58172, -0.034418,  -0.41518,    0,      0,     -1,       0.2286,    0,      0,          0; 
                 0,       0,        0,          0,          0.18,   0.77,   0.05,   -1,         0,      0,          0;
                -1,       0,        0,          0,          0,      0,      0,       0,         0,      V_liq/V_gas,0; 
                 0,      -1,        0,          0,          0,      0,      0,       0,         0,      0,          V_liq/V_gas].';
             
    production    = petersen * rate;
    convectionIn  = qV .* xIn;
    convectionOut = qV .* (-x);
    dx_dt         = convectionIn + convectionOut + production;
end
