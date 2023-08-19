function dx = BMR4_AB_frac_ode_soeren(t,x,physParams,input,modParams)
% delivers the right-hand side of the homogeneous ODE system of the
% BMR4+AB_frac (incl. balancing of ash & 2 CH fractions) with constant 
% values for u and xin using Manuels/Sörens petersen matrix
    
    % expand parameter vectors for better understanding
    % physical parameters first:
    V_liq = physParams(1); 
    V_gas = physParams(2); 
    p_atm = physParams(3); 

    % Model parameters 
    K_H_ch4 =   modParams(1); 
    K_H_co2 =   modParams(2); 
    R =         modParams(3); 
    T =         modParams(4); 
    k_La =      modParams(5); 
    k_chF =     modParams(6); 
    k_chS =     modParams(7); 
    k_dec =     modParams(8); 
    k_li =      modParams(9); 
    k_p =       modParams(10); 
    k_pr =      modParams(11); 
    p_h2o =     modParams(12);  
    fracChFast= modParams(13);

    n = length(x);
    u = input(1);           % volume flow of substrate feed
    xIn = zeros(n,1);       % make sure xIn is always a column vector 
    xIn(:) = input(2:end);  % inlet concentrations

    % subdivide carbohydrates into fast and slow acc. to fraction parameter
    xInChTot = xIn(5); 
    xInChF = fracChFast*xInChTot; 
    xInChS = (1-fracChFast)*xInChTot; 
    xIn([5,6]) = [xInChF,xInChS];
    
    % 4 Algebraic equations
    p_ch4 = R*T * x(n-1)/16;
    p_co2 = R*T * x(n)/44;
    p_gas = p_ch4 + p_co2 + p_h2o;
    q_gas = k_p*(p_gas - p_atm)*p_gas/p_atm;
    
    % liquid and gas volume flow:
    qV = zeros(n,1);    % ensure that qV is always a column vector
    qV(1:n-2,1) = u/V_liq;
    qV(n-1:n,1) = q_gas/V_gas;
    
    % extract states for better readability:
    xCell = num2cell(x([1,2,5:9])); 
    [Sch4,SIC,XchF,XchS,Xpr,Xli,Xbac] = xCell{:};
    
    % 6 Rate equations
    rate = [k_chF * XchF;
            k_chS * XchS;
            k_pr * Xpr;
            k_li * Xli;
            k_dec* Xbac;
            k_La *(Sch4 - 16*K_H_ch4*p_ch4);
            k_La *(SIC - 44*K_H_co2*p_co2)];

    % Petersen-Matrix ADM1-R4 nach Weinrich 2017, S. 211 mit Korrektur bei
    % alpha_8 und unter Berücksichtigung von X_ash (inert)
    %           S_ch4    S_IC      S_IN        S_h20       X_chF   X_chS   X_pr    X_li     X_bac      X_ash   S_ch4,g     S_co2,g  
    petersen =  [0.2482,  0.6809,  -0.0207,    -0.0456,    -1,      0,      0,      0,       0.1372,    0,      0,          0; 
                 0.2482,  0.6809,  -0.0207,    -0.0456,     0,     -1,      0,      0,       0.1372,    0,      0,          0;
                 0.3221,  0.7954,   0.1689,    -0.4588,     0,      0,     -1,      0,       0.1723,    0,      0,          0; 
                 0.6393,  0.5817,  -0.0344,    -0.4152,     0,      0,      0,     -1,       0.2286,    0,      0,          0; 
                 0,       0,        0,          0,          0.18,   0,      0.77,   0.05,   -1,         0,      0,          0;
                -1,       0,        0,          0,          0,      0,      0,      0,       0,         0,      V_liq/V_gas,0; 
                 0,      -1,        0,          0,          0,      0,      0,      0,       0,         0,      0,          V_liq/V_gas]';
             
    production    = petersen * rate;
    convectionIn  = qV .* xIn;
    convectionOut = qV .* (-x);
    dx = convectionIn + convectionOut + production;    
end 
