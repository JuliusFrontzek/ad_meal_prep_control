function dx = BMR4_AB_ess_ode(t,x,physParams,input,modParams)
% delivers the right-hand side of the homogeneous ODE system of the
% BMR4+AB (incl. balancing of ash) with constant values for u and xin 
    
    n = length(x);
    u = input(1);           % volume flow of substrate feed [L/d]
    xIn = zeros(n,1);       % make sure xIn is always a column vector 
    xIn(:) = input(2:end);     % inlet concentrations

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
    k_ch =      modParams(6); 
    k_dec =     modParams(7); 
    k_li =      modParams(8); 
    k_p =       modParams(9); 
    k_pr =      modParams(10); 
    p_h2o =     modParams(11);  
    
    % 4 Algebraic equations
    p_ch4 = R*T * x(n-1)/16;
    p_co2 = R*T * x(n)/44;
    p_gas = p_ch4 + p_co2 + p_h2o;
    q_gas = k_p*(p_gas - p_atm)*p_gas/p_atm;
    
    % liquid and gas volume flow:
    qV = zeros(n,1);    % ensure that qV is always a column vector
    qV(1:n-2) = u/V_liq;
    qV(n-1:n) = q_gas/V_gas;
    
    % extract states for better readability:
    xCell = num2cell(x([1,2,3:6])); 
    [Sch4,SIC,Xch,Xpr,Xli,Xbac] = xCell{:};
    
    % 6 Rate equations
    rate = [k_ch * Xch;
            k_pr * Xpr;
            k_li * Xli;
            k_dec* Xbac;
            k_La *(Sch4 - 16*K_H_ch4*p_ch4);
            k_La *(SIC - 44*K_H_co2*p_co2)];

    % Petersen-Matrix ADM1-R4-ess nach Weinrich 2017
    %           S_ch4     S_IC      X_ch    X_pr    X_li     X_bac      S_ch4,g     S_co2,g  
    petersen = [ 0.2482,  0.6809,  -1,      0,      0,       0.1372,    0,          0; 
                 0.3221,  0.7954,   0,     -1,      0,       0.1723,    0,          0; 
                 0.6393,  0.5817,   0,      0,     -1,       0.2286,    0,          0; 
                 0,       0,        0.18,   0.77,   0.05,   -1,         0,          0;
                -1,       0,        0,      0,      0,       0,         V_liq/V_gas,0; 
                 0,      -1,        0,      0,      0,       0,         0,          V_liq/V_gas]';
             
    production    = petersen * rate;
    convectionIn  = qV .* xIn;
    convectionOut = qV .* (-x);
    dx = convectionIn + convectionOut + production;    
end 







