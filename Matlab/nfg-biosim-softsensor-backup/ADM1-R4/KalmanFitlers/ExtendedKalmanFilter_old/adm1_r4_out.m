% Take differential states -> calculate algebraic states and measurements
function y = adm1_r4_out(x, physParms, modParms)
% x expected as matrix of dim [nSamples,nStates] --> time as col vector

    % 11 Process Parameters
%     K_H_ch4     = modParms(1);
%     K_H_co2     = modParms(2);
    R           = modParms(3);
    T_op        = modParms(4); 
%     k_La        = modParms(5);
%     k_ch        = modParms(6); 
%     k_dec       = modParms(7);
%     k_li        = modParms(8);
    k_p         = modParms(9);
%     k_pr        = modParms(10);
    p_h2o       = modParms(11);
    
    % physical parameters
%     V_liq       = physParms(1);
%     V_gas       = physParms(2);
    p_atm       = physParms(3);
    
    % 4 Algebraic equations
    [~,n] = size(x);
    p_ch4       = R*T_op.*x(:,n-1)./16;
    p_co2       = R*T_op.*x(:,n)./44;
    p_gas       = p_ch4 + p_co2 + p_h2o;
    q_gas       = k_p.*(p_gas - p_atm).*p_gas./p_atm;
    
    % Differential States
    y           = x;
    % Algebraic States
    y(:, n+1:n+4) = [p_ch4, p_co2, p_gas, q_gas];
end
