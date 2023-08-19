% forces a matrix to be diagonal positive definite with a lower limit
function P_pos = force_P_pos_definite(POld,eps)
    
    % filter out only the diagonal of P: 
    P_pos = diag( diag(POld));

    nStates = length(diag(P_pos)); % number of states corresponding to POld
    
    % iterate through all n diagonal entries:
    for k=1:nStates
        % exclude NAN
        if isnan(P_pos(k,k))
            P_pos(k,k) = eps;
        end
        % make real
        if imag(P_pos(k,k)) ~= 0
            P_pos(k,k) = abs(P_pos(k,k));
        end
        % saturate extreme unplausible values
        if P_pos(k,k) < eps  
            P_pos(k,k) = eps;   % extremely small values 
        elseif P_pos(k,k) > 1/eps
            P_pos(k,k) = 1/eps; % extremely large values
        end
    end
end
