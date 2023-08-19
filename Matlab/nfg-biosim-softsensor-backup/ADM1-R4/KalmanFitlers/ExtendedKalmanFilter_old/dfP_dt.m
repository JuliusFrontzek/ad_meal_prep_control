%% Version
% (R2022b) Update 5
% Erstelldatum: 21.4.2023
% Autor: Simon Hellmann

function dxPdt = dfP_dt(xP,u,xi,params,Q,f,dfdx)
% compute right-hand side of ODE of both states and state error covariance 
% matrix P. 
% dxPdt - right-hand side of ODE (of states and dynamics of state error
% covariance matrix)
% xP - states and covariance matrix (reshaped as vector), stacked above
% each other as column vector
% u - feedVolumeFlow
% xi - inlet concentrations (nStates,1)
% params - struct with stoichiometric coefficients a, aggregated constant
% parameters c and time-variant parameters th
% Q - power spectral density matrix of process noise
% f - function handle of ODEs of system equations
% dfdx - function handle of partial derivatives of df/dx 

    nStates = length(xi); 
    
    % extract constant parameters out of struct: 
    th = params.th; 
    c = params.c; 
    a = params.a;
    
    dxPdt = zeros(size(xP));    % placeholder
    
    % separate states and covariance matrix:
    x = xP(1:nStates); 
    P = reshape(xP(nStates+1:end),[nStates,nStates]);
    
    % clipping of negative concentrations:
    % x(x < 0) = 0; 
    
    %% ODEs of states:
    dxdt = f(x,u,xi,th,c,a);
    dxPdt(1:nStates) = dxdt; 
    
    %% ODEs of state error covariance matrix:
    
    % partial derivatives of the right-hand side w.r.t. states, evaluated
    % at current estimate x (which is actually xHat):
    F = dfdx(x,u,th,c,a);
    dPdt = F*P + P*F' + Q;  % dynamics of state error covariance matrix
    % reshape matrix as long column vector and append values dxPdt:
    dxPdt(nStates+1:end) = reshape(dPdt,[],1);

end