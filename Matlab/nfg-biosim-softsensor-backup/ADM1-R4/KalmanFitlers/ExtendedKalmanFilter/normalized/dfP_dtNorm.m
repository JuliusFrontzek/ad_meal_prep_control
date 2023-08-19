%% Version
% (R2022b) Update 5
% Erstelldatum: 5.6.2023
% Autor: Simon Hellmann

function dxPNormdt = dfP_dtNorm(xPNorm,uNorm,xiNorm,params,QNorm,...
                                fNorm,dfdxNorm,TxNum,TuNum)
% compute right-hand side of ODEs of both normalized states and normalized
% state error covariance matrix P (normalized version)

% dxPNormdt - normalized right-hand side of ODE (of states and dynamics of state error
% covariance matrix)
% xPNorm - normalized states and covariance matrix (reshaped as vector), 
% stacked above each other as column vector
% uNorm - normalized feedVolumeFlow
% xiNorm - normalized inlet concentrations (nStates,1)
% params - struct with stoichiometric coefficients a, aggregated constant
% parameters c and time-variant parameters th
% QNorm - normalized power spectral density matrix of process noise (tuning
% matrix)
% fNorm - function handle of normalized ODEs of system equations
% dfdxNorm - function handle of partial derivatives of df/dx, both in 
% normalized coordinates 
% TxNum - normalization matrix of states
% TuNum - normalization matrix of inputs

    nStates = length(xiNorm); 
    
    % extract constant parameters out of struct: 
    th = params.th; 
    c = params.c; 
    a = params.a;
    
    dxPNormdt = zeros(size(xPNorm));    % placeholder
    
    % separate states and covariance matrix:
    xNorm = xPNorm(1:nStates); 
    PNorm = reshape(xPNorm(nStates+1:end),[nStates,nStates]);
    
    % clipping of negative concentrations:
    % xNorm(xNorm < 0) = 0; 
    
    %% ODEs of states:
    dxNormdt = fNorm(xNorm,uNorm,xiNorm,th,c,a,TxNum,TuNum);    
    dxPNormdt(1:nStates) = dxNormdt; 
    
    %% ODEs of state error covariance matrix P:
    % partial derivatives of the right-hand side w.r.t. states, evaluated
    % at current estimate x (which is actually xHat):
    FNorm = dfdxNorm(xNorm,uNorm,th,c,a,TxNum,TuNum);
    %     Q = zeros(nStates);     % XY Rania
    dPNormdt = FNorm*PNorm + PNorm*FNorm' + QNorm;  % dynamics of normalized state error covariance matrix
    % reshape matrix as long column vector and append values dxPdt:
    dxPNormdt(nStates+1:end) = reshape(dPNormdt,[],1);

end