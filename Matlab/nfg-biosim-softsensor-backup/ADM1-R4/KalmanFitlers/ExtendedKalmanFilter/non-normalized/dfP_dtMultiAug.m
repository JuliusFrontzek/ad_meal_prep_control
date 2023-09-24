%% Version
% (R2022b) Update 5
% Erstelldatum: 24.09.2023
% Autor: Simon Hellmann

function dxPdt = dfP_dtMultiAug(xP,u,xi,params,Q,f,dfdx,nStates,nAug)
% compute right-hand side of ODE of both states and state error covariance 
% matrix P. Allow for multiple augmentations

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
% nStates - # states (without sample state augmentation)
% nAug - # state augmentations
  
    % extract constant parameters out of struct: 
    th = params.th; 
    c = params.c; 
    a = params.a;
    
    dxPdt = zeros(size(xP));    % placeholder
    
    % XY hier klar überprüfen für multiple augmentation!
    % separate states and covariance matrix:
    x = xP(1:(1+nAug)*nStates); 
    P = reshape(xP((1+nAug)*nStates+1:end),[(1+nAug)*nStates,(1+nAug)*nStates]);
    
%     if flagAugmented == 1
%         xAug = xP(1:2*nStates); 
%         PAug = reshape(xP(2*nStates+1:end),[2*nStates,2*nStates]);
%         x = xAug(1:nStates); % ignore the sample state for further computations
%         P = PAug; 
%     else
%         x = xP(1:nStates); 
%         P = reshape(xP(nStates+1:end),[nStates,nStates]);
%     end    

    % clipping of negative concentrations:
    % x(x < 0) = 0; 
    
    %% ODEs of states:
    dxdt = f(x,u,xi,th,c,a);
    
    % XY hier klar überprüfen für multiple augmentation!
    dxsdt = zeros(nStates,1); % maintain zeros for sample state(s)
    dxAugdt = [dxdt;repmat(dxsdt,nAug,1)];  % propagated augmented state
%     dxPdt(1:2*nStates) = dxAugdt; 
    dxPdt(1:(1+nAug)*nStates) = dxAugdt;    

%     if flagAugmented == 1
%         dxsdt = zeros(nStates,1); % maintain sample state
%         dxAugdt = [dxdt;dxsdt]; 
%         dxPdt(1:2*nStates) = dxAugdt; 
%     else 
%         dxPdt(1:nStates) = dxdt; 
%     end
    
    %% ODEs of state error covariance matrix:
    
    % partial derivatives of the right-hand side w.r.t. states, evaluated
    % at current estimate x (which is actually xHat):
    F = dfdx(x,u,th,c,a);
    
    Fs = zeros(nAug*nStates);    % for sample-states
    Fa = blkdiag(F,Fs);     % augmented F-matrix
%         Q = zeros(nStates);   % XY Rania
    Qs = zeros(nAug*nStates);    % for sample-states
    Qa = blkdiag(Q,Qs);     % augmented Q-matrix

    dPAugdt = Fa*P + P*Fa' + Qa;  % dynamics of state error covariance matrix
    % reshape matrix as long column vector and append values dxPdt:
    dxPdt((1+nAug)*nStates+1:end) = reshape(dPAugdt,[],1);
    
    % XY hier klar überprüfen für multiple augmentation!
%     if flagAugmented == 1
%         % augment F- and Q-matrix:
%         Fs = zeros(nStates);    % for sample-states
%         Fa = blkdiag(F,Fs);     % augmented F-matrix
% %         Q = zeros(nStates);   % XY Rania
%         Qs = zeros(nStates);    % for sample-states
%         Qa = blkdiag(Q,Qs);     % augmented Q-matrix
%         
%         dPAugdt = Fa*PAug + PAug*Fa' + Qa;  % dynamics of state error covariance matrix
%         % reshape matrix as long column vector and append values dxPdt:
%         dxPdt(2*nStates+1:end) = reshape(dPAugdt,[],1);
%     else 
% %         Q = zeros(nStates);   % XY Rania
%         dPdt = F*P + P*F' + Q;  % dynamics of state error covariance matrix
%         % reshape matrix as long column vector and append values dxPdt:
%         dxPdt(nStates+1:end) = reshape(dPdt,[],1);
%     end
    
end