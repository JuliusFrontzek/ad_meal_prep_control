function [xPlus,PPlus,Kv] = extended_kalman_filter(xOld,u,yMeas,tSpan,POld,Q,R)

% Parametersatz I
p_KF = [0.1,0.2,0.6];

% Parametersatz II
% p_KF = [0.11,0.205,0.59];

% Tuning für Rauschen auf x0 in Bioprocess

%% Time Update
xPOld = [xOld;reshape(POld,[],1)];  % combined vector of x_hat and P

% integriere Sytemverhalten (x und P) im Interval t_span:
[t,xPMinus] = ode45(@dfdt_P,tSpan,xPOld,[],u,p_KF,Q);

% extrahiere Zustände und Kovarianzmatrix kurz vor dem "neuen" Messzeitpunkt:
xMinus = xPMinus(end,1:3)';
PMinus = reshape(xPMinus(end,4:end),[3,3]);

%% Measurement Update
H = dhdx(xMinus);       % partial derivative of output, evaluated at xMinus
S = H*PMinus*H' + R;    % auxiliary matrix
K = PMinus*H'*inv(S);   % Kalman Gain matrix

h = messgleichung(xMinus);  % output
Kv = K*(yMeas - h);   % effective correction of Kalman Gain on state estimate (n,1); 
xPlus = xMinus + Kv;    % updated state estimation

PPlus = PMinus - K*H*PMinus; % updated covariance of estimation error

% hard correction of negative state estimations: 
xPlus(xPlus<0) = 0; 
end




