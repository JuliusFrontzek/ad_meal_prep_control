function [x_hat,P0,Kv] = extended_kalman_filter(x_hat,y,t_span,u,P0)

%% Tuning

% Initialisierung
Q = eye(3);
R = eye(2);

% Parametersatz I
p_KF = [0.1,0.2,0.6];

Q = 1e-2*eye(3); 
Q(2,2) = 1e-2;
Q(3,3) = 1e-5;
Q = 0.01*diag([10,75,5]);
R = [1.15^2 0; 0 0.5^2];

% Parametersatz II
% p_KF = [0.11,0.205,0.59];
% 
% Q = eye(3); 
% Q(1,1) = 20; 
% Q(2,2) = 20;
% Q(3,3) = 0.005;
% Q(2,1) = 1;
% Q(1,2) = Q(2,1);

% R = [1.15^2 0; 0 0.5^2];

% Tuning für Rauschen auf x0 in Bioprocess

% Q = 1e1*eye(3);
% Q(2,1) = 10;
% Q(1,2) = Q(2,1);
% Q(3,3) = 10;
% 
% R = [1.15^2 0; 0 0.5^2];



%% Time Update
x_hat
x_hat0 = x_hat';
x_hat0(4:12,1) = P0(1:3^2)';

[t,x] = ode45(@dfdt_P,t_span,x_hat0,[],u,Q,p_KF);

x_hat = x(end,1:3)';
P0 = reshape(x(end,4:end),[3,3]);


%% Measurement Update

H = dhdx(x_hat);

y_hat = messgleichung(x_hat);

% Innovation und Residualkovarianz
v = y-y_hat;
S = H*P0*H'+ R;

% Kalman-Faktor
K = P0*H'*inv(S);

Kv = K*v;

% x_hat(k+1|k+1) und P(k+1|k+1)
x_hat = x_hat+K*v;

P0 = P0-K*H*P0
% P0 = P0-K*H*P0-P0*H'*K'+K*S*K'

x_hat(x_hat<0) = 0;

x_hat = x_hat';
end




