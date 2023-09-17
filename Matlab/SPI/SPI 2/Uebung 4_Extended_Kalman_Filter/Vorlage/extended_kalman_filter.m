function [] = extended_kalman_filter()

%% Tuning

% Initialisierung


% Parametersatz I
p_KF = [0.1,0.2,0.6];


% Parametersatz II
% p_KF = [0.11,0.205,0.59];


% Tuning für Rauschen auf x0 in Bioprocess




%% Time Update

x_hat0 = x_hat';

[t,x] = ode45(@dfdt_P,t_span,x_hat0,[],u,p_KF);

x_hat = x(end,1:3)';


%% Measurement Update




x_hat = x_hat';
end




