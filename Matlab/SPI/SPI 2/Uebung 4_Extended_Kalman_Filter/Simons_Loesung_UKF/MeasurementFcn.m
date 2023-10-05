function y = MeasurementFcn(x,testParam)
% testParam wurde nur hinzugef�gt, um zu testen, ob die Messgleichung auch
% weitere Argumente haben kann

%% Definition der Zust�nde

mX = x(1,:);
% mS = x(2,:);
V = x(3,:);

%% Messgleichungen

y(1,:) = mX./V + testParam;
y(2,:) = V;