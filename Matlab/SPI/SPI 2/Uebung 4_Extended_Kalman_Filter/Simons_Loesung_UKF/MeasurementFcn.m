function y = MeasurementFcn(x,testParam)
% testParam wurde nur hinzugefügt, um zu testen, ob die Messgleichung auch
% weitere Argumente haben kann

%% Definition der Zustände

mX = x(1,:);
% mS = x(2,:);
V = x(3,:);

%% Messgleichungen

y(1,:) = mX./V + testParam;
y(2,:) = V;