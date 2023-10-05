function y = messgleichung(x)

%% Definition der Zustände

mX = x(1,:);
% mS = x(2,:);
V = x(3,:);

%% Messgleichungen

y(1,:) = mX./V;
y(2,:) = V;