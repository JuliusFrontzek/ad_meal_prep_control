clc
clear all
close all

% define symbolic variables:
x = sym('x', [12 1]);   % states
xx = sym('x',[3,1]);    % states for smaller model
syms u real             % input
th = sym('th', [1,6]);  % 6 unknown time-variant parameters (theta)
c = sym('c', [1,21]);   % known & time-invariant parameters
xi = sym('xi', [10,1]); % 10 (assumed) known inlet concentrations: 

% define petersen matrix symbolically:
a = sym('a', [9,5]);
% a = sym('a%d%d', [9,5]);    % for nice matrix entries a11, a12, etc (instead of a1_1, etc.)
aSimple = sym('a%d%d',[3,2]);  % for proper a11, a12, etc. in matrix entries

% beziehe symbolische, vektorwertige Funktion:
dynamics = BMR4_AB_frac_ode_sym(x, u, xi, th, c, a);
dynamicsTest = odeFunSymTrial(xx,u,aSimple); 

% berechne Jacobimatrizen symbolisch (nur zum Test)
A = jacobian(dynamics, x); 
B = jacobian(dynamics, u); 

% mache daraus numerische Funktion. Beachte, dass du die abhängigen
% Variablen definieren und dann bei der Verwendung auch beibehalten musst.
f = matlabFunction(dynamics, 'Vars', {x, u, xi, th, c, a}); 
fSimple = matlabFunction(dynamicsTest, 'Vars', {xx,u,aSimple}); 
Am = matlabFunction(A, 'Vars', {x, u, th, c, a}); 

% teste numerische Funktion, indem du numerische Eingangs-Variablen
% definierst (beachte, dass Dimensionen mit symbolischen Variablen 
% übereinstimmen müssen):
xn = ones(12, 1); 
xin = ones(10,1); 
uN = 1.0;
tn = ones(1, 6); 
cn = ones(1, 21); 
aN = ones(9, 5); 

% define numeric variables:
xxN = ones(3,1);
xx0 = ones(3,1); 
aNSimple = ones(3,2); 

% test if function handles works: 
f(xn, uN, xin, tn, cn, aN)
Am(xn, uN, tn, cn, aN)
fSimple(xxN, uN, aNSimple)

% create function handles:
myOdeFunX = @(x) fSimple(x, uN, aNSimple);    % only x-dependent
myOdeFunTX = @(t,x) fSimple(x, uN, aNSimple); % x and t-dependent

% simulate dynamic system: 
tSpan = [1,10]; 
[tSim,xSim] = ode45(myOdeFunTX,tSpan,xx0);

%% plot results
figure
plot(tSim,xSim)
legend('x1', 'x2', 'x3')

%% call other function with function handle as input
callFunctionHandle(fSimple,xxN,uN,aNSimple); 
