% -------------------------------------------------------------------------
%
% Run different simplifications of a mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% Example: Anaerobic co-digestion of maize silage and cattle manure
%
% Implementation in Matlab 
% R2019b (The MathWorks, Inc.)
%
% Matlab ODE function
%
% Version 1.0
%
% https://github.com/soerenweinrich/ADM1
%
% Copyright (c) 2021 Sören Weinrich
% E-Mail: soeren.weinrich@dbfz.de
%
% Additional information (citation example):
%
% Weinrich, S.; Nelles, M. (2021).
% Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
% Model development and stoichiometric analysis. Bioresource Technology. 
% In press. https://doi.org/10.1016/j.biortech.2021.125124.
%
% -------------------------------------------------------------------------

close all
clear all
clc

%% Initialise model
% Load standard model parameters
load('Model_data\ADM1_parameters.mat')
% Load experimental data
load('Model_data\ADM1_input_data.mat')
% Select model type: ADM1, ADM1_R1, ADM1_R2, ADM1_R3 or ADM1_R4
model = 'ADM1_R3';

%% Run selected model 
systemParameters = system.Variables; 
p_atm = 1.0133; % atmospheric pressure [bar] 
systemParameters(3) = p_atm; % XY: Achtung: Typo im GitHub 
systemInput = input.ADM1_R3.Variables;
systemInputSS = systemInput(2,:);   % only take the second feeding, which perturbes the system out of its steady state
systemInputSS(1) = 0;   % set initial feeding time to zero
modelParameters = parameters.ADM1_R3.Variables; 

% change scale from lab to FBGA: 
scaleFactor = 1; 
systemParameters(1) = scaleFactor*systemParameters(1); 
systemParameters(2) = 0.1*systemParameters(1); % maintain Vg/Vl ratio
systemInputSS(2) = scaleFactor*systemInputSS(2); 

% Solve ODE of mass-based ADM1-R3
% odeFunR3 = @(t,x) ADM1_R3_mass(t,x,systemParameters,systemInputSS,modelParameters); 
odeFunR3 = @(t,x) ADM1_R3_mass_edit(t,x,systemParameters,systemInputSS,modelParameters); 
x0 = initial.ADM1_R3.Variables; % x0
tSS = 300; 
time_range = [0 tSS];           % time range in days
% simulate system and save results as struct: 
odeObj = ode15s(odeFunR3,time_range,x0);
% note that the struct odeObj has fields x and y instead of t and y. So the
% field x represents the time steps chosen by the integrator

% Calculate model ouput
nSteps = size(odeObj.x,2);
for i = 1:nSteps
    [t(i),y(i,:)] = ADM1_R3_mass_output(odeObj.x(1,i),odeObj.y(:,i),systemParameters,modelParameters);
end

%% Set model output as table
output=output.(model);      % the struct output is loaded at top of script
output{1:nSteps,:} = [t' y];% add simulation results to pre-existing table
tVec = output{:,1};             % extract table content as array
modelOutput = output{:,2:end};  % dito

%% Plot model output
plot(tVec,modelOutput);
% plotbrowser('on');
% Set legend
l = legend(output.Properties.VariableNames(2:end));
set(l,'Interpreter','none','visible','off');
% Set title
t = title(['Simulation results of the ',model]);
set(t,'Interpreter','none');
% Set axis labels
xlabel('Time [d]');
ylabel('Model output'); 

%% extract, interpolate and plot measurements only
measurements = modelOutput(:,18:end); % [S_nh4,S_co2,Phi,SHPlus,pH,p_co2,p_ch4,p_gas,q_gas]
% set fixed time vector to interpolate along: 
dt = 1;     % sample time [d]
tVecIntoSS = 0:dt:tSS;  
measurementsInt = interp1(tVec,measurements,tVecIntoSS); 

figure()
plot(tVecIntoSS,measurementsInt); 
legend('Snh4','Sco2','Phi','SHPlus','pH','pco2','pch4','pgas','qgas')

%% determine steady-state outputs in sequence acc. to arXiv-paper: 
ySSPre = output{end,19:end}; 
% take that vector appart:
Snh4 = ySSPre(1); 
Sco2 = ySSPre(2); 
phi = ySSPre(3); 
ShPlus = ySSPre(4); 
pH = ySSPre(5); 
pch4 = ySSPre(6); 
pco2 = ySSPre(7); 
pgas = ySSPre(8); 
volFlow = ySSPre(9); 

% obtain required states; 
xSS = odeObj.y(:,end);  % steady state
SacSS = xSS(1); 
SINSS = xSS(4); 

% compute TS (VS is missing because no ash considered so far): 
rho = 1000;     % mass densitiy [g/l]
TS = 1 - xSS(5)/rho; 

% re-order outputs: 
ySS = [volFlow, pch4, pco2, pH, SINSS, TS, nan, SacSS]; % leave VS nan for now

%% save results: 
resultsSoeren = struct; % empty struct
resultsSoeren.x0 = x0; 
resultsSoeren.xSS = xSS; 
resultsSoeren.ySS = ySS; 
resultsSoeren.tVecIntoSS = tVecIntoSS; % sample time of 1 d
resultsSoeren.y = measurementsInt; % measurements [S_nh4,S_co2,Phi,SHPlus,pH,p_co2,p_ch4,p_gas,q_gas] interpolated to time grid of tVecIntoSS
resultsSoeren.input = systemInputSS; 

save('SteadyState_ADM1-R3_Soeren.mat', 'resultsSoeren'); 

%% Clear variables 
clearvars i num_t t y l time_range