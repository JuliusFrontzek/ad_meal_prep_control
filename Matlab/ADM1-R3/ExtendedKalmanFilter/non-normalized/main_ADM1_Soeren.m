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

% add folders to search path: 
addpath("SoerensFiles\")

%% Initialise model
% Load standard model parameters
load('SoerensFiles\Model_data\ADM1_parameters.mat')
% Load experimental data
load('SoerensFiles\Model_data\ADM1_input_data.mat')
% Select model type: ADM1, ADM1_R1, ADM1_R2, ADM1_R3 or ADM1_R4
model = 'ADM1_R3';

%% Run selected model 
systemParameters = system.Variables; 
systemInput = input.ADM1_R3.Variables;
systemInputSS = systemInput(2,:);   % only take the second feeding, which perturbes the system out of its steady state
systemInputSS(1) = 0;   % set initial feeding time to zero
modelParameters = parameters.ADM1_R3.Variables; 

% Solve ODE of mass-based ADM1-R3
% odeFunR3 = @(t,x) ADM1_R3_mass(t,x,systemParameters,systemInputSS,modelParameters); 
odeFunR3 = @(t,x) ADM1_R3_mass_edit(t,x,systemParameters,systemInputSS,modelParameters); 
x0 = initial.ADM1_R3.Variables; % x0
time_range = [0 300];           % time range in days
% simulate system and save results as struct: 
odeObj = ode15s(odeFunR3,time_range,x0);
% note that the struct odeObj has fields x and y instead of t and y. So the
% field x represents the time steps chosen by the integrator

% Calculate model ouput
nSteps = size(odeObj.x,2);
for i = 1:nSteps
    [t(i),y(i,:)] = ADM1_R3_mass_output(odeObj.x(1,i),odeObj.y(:,i),systemParameters,modelParameters);
end

%% Set model output
output=output.(model);
output{1:nSteps,:} = [t' y];

%% Plot model output
plot(output{:,1},output{:,2:end});
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
Snh4SS = xSS(4) - xSS(15);

% compute TS (VS is missing because no ash considered so far): 
rho = 1000;     % mass densitiy [g/l]
TS = 1 - xSS(5)/rho; 

% re-order outputs: 
ySS = [volFlow, pch4, pco2, pH, Snh4SS, TS, nan, SacSS]; % leave VS nan for now

%% save results: 
resultsSoeren = struct; % empty struct
resultsSoeren.x0 = x0; 
resultsSoeren.xSS = xSS; 
resultsSoeren.ySS = ySS; 
resultsSoeren.input = systemInputSS; 

% create sub-folder (if non-existent yet) and save results there
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedOutput');
mkdir(pathToResults);   % create subfolder (gives warning if exists yet)
fileName = 'SteadyState_ADM1-R3_Soeren.mat'; 
save(fullfile(pathToResults,fileName), 'resultsSoeren')

%% Clear variables 
clearvars i num_t t y l time_range