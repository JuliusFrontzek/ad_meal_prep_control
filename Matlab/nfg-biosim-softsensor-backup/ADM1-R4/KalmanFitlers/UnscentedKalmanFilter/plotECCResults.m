%% Version
% (R2022b) Update 6
% Erstelldatum: 18.10.2023
% Autor: Simon Hellmann

close all 
clear 
clc

% create summarizing plots for ECC-2024 paper

%% create a bar plot with the improvements in computation time
%   

x      = 100;   % Screen position
y      = 200;   % Screen position
width  = 600;   % Width of figure (by default in pixels)
heightFactor = 0.3;
height = width*heightFactor; % Height of figure
figBarPlot = figure('Position', [x y width height]);

noiseVersion = categorical({'add','aug','fully-aug'});
tRun = [54.17 48.42 28.77 5.12; ... %   add                     
        129.57 97.59 56.24 17.17; ... % aug         
        158.06 147.55 66.21 19.36]; %   fully-aug
b = bar(noiseVersion,tRun); % "CData", colorVector
% change color of bars (1 color for 1 setup)
colPalMagma = ["#fc8961", "#b73779", "#51127c", "#000004"]; % 5 colors magma
for k = 1:4
    b(k).FaceColor = colPalMagma(k);
end
ylim([0,165])
ylabel('run time in [s]', 'Interpreter','latex')
% set up legend:
legend('NLP', 'NLP-grad', 'NLP-grad-hess', 'QP','Position',[0.18 0.6 0.2 0.2],'Interpreter','latex')
fontsize(figBarPlot,14,"point")

% prepare plot to be saved as pdf:
set(gcf, 'PaperPosition', [0 0 16 16*heightFactor]); % Position plot at left hand corner with specified width and height.
set(gcf, 'PaperSize', [16 16*heightFactor]); % Set the paper to have specified width and height (excess is cut off)

% % create sub-folder (if non-existent yet) and save plot there:
% currPath = pwd; 
% pathToResults = fullfile(currPath,'generatedPlots');
% if ~exist(pathToResults, 'dir')
%     mkdir(pathToResults)
% end
% plotName = 'barPlotcUKF'; 
% fileName = fullfile(pathToResults,plotName); 
% saveas(gcf, fileName, 'pdf') %Save figure


%% create a plot of nRMSE over tCalc
%                         1          2         3           4         5         6        7         8       
eccColorPaletteRMSE = ["#fcfdbf","#febb81", "#f8765c", "#d3436e","#982d80","#5f187f","#221150","#000004"]; 
markerShapes =        {'diamond', '>',        '^',      '<',     'square',    'o',  'pentagram','hexagram'}; 
% markerShapes =        {'o',   'square', 'hexagram',  'o',      '*',   'diamond', 'square',  '<',      '>',      'v',      '^'}; 
%               1                     2                    3               4               5          6             7                8                   
myLabels = {'UKF-SR-$\gamma$','UKF-add-$\gamma$','UKF-aug-$\gamma$','UKF-fully-aug','UKF-sysID','cUKF-add-QP','cUKF-aug-QP','cUKF-fully-aug-QP'}; 
stdSz = 50; % >standard size of markers
newStdSz = 10; % new standard size of markers
magnifier = [1.5, 2, 3];
myLineWidth = 1.5; 

% tCalc = [1,2,2,3,4,5,50,45,30,20]; 
tCalc = [2.18,  2.15,   3.84,   4.14,   1.73,   5.12,   17.17,  19.36];
nRMSE = [0.0118,0.0118, 0.0116, 0.0125, 0.0158, 0.0105, 0.0206, 0.0194]; 

x      = 100;   % Screen position
y      = 200;   % Screen position
width  = 700;   % Width of figure (by default in pixels)
heightFactor2 = 0.6;
height = width*heightFactor2; % Height of figure
figure_t_RMSE = figure('Position', [x y width height]);
% create an additional (empty) subplot window with just for the legend: 
h1 = subplot(2,1,1); 
plot(1, nan, 1, nan, 'r'); % just plot anything for fun (will be invisible later)
set(h1, 'OuterPosition', [0,0.5, 1, .21]); % [x0, y0, width, height]
set(h1, 'Visible', 'off'); % make plot invisible
h2 = subplot(2,1,2); 
% scatter(tCalc(1),nRMSE(1),stdSz,markerShapes{1},'filled','MarkerFaceColor',eccColorPaletteRMSE(1),...
%     'MarkerEdgeColor',eccColorPaletteRMSE(1),'DisplayName',myLabels{1})
plot(tCalc(1),nRMSE(1),'w','DisplayName',myLabels{1},'Marker',markerShapes{1},...
    'MarkerSize',newStdSz*magnifier(1), 'MarkerFaceColor',eccColorPaletteRMSE(1),'MarkerEdgeColor','k','LineWidth',myLineWidth);
hold on 
plot(tCalc(2),nRMSE(2),'w','DisplayName',myLabels{2},'Marker',markerShapes{2},...
    'MarkerSize',newStdSz, 'MarkerFaceColor',eccColorPaletteRMSE(3)); 
plot(tCalc(3),nRMSE(3),'w','DisplayName',myLabels{3},'Marker',markerShapes{3},...
    'MarkerSize',newStdSz, 'MarkerFaceColor',eccColorPaletteRMSE(4)); 
plot(tCalc(4),nRMSE(4),'w','DisplayName',myLabels{4},'Marker',markerShapes{4},...
    'MarkerSize',newStdSz, 'MarkerFaceColor',eccColorPaletteRMSE(5)); 
plot(tCalc(5),nRMSE(5),'w','DisplayName',myLabels{5},'Marker',markerShapes{5},...
    'MarkerSize',newStdSz, 'MarkerFaceColor',eccColorPaletteRMSE(2)); 
plot(tCalc(6),nRMSE(6),'w','DisplayName',myLabels{6},'Marker',markerShapes{6},...
    'MarkerSize',newStdSz, 'MarkerFaceColor',eccColorPaletteRMSE(6)); 
plot(tCalc(7),nRMSE(7),'w','DisplayName',myLabels{7},'Marker',markerShapes{7},...
    'MarkerSize',newStdSz*magnifier(1), 'MarkerFaceColor',eccColorPaletteRMSE(7)); 
plot(tCalc(8),nRMSE(8),'w','DisplayName',myLabels{8},'Marker',markerShapes{8},...
    'MarkerSize',newStdSz*magnifier(1), 'MarkerFaceColor',eccColorPaletteRMSE(8)); 
xlabel('run time [s]','Interpreter','latex')
xlim([0,20]); 
ylabel('NRMSE','Interpreter','latex')
ylim([0,0.025])
fontsize(figure_t_RMSE,15,'points')
leg = legend('Interpreter','latex','Position',[0.38 0.47 0.5 0.2]); 
leg.NumColumns = 2; %   2 column legend

% prepare plot to be saved as pdf:
heightFactorExport = 0.42; 
set(gcf, 'PaperPosition', [0 0 18 18*heightFactor2]); % Position plot at left hand corner with specified width and height
set(gcf, 'PaperSize', [18 18*heightFactorExport]); % Set the paper to have specified width and height (excess is cut off).
% create sub-folder (if non-existent yet) and save plot there:
currPath = pwd; 
pathToResults = fullfile(currPath,'generatedPlots');
if ~exist(pathToResults, 'dir')
    mkdir(pathToResults)
end
plotName = 'nRMSE_vs_runTime'; 
fileName = fullfile(pathToResults,plotName); 
saveas(gcf, fileName, 'pdf') %Save figure