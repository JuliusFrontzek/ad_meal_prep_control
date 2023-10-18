%% Version
% (R2022b) Update 6
% Erstelldatum: 18.10.2023
% Autor: Simon Hellmann

% create summarizing plots for ECC-2024 paper

%% create a bar plot with the improvements in computation time
%   
close 
figBarPlot = figure(); 
noiseVersion = categorical({'add','aug','fully-aug'});
tRun = [54.17 48.42 28.77 5.12; ... %   add                     
        129.57 97.59 56.24 17.17; ... % aug         
        158.06 147.55 66.21 19.36]; %   fully-aug
b = bar(noiseVersion,tRun); % "CData", colorVector
% change color of bars (1 color for 1 setup)
colorPaletteViridis = ["#fc8961", "#b73779", "#51127c", "#000004"]; % 5 colors magma
for k = 1:4
    b(k).FaceColor = colorPaletteViridis(k);
end
ylim([0,160])
ylabel('run time in [s]')
legend('NLP', 'NLP-grad', 'NLP-grad-hess', 'QP','Location','northwest')
fontsize(figBarPlot,20,"point")


%%
y = [2 2 3; 2 5 6; 2 8 9; 2 11 12];
h = bar(y); 


% %% create a plot of nRMSE over tCalc
% %                         1        2           3         4         5           6         7        8         9        10        11       
% eccColorPaletteRMSE = ["#000004","#b73779", "#721f81", "#932b80","#51127c","#3b528b","#287c8e","#fed799","#feb078","#fc8961","#f1605d"]; 
% markerShapes =        {'o',   'square', 'hexagram',  'o',      '*',   'diamond', 'square',  '<',      '>',      'v',      '^'}; 
% %               1               2              3                 4         5                6                  7                8           9                 10           11
% myLabels = {'UKF-sysID','UKF-add','UKF-add-\gamma','SR-UKF','SR-UKF-\gamma','UKF-aug','UKF-fully-aug','cUKF-NLP','cUKF-NLP-grad','cUKF-NLP-grad-hess','cUKF-QP'}; 
% stdSz = 50; % standard size of markers
% magnifier = [1.5, 2, 3];
% myLineWidth = 1.5; 
% 
% % tCalc = [1,2,2,3,4,5,50,45,30,20]; 
% tCalc = [1.7288, 1.5909, 1.4433, 1.4735, 1.5681, 2.7592, 3.5497, 54.1704, 48.4247, 28.7659, 5.1166];
% nRMSE = [0.0158, 0.0261, 0.0261, 0.0118, 0.0118, 0.0542, 0.0666, 0.0105, 0.0105, 0.0105, 0.0105]; 
% 
% % close
% figure_t_RMSE = figure;
% scatter(tCalc(1),nRMSE(1),stdSz,markerShapes{1},'filled','MarkerFaceColor',eccColorPaletteRMSE(1),...
%     'MarkerEdgeColor',eccColorPaletteRMSE(1),'DisplayName',myLabels{1})
% hold on 
% scatter(tCalc(2),nRMSE(2),stdSz*magnifier(3),markerShapes{2},'LineWidth', myLineWidth,...
%     'MarkerEdgeColor',eccColorPaletteRMSE(2),'DisplayName',myLabels{2})
% scatter(tCalc(3),nRMSE(3),stdSz*magnifier(1),markerShapes{3},'filled','MarkerFaceColor',eccColorPaletteRMSE(3),...
%     'MarkerEdgeColor',eccColorPaletteRMSE(3),'DisplayName',myLabels{3})
% scatter(tCalc(4),nRMSE(4),stdSz*magnifier(3),markerShapes{4},'LineWidth',myLineWidth, ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(4),'DisplayName',myLabels{4})
% scatter(tCalc(5),nRMSE(5),stdSz*magnifier(3),markerShapes{5},'LineWidth', 1, ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(5),'DisplayName',myLabels{5})
% scatter(tCalc(6),nRMSE(6),stdSz,markerShapes{6},'filled','MarkerFaceColor',eccColorPaletteRMSE(6), ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(6),'DisplayName',myLabels{6})
% scatter(tCalc(7),nRMSE(7),stdSz*magnifier(2),markerShapes{7},'filled','MarkerFaceColor',eccColorPaletteRMSE(7), ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(7),'DisplayName',myLabels{7})
% scatter(tCalc(8),nRMSE(8),stdSz,markerShapes{8},'filled','MarkerFaceColor',eccColorPaletteRMSE(8), ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(8),'DisplayName',myLabels{8})
% scatter(tCalc(9),nRMSE(9),stdSz,markerShapes{9},'filled','MarkerFaceColor',eccColorPaletteRMSE(9), ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(9),'DisplayName',myLabels{9})
% scatter(tCalc(10),nRMSE(10),stdSz,markerShapes{10},'filled','MarkerFaceColor',eccColorPaletteRMSE(10), ...
%     'MarkerEdgeColor',eccColorPaletteRMSE(10),'DisplayName',myLabels{10})
% scatter(tCalc(11),nRMSE(11),stdSz,markerShapes{11},'filled','MarkerFaceColor',eccColorPaletteRMSE(11),...
%     'MarkerEdgeColor',eccColorPaletteRMSE(11),'DisplayName',myLabels{11})
% xlabel('average run time [s]')
% xlim([0,60]); 
% ylabel('average nRMSE')
% ylim([0,0.08])
% legend()
