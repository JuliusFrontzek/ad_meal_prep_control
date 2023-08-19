%% Version
% (R2022b) Update 5
% Erstelldatum: 13.07.2023
% Autor: Simon Hellmann

function plotOutputs(flagModel,flagFrac,yClean,yMeas,feedVolFlow,tGrid,tEvents)
% plots outputs based on specified model

% flagModel -   3: ADM1-R3; 4: ADM1-R4
% flagFrac -    0: no -frac (only 1 CH-fraction); 1: -frac (second CH-fraction)
% yClean -      noise-free model outputs
% yMeas -       noisy synthetic measurements
% feedVolFlow - feeding volume flow [l/d]
% tGrid -       fine online time grid
% tEvents -     time instance when change in feeding happens (on/off)

fracAddOn = '-frac'; 
colorPaletteHex = ["#003049","#d62828","#f77f00","#02C39A","#219ebc"]; 
% plot and compare the clean results with noisy measurements: 
figure()

switch flagModel
    case 3
        modelName = 'R3'; 
        if flagFrac == 1
            % add extension '-frac' to model name:
            modelName = [modelName,fracAddOn]; 
        end
        
        % gas volume flow: 
        subplot(4,2,1)
        scatter(tGrid,yMeas(:,1)/24,'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,1)/24,'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
        ylabel('gas vol flow [l/h]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        
        % pch4: 
        subplot(4,2,2)
        scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
        hold on
        plot(tGrid,yClean(:,2),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('p_{ch4} [bar]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % pco2:
        subplot(4,2,3)
        scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,3),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('p_{co2} [bar]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % pH:  
        subplot(4,2,4)
        scatter(tGrid,yMeas(:,4),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
        hold on; 
        plot(tGrid,yClean(:,4),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('pH value [-]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % Snh4:  
        subplot(4,2,5)
        scatter(tGrid,yMeas(:,5),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,5),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('NH4-N [g/l]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % TS:  
        subplot(4,2,6)
        scatter(tGrid,yMeas(:,6),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,6),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('total solids [-]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % VS:  
        subplot(4,2,7)
        scatter(tGrid,yMeas(:,7),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
        hold on; 
        plot(tGrid,yClean(:,7),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('volatile solids [-]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % Sac:  
        subplot(4,2,8)
        scatter(tGrid,yMeas(:,8),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
        hold on; 
        plot(tGrid,yClean(:,8),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('acetic acid [g/l]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')

    case 4
        modelName = 'R4'; 
        if flagFrac == 1
            % add extension '-frac' to model name:
            modelName = [modelName,fracAddOn]; 
        end

        % gas volume flow: 
        subplot(3,2,1)
        scatter(tGrid,yMeas(:,1)/24,'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,1)/24,'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5); 
        ylabel('gas vol flow [l/h]')
        ylim([0,16])
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        
        % pch4: 
        subplot(3,2,2)
        scatter(tGrid,yMeas(:,2),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5);
        hold on
        plot(tGrid,yClean(:,2),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('p_{ch4} in bar')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % pco2:
        subplot(3,2,3)
        scatter(tGrid,yMeas(:,3),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,3),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('p_{co2} in bar')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % SIN:  
        subplot(3,2,4)
        scatter(tGrid,yMeas(:,4),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,4),'r',...
             'LineWidth',1.2,'DisplayName','clean')
        ylabel('inorg. nitrogen in g/L')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % TS:  
        subplot(3,2,5)
        scatter(tGrid,yMeas(:,5),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5); 
        hold on; 
        plot(tGrid,yClean(:,5),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('total solids [-]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')
        
        % VS:  
        subplot(3,2,6)
        scatter(tGrid,yMeas(:,6),'DisplayName','noisy',...
                'Marker','.', 'Color', colorPaletteHex(1), 'LineWidth',1.5)
        hold on; 
        plot(tGrid,yClean(:,6),'DisplayName','clean',...
             'LineStyle','-.', 'Color', colorPaletteHex(2), 'LineWidth',1.5);
        ylabel('volatile solids [-]')
        xlabel('time [d]')
        yyaxis right
        stairs(tEvents, feedVolFlow/24, 'DisplayName','feeding',...
               'LineStyle','-', 'Color', colorPaletteHex(4), 'LineWidth',1.5); 
        set(gca, "YColor", 'k')     % make right y-axis black 
        ylabel('feed vol flow [l/h]')
        legend('Location','NorthEast'); 
        set(gca, "YColor", 'k')

end % switch

sgtitle(['Clean and noisy simulation outputs ',modelName])

end % fun