close all

% Sample data
data1 = rand(100, 1) + 1; % Mais
data2 = rand(100, 1) + 2; % Gras
data3 = rand(100, 1) + 3; % Rübe

% Combine data into a cell array
data = {data1, data2, data3};

% Create a figure
figure;

% Loop through each dataset to create stacked boxplots
hold on; % Keep the current plot
for i = 1:length(data)
    % Create a horizontal boxplot for each dataset
    boxplot(data{i}, 'Orientation', 'horizontal', 'Positions', i);
end
hold off;

% Set labels and title
xlabel('Kohlenhydrate');
ylabel('Substrate');
yticks(1:length(data)); % Set y-ticks to match the number of groups
yticklabels({'Mais', 'Gras', 'Rübe'}); % Custom y-tick labels
title('Stacked Horizontal Boxplots');