close all

% Sample data
data = rand(100, 3); % Generate random data for three groups

% Create a boxplot
h = boxplot(data, 'ColorGroup', [1; 2; 3]); % Assign color groups

% Set colors for each box
colors = [1 0 0; 0 1 0; 0 0 1]; % Red, Green, Blue
for i = 1:width(h)
    set(h(i), 'Color', colors(i, :)); % Set color for each box
end

xlabel('Values');
ylabel('Groups');
title('Boxplot with Different Colors');