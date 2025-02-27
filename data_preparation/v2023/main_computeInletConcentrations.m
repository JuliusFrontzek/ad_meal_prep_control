%% Version: Matlab R2022b, Update 6
% Author: Simon Hellmann
% first created: 2024-08-25

% compute inlet concentrations from lab measurements of substrates
close all
clear
clc

% specify substrate type: 
agriculturalSubstrates = {'Maissilage', 'Grassilage', 'Getreidestroh', 'Zuckerruebensilage', ...
                          'Rinderguelle', 'Schweineguelle', 'HTK'}; 
% BMP values (source: pig manure/HTK/Stroh: KTBL Gasausbeute in
% landwirtsch. Biogasanlagen (2021), Mais/Gras/Rindergülle: DBFZ batch tests, 
% Rübe: Vazifehkhoran (2016))
bmpAgriculturalSubstrates = [357, 372, 240, 389, 246, 230, 272]; % [L_N/kg_oTS]

% compute inlet concentrations for all 7 agricultural substrates: 
nSubstrates = numel(bmpAgriculturalSubstrates); 
% allocate memory: 
xInAgrSubstrates = nan([18,nSubstrates]); 
TSAgrSubstrates = nan([nSubstrates,1]);
bk_numbers = cell(50,nSubstrates); 
nominal_values = nan([nSubstrates,5]); 
S_macroNutrients_tables = struct;% struct with tables of all realized macro nutrients

for k = 1:nSubstrates       % iterate over all substrates
    substrateNumber = k;    % 1...7 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    BMP = bmpAgriculturalSubstrates(substrateNumber); 
    BMP_stoich = 420;   % max. BMP for agricultural substrates [L_N/kg_FoTS], Weißbach (2009)
    
    %% get all tables (respect the appropriate sheets): 
    addpath('../../data/data_in/dbfz_substrate_data_2017_2025/')
    T2018 = readtable('2018_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2019 = readtable('2019_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2020 = readtable('2020_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2021 = readtable('2021_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2022 = readtable('2022_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    colNamesOrig = T2022.Properties.VariableNames;      % get column names (the same for all tables)
    colNamesEdit = {'Probeneingangsnummer','TS','water','NH4N','XA','XP','XL','XC','AcMeas','dilutionFactor','Ac','Kommentar'}; 

    % transform tables into arrays without Probeneingangsnummer and Kommentar, unify to 1 array:
    labDataRaw2018 = T2018(:,2:11); 
    labDataRaw2019 = T2019(:,2:11); 
    labDataRaw2020 = T2020(:,2:11); 
    labDataRaw2021 = T2021(:,2:11); 
    labDataRaw2022 = T2022(:,2:11);
    labDataRaw = [labDataRaw2018;labDataRaw2019;labDataRaw2020;labDataRaw2021;labDataRaw2022]; 
    
    % remove all rows that contain nans (i.e. are not complete): 
    mask_not_nan = ~any(isnan(table2array(labDataRaw)),2); 
    labDataEffPre = labDataRaw(mask_not_nan,:); % effective lab data
    % remove all rows containing negative values: 
    mask_positive = all(table2array(labDataEffPre) > 0,2); % ignore Probeneingangsnummer 
    labDataEff = labDataEffPre(mask_positive,:);    
    
    % insert special inlet specifications for silages acc. to following sources:
    % acetic acids based on Weissbach & Strubelt (2008), use equivalent acid values.
    % pH based on average value in DBFZ lab analyses of substrates. 
    % NH4-N from Schlattmann (2011), Tab. 26; and Demirel & Scherer (2008),
    % Tab. 1 for sugar beet silage:
    switch substrate 
        case "Maissilage"
            substrate_lit_values.pH = 3.8;      % [-] (Weißbach 2008)
            substrate_lit_values.NH4N = 0.764;  % [g/L]
            substrate_lit_values.ACeq = 10.32;  % [g/kg FM] = [g/L]
        case "Grassilage"
            substrate_lit_values.pH = 4.8;      % Weißbach (2008)
            substrate_lit_values.NH4N = 1.574;
            substrate_lit_values.ACeq = 10.44;
        case "Zuckerruebensilage"
            substrate_lit_values.pH = 3.9;      % Weißbach (2008)
            substrate_lit_values.NH4N = 0.070;
            substrate_lit_values.ACeq = 8.17;
        case "Rinderguelle"
            substrate_lit_values = []; % refresh blank
            substrate_lit_values.pH = 8.54;     % Fisgativa (2020)
        otherwise 
            substrate_lit_values = []; 
    end
    
    % allocate memory:
    nSamples = size(labDataEff,1);  % # complete samples per substrate
    xInMat = nan(18,nSamples);  
    TSVec = nan(1,nSamples); 

    % compute xIn for ADM1-R3-frac from all samples in labDataEff 1 by 1:
    for kk = 1:nSamples
        [xInMat(:,kk),TSVec(kk)] = computeX_inR3FracFromLabMeasurements(labDataEff{kk,:},BMP,BMP_stoich,substrate_lit_values); 
    end 
    
    % compute mean of all xIn_i and TS from the for loop and save in
    % separate structs: 
    xIn_mean = mean(xInMat,2);
    TS_mean = mean(TSVec);
    xInAgrSubstrates(:,k) = xIn_mean; 
    TSAgrSubstrates(k) = TS_mean;  
%     % save all full realizations of inlet concentrations for later analysis:
%     xInMatSubstrates.(substrate) = xInMat; 
    % create table of all available data points of 3 macronutrients:
    Xch = xInMat(6,:)';
    Xpr = xInMat(8,:)';
    Xli = xInMat(9,:)';
    S_macroNutrients_tables.(substrate) = table(Xch, Xpr, Xli); 
    
    %% save mean nominal values in table
    nominal_values(k,1) = TS_mean*100; % [% FM]
    nominal_values(k,2) = BMP;
    raw_ash_pro_li = mean(labDataEff(:,4:6)); % [g/kg TS]
    nominal_values(k,3:5) = table2array(raw_ash_pro_li)/10; % [% TS] 

    %% get BK numbers of surviving samples: 
    bk_numbers_2018 = T2018(:,1); 
    bk_numbers_2019 = T2019(:,1); 
    bk_numbers_2020 = T2020(:,1); 
    bk_numbers_2021 = T2021(:,1); 
    bk_numbers_2022 = T2022(:,1);
    bk_numbers_init = [bk_numbers_2018;bk_numbers_2019;bk_numbers_2020;bk_numbers_2021;bk_numbers_2022]; 
    
    % apply the same filters as above: 
    bk_numbers_pre = bk_numbers_init(mask_not_nan,:); % effective lab data
    bk_numbers_eff = bk_numbers_pre(mask_positive,:); 
    
    bk_numbers(1:nSamples,k) = bk_numbers_eff.Variables; 

end % for

%% build and save tables: 
bk_numbers_tab = cell2table(bk_numbers, 'VariableNames',agriculturalSubstrates); % turn into table

% nominal values for uncertain scenarios:
float_tab_nominal = array2table(nominal_values, 'VariableNames', ...
    {'TS [%FM]','BMP [L-norm CH4/kg VS]','XA [%DM]','XP [%DM]','XL [%DM]'});
% get table from relevant substrate names: 
string_tab = cell2table(agriculturalSubstrates', 'VariableNames', {'Substrate'});
% Combine the two tables horizontally
nominal_values_tab = [string_tab,float_tab_nominal];

% resulting mean inlet concentrations [g/L]:
inlet_concentrations_tab = array2table(xInAgrSubstrates, 'VariableNames', agriculturalSubstrates);

% save
writetable(bk_numbers_tab,'../../data/data_out/v2023/used_bk_numbers.xlsx');
writetable(nominal_values_tab,'../../data/data_out/v2023/nominal_values.xlsx');
writetable(inlet_concentrations_tab,'../../data/data_out/v2023/nominal_inlet_concentrations.xlsx');

% save macronutrient concentrations as excel with different tabs: 
fileName = '../../data/data_out/v2023/macro_nutrients_concentrations.xlsx';
for k = 1:nSubstrates
    % get current table: 
    substrate = agriculturalSubstrates{k}; 
    curr_macro_nutrients_tab = S_macroNutrients_tables.(substrate);  

    % save table as new tab:
    writetable(curr_macro_nutrients_tab,fileName, 'Sheet',substrate); 
end

%% compute descriptive statistics of all substrates 
% (# samples, min, max, mean, median, q25, q75, IQR)
% create a struct that has a field for each substrate; each field contains
% a table with the statistical values (columns) for ch/pr/li (rows); then
% saved as excel file with different tabs
macroNutrients_descr_stats = struct(); % allocate memory

for kk = 1:nSubstrates
    macros_descr_stats_kk = nan(8,3); % allocate memory (8 measures, 3 macro nutrients)
    
    % get data: 
    substrateNumber = kk;    % 1...7 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    macroNutrientsValues = table2array(S_macroNutrients_tables.(substrate)); % macronutrient values for one substrate
    
    % compute statistical measures: 
    nSamples = size(macroNutrientsValues,1);
    minVals = min(macroNutrientsValues); 
    maxVals = max(macroNutrientsValues); 
    meanVals = mean(macroNutrientsValues);
    medianVals = median(macroNutrientsValues); 
    q25 = quantile(macroNutrientsValues,0.25); % 25 percent quantile
    q75 = quantile(macroNutrientsValues,0.75); % 75 percent quantile
    iqr = q75 - q25; % inter quartile range

    % put statistical measures into array: 
    macros_descr_stats_kk(1,:) = repmat(nSamples,[1,3]); 
    macros_descr_stats_kk(2:end,:) = [minVals;maxVals;meanVals;medianVals;q25;q75;iqr];  

    % turn into table and save in struct: 
    macro_col_names = {'X_ch','X_pr','X_li'}; 
    macro_row_names = {'n_samples','min','max','mean','median','q25','q75','iqr'}; 
    macros_descr_stats_tab_kk = array2table(macros_descr_stats_kk, ...
        'VariableNames',macro_col_names, 'RowNames',macro_row_names);
    macroNutrients_descr_stats.(substrate) = macros_descr_stats_tab_kk; 

    % save as excel: 
    fileName = '../../data/data_out/v2023/macro_nutrients_descriptive_statistics.xlsx'; 
    writetable(macros_descr_stats_tab_kk,fileName, 'Sheet',substrate, 'WriteRowNames',true)
end

%% draw boxplots of macro nutrients
close all
do_mpc_substrates = agriculturalSubstrates([1,2,4,5]); 

figure
tiledlayout(3,1)

% carbohydrates: 
nexttile
hold on
for k = 1:4
    macro_table = S_macroNutrients_tables.(do_mpc_substrates{k}); 
    carbs = table2array(macro_table(:,1)); 
    boxplot(carbs, 'orientation','horizontal','Positions',k)
end
xlabel('\xi_{ch} [g/L]'); % Label for the x-axis
ylabel('Substrate');
yticks(1:4); 
yticklabels({'maize','grass','sugar beet','manure'}); % Custom y-tick labels

% proteins: 
nexttile
hold on
for k = 1:4
    macro_table = S_macroNutrients_tables.(do_mpc_substrates{k}); 
    proteins = table2array(macro_table(:,2)); 
    boxplot(proteins, 'orientation','horizontal','Positions',k)
end
xlabel('\xi_{pr} [g/L]'); % Label for the x-axis
ylabel('Substrate');
yticks(1:4); 
yticklabels({'maize','grass','sugar beet','manure'}); % Custom y-tick labels

% lipids: 
nexttile
hold on
for k = 1:4
    macro_table = S_macroNutrients_tables.(do_mpc_substrates{k}); 
    lipids = table2array(macro_table(:,3)); 
    boxplot(lipids, 'orientation','horizontal','Positions',k)
end
xlabel('\xi_{li} [g/L]'); % Label for the x-axis
ylabel('Substrate');
yticks(1:4); 
yticklabels({'maize','grass','sugar beet','manure'}); % Custom y-tick labels
sgtitle('Distribution of computed macronutrient inlet concentrations')
% save: 
pdf_path_name = '../../results/plots/boxplots_influent_macronutrients.pdf'; 
exportgraphics(gcf,pdf_path_name, 'ContentType','vector');
