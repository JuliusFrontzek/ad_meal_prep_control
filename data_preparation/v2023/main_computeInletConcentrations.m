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
% create structs with available substrates as field names for all desired results: 
xInAgrSubstrates = struct;  % for the average inlet concentration vector
macroNutrientsTabs = struct;% struct with tables of all realized macro nutrients
TSAgrSubstrates = struct;   % for TS values per substrate
bk_numbers = cell(50,nSubstrates); % allocate memory

for k = 1:nSubstrates       % iterate over all substrates
    substrateNumber = k;    % 1...7 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    BMP = bmpAgriculturalSubstrates(substrateNumber); 
    BMP_stoich = 420;   % max. BMP for agricultural substrates [L_N/kg_FoTS], Weißbach (2009)
    
    % get all tables (respect the appropriate sheets): 
    addpath('../../data/data_in/dbfz_substrate_data_2017_2025/')
    T2018 = readtable('2018_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2019 = readtable('2019_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2020 = readtable('2020_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2021 = readtable('2021_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2022 = readtable('2022_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    colNamesOrig = T2022.Properties.VariableNames;      % get column names (the same for all tables)
    colNamesEdit = {'Probeneingangsnummer','TS','water','NH4N','XA','XP','XL','XC','AcMeas','dilutionFactor','Ac','Kommentar'}; 

    % transform tables into arrays without Probeneingangsnummer and Kommentar, unify to 1 array:
    labDataRaw2018 = T2018{:,2:11}; 
    labDataRaw2019 = T2019{:,2:11}; 
    labDataRaw2020 = T2020{:,2:11}; 
    labDataRaw2021 = T2021{:,2:11}; 
    labDataRaw2022 = T2022{:,2:11};
    labDataRaw = [labDataRaw2018;labDataRaw2019;labDataRaw2020;labDataRaw2021;labDataRaw2022]; 
    
    % remove all rows that contain nans (i.e. are not complete): 
    mask_not_nan = ~any(isnan(labDataRaw),2); 
    labDataEffPre = labDataRaw(mask_not_nan,:); % effective lab data
    % remove all rows containing negative values: 
    mask_positive = all(labDataEffPre > 0,2); % ignore Probeneingangsnummer 
    labDataEff = labDataEffPre(mask_positive,:); 
    
    % allocate memory:
    nSamples = size(labDataEff,1);  % # complete samples per substrate
    xInMat = nan(18,nSamples);  
    TSVec = nan(1,nSamples);    
    
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

    % pick rows from labDataEff 1 by 1 and compute xIn, save in xInMat:
    for kk = 1:nSamples
        labData = labDataEff(kk,:); 
        [xInMat(:,kk),TSVec(kk)] = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich,substrate_lit_values); 
    end % for
    
    % compute mean of all xIn_i and TS from the for loop and save in
    % separate structs: 
    xInSubstrate = mean(xInMat,2);
    TSSubstrate = mean(TSVec);
    xInAgrSubstrates.(substrate) = xInSubstrate; 
    TSAgrSubstrates.(substrate) = TSSubstrate;  
%     % save all full realizations of inlet concentrations for later analysis:
%     xInMatSubstrates.(substrate) = xInMat; 
    % create table of all available data points of 3 macronutrients:
    Xch = xInMat(6,:)';
    Xpr = xInMat(8,:)';
    Xli = xInMat(9,:)';
    macroNutrientsTabs.(substrate) = table(Xch, Xpr, Xli); 

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

fileName = '../../data/data_out/v2023/xIn_TS_AgrSubstrates.mat'; 
save(fileName, 'xInAgrSubstrates', 'TSAgrSubstrates', 'macroNutrientsTabs')
bk_numbers_tab = cell2table(bk_numbers, 'VariableNames',fieldnames(xInAgrSubstrates)); % turn into table
writetable(bk_numbers_tab,'../../data/data_out/v2023/used_bk_numbers.xlsx');

%% create csv with min, max, mean, median, q25, q75 and n of all substrates:
% XY: hier muss eine CSV entstehen, die den Substratnahmen in der ersten
% Zeile hat und dann untereinander die kondensierten statistischen Größen!
% Bisher ist das noch falsch und die CSV nicht brauchbar!
macNutrientsResults = nan(7,3*nSubstrates); % allocate memory
for kk = 1:nSubstrates
    substrateNumber = kk;    % 1...7 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    currMacroNutrientsTab = macroNutrientsTabs.(substrate); % get table of macronutrient values for one substrate
    currMacroNutrients = currMacroNutrientsTab{:,:};        % table -> array   
    nSamples = size(currMacroNutrients,1);
    minVals = min(currMacroNutrients); 
    maxVals = max(currMacroNutrients); 
    meanVals = mean(currMacroNutrients);
    medianVals = median(currMacroNutrients); 
    q25 = quantile(currMacroNutrients,0.25); % 25 percent quantile
    q75 = quantile(currMacroNutrients,0.75); % 75 percent quantile

    % put all statistical measures into array: 
    macNutrientsResults(:,kk*3-2:kk*3) = [minVals;maxVals;meanVals;medianVals;q25;q75;ones(1,3)*nSamples];
end

% save statistical results as csv: 
csvFileName = 'macroNutrientsResults.xlsx'; 
% save(csvFileName, 'macNutrientsResults')



 