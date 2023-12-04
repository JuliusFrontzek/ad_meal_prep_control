%% Version
% (R2022b) Update 5
% Erstelldatum: 23.08.2023
% last modified: 04.12.2023
% Autor: Simon Hellmann

% compute inlet concentrations from lab measurements of substrates
close all
clear
clc

% specify substrate type: 
agriculturalSubstrates = {'Maissilage', 'Grassilage', 'Getreidestroh', 'Zuckerruebensilage', ...
                          'Rinderguelle', 'Schweineguelle', 'HTK'}; 
% BMP values from KTBL Gasausbeute in landwirtsch. Biogasanlagen (2021): 
bmpAgriculturalSubstrates = [357, 315, 240, 349, 230, 230, 272]; % [L_N/kg_oTS]

% compute inlet concentrations for all 7 agricultural substrates: 
nSubstrates = numel(bmpAgriculturalSubstrates); 
% create struct with available substrates as field names, save xInSubstrate
% therein and save the entire struct:
xInAgrSubstrates = struct; 
TSAgrSubstrates = struct; 
for k = 1:nSubstrates
    substrateNumber = k; % 1...7 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    BMP = bmpAgriculturalSubstrates(substrateNumber); 
    BMP_stoich = 420;   % max. BMP for agricultural substrates [L_N/kg_FoTS]
    
    % pull all tables from absolute path, respect the appropriate sheet: 
    addpath('\\dbfz-user.leipzig.dbfz.de\user$\shellmann\Notizen & Unterlagen\Messdaten\Labordaten 2018_2022') 
    T2018 = readtable('2018_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2019 = readtable('2019_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2020 = readtable('2020_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2021 = readtable('2021_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    T2022 = readtable('2022_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate,'VariableNamingRule','preserve'); 
    colNamesOrig = T2022.Properties.VariableNames;      % get column names (the same for all tables)
    colNamesEdit = {'Probeneingangsnummer','TS','water','NH4N','XA','XP','XL','XC','AcMeas','dilutionFactor','Ac','Kommentar'}; 
%     colNamesMeasurements = colNamesEdit(2:end); 

    % transform tables into arrays without Probeneingangsnummer and Kommentar, unify to 1 array:
    labDataRaw2018 = T2018{:,2:11}; 
    labDataRaw2019 = T2019{:,2:11}; 
    labDataRaw2020 = T2020{:,2:11}; 
    labDataRaw2021 = T2021{:,2:11}; 
    labDataRaw2022 = T2022{:,2:11};
    labDataRaw = [labDataRaw2018;labDataRaw2019;labDataRaw2020;labDataRaw2021;labDataRaw2022]; 
    
    % remove all rows which are not complete (having NaNs)
    labDataEffPre = labDataRaw(~any(isnan(labDataRaw),2),:); % effective lab data
    % remove all rows containing negative values: 
    labDataEff = labDataEffPre(all(labDataEffPre > 0,2),:); 
    
    % pick rows from labDataEff 1 by 1 and compute xIn, save in xInMat:
    nSamples = size(labDataEff,1);  % # complete samples per substrate
    xInMat = nan(18,nSamples);  % input vector size for ADM1-R3-frac
    TSVec = nan(1,nSamples);    % 
    for kk = 1:nSamples
        labData = labDataEff(kk,:); 
        [xInMat(:,kk),TSVec(kk)] = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich); 
    end % for
    
    % compute mean of all xIn_i and TS from the for loop and save in
    % separate structs: 
    xInSubstrate = mean(xInMat,2);
    TSSubstrate = mean(TSVec);
    xInAgrSubstrates.(substrate) = xInSubstrate; 
    TSAgrSubstrates.(substrate) = TSSubstrate;  
end % for

fileName = 'xIn_TS_AgrSubstrates.mat'; 
save(fileName, 'xInAgrSubstrates', 'TSAgrSubstrates')
 