%% Version
% (R2022b) Update 5
% Erstelldatum: 23.08.2023
% Autor: Simon Hellmann

% compute inlet concentrations from lab measurements of substrates
clear all
clc

% specify substrate type: 
agriculturalSubstrates = {'Maissilage', 'Grassilage', 'Getreidestroh',...
                          'Rinderguelle', 'Schweineguelle', 'HTK'}; 
% BMP values from KTBL Gasausbeute in landwirtsch. Biogasanlagen (2021): 
bmpAgriculturalSubstrates = [357, 315, 240, 230, 230, 272]; % [L_N/kg_oTS]

% compute inlet concentrations for all 6 agricultural substrates: 
nSubstrates = numel(bmpAgriculturalSubstrates); 
% create struct with available substrates as field names, save xInSubstrate
% therein and save the entire struct:
xInAgrSubstrates = struct; 
for k = 1:nSubstrates
    substrateNumber = k; % 1...6 for respective substrate
    substrate = agriculturalSubstrates{substrateNumber}; 
    BMP = bmpAgriculturalSubstrates(substrateNumber); 
    BMP_stoich = 420;   % max. BMP for agricultural substrates [L_N/kg_FoTS]
    
    % pull table from absolute path, respect the appropriate sheet: 
    addpath('\\dbfz-user.leipzig.dbfz.de\user$\shellmann\Notizen & Unterlagen\Messdaten\Labordaten 2021_2022 final') 
    T2021 = readtable('2021_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate); 
    T2022 = readtable('2022_Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate); 
    colNamesOrig = T2022.Properties.VariableNames;      % get column names 
    colNamesEdit = {'Probeneingangsnummer','TS','water','NH4N','XA','XP','XL','XC','AcMeas','dilutionFactor','Ac'}; 
    colNamesMeasurements = colNamesEdit(2:end); 
    % transform tables into arrays without Probeneingangsnummer, unify to 1 array:
    labDataRaw2021 = T2021{:,2:end}; 
    labDataRaw2022 = T2022{:,2:end};
    labDataRaw = [labDataRaw2021;labDataRaw2022]; 
    
    % remove all rows which are not complete (having NaNs)
    labDataEffPre = labDataRaw(~any(isnan(labDataRaw),2),:); % effective lab data
    % remove all rows containing negative values: 
    labDataEff = labDataEffPre(all(labDataEffPre > 0,2),:); 
    
    % pick rows from labDataEff 1 by 1 and compute xIn, save in xInMat:
    nSamples = size(labDataEff,1);  % # complete samples per substrate
    xInMat = nan(18,nSamples); % input vector size for ADM1-R3-frac
    for kk = 1:nSamples
        labData = labDataEff(kk,:); 
        xInMat(:,kk) = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich); 
    end % for
    
    % compute mean of all xIn_i from the for loop as final result of substrate:
    xInSubstrate = mean(xInMat,2); 
    xInAgrSubstrates.(substrate) = xInSubstrate; 
end % for

fileName = 'xInAgrSubstrates.mat'; 
save(fileName, 'xInAgrSubstrates')
 