%% Version
% (R2022b) Update 5
% Erstelldatum: 23.08.2023
% Autor: Simon Hellmann

% compute inlet concentrations frm lab measurements of substrates
clear all
clc

% specify substrate type: 
agriculturalSubstrates = {'Maissilage', 'Grassilage', 'Getreidestroh',...
                          'Rinderguelle', 'Schweineguelle', 'Huehertrockenkot'}; 
% BMP values from KTBL Gasausbeute in landwirtsch. Biogasanlagen (2021): 
bmpAgriculturalSubstrates = [357, 315, 240, 230, 230, 272]; % [L_N/kg_oTS]
substrateNumber = 1; % pick 1...6 for respective substrate
substrate = agriculturalSubstrates{substrateNumber}; 
BMP = bmpAgriculturalSubstrates(substrateNumber); 
BMP_stoich = 420;   % max. BMP for agricultural substrates [L_N/kg_FoTS]

% pull table from absolute path, respect the appropriate sheet: 
addpath('\\dbfz-user.leipzig.dbfz.de\user$\shellmann\Notizen & Unterlagen\Messdaten\Labordaten 2022 final') 
T = readtable('Labordaten_landwirtschaftliche_Substrate.xlsx','Sheet',substrate); 
colNamesOrig = T.Properties.VariableNames;      % get column names 
colNamesEdit = {'Probeneingangsnummer','TS','water','NH4N','XA','XP','XL','XC','Ac'}; 
colNamesMeasurements = colNamesEdit(2:end); 
labDataRaw = T{:,2:end}; % transform table into array without Probeneingangsnummer

% remove all rows which are not complete (having NaNs)
labDataEffPre = labDataRaw(~any(isnan(labDataRaw),2),:); % effective lab data
% remove all rows containing negative values: 
labDataEff = labDataEffPre(all(labDataEffPre > 0,2),:); 

% pick rows from labDataEff 1 by 1 and compute xIn, save in xInMat:
nSamples = size(labDataEff,1);  % # complete samples per substrate
xInMat = nan(18,nSamples); 
for k = 1:nSamples
    labData = labDataEff(k,:); 
    xInMat(:,k) = computeX_inR3FracFromLabMeasurements(labData,BMP,BMP_stoich); 
end

% compute mean of all xIn_i from the for loop as final result of substrate:
xInSubstrate = mean(xInMat,2); 

% create struct with available substrates as field names, save xInSubstrate
% therein and save the entire struct:
xInAgrSubstrates = struct; 
xInAgrSubstrates.(substrate) = xInSubstrate; 
fileName = 'xInAgrSubstrates.mat'; 
save(fileName, 'xInAgrSubstrates')
 