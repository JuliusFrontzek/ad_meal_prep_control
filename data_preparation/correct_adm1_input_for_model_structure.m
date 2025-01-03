%% Version: Matlab R2022b, Update 6
% Author: Simon Hellmann
% first created: 2024-12-14

function xIn_after = correct_adm1_input_for_model_structure(xIn_before,model_structure)
% Info: adapt vector of inlet concentrations to model structures other 
% than ADM1-R3

% Outputs: 
% - xIn_after:  vector of inlet concentrations after adaptation

% Inputs: 
% - xIn_after:  vector of inlet concentrations before adaptation (based on
%               ADM1-R3)
% - model_structure (string): (string) 'R3' (with ash, water and ion), 
%   'R3-Core': (no water/ash/ions), 'R3-frac': (2nd CH fraction), 'R4'

switch model_structure
    case "R3"
        xIn_after = xIn_before;     % R3: no change necessary
    case "R3-frac"
        xIn_after = nan(numel(xIn_before)+1,1);   % R3-frac: add second CH fraction
        xIn_after(1:6) = xIn_before(1:6); 
        xIn_after(8:end) = xIn_before(7:end); 
        xIn_after(7) = 0; % all CH to fast fraction
    case "R3-Core"
        xIn_after = xIn_before([1:4,6:10,13:17]); % R3-Core: omit h2o, ash and S_ion
    case "R4"
        xIn_after = nan(11,1);  % allocate memory
        xIn_after([1:4,6,7,9:11]) = xIn_before([2:5,7,8,11,16,17]); % all ecxept X_ch & X_bac
        xIn_after(5) = xIn_before(1) + xIn_before(6);   % add acetic acid to carbohydrates
        xIn_after(8) = xIn_before(9) + xIn_before(10);  % combine X_bac and X_ac
end