%% Version: Matlab R2024a, Update 3
% Author: Simon Hellmann (thanks youChat)
% first created: 2025-01-03

clc
clear
close all

%% user settings: 
model_structure = 'R3-frac'; 
flag_delete_rows_with_nan = 1; % 1: do, 0: dont
flag_CH_method = 1; % 1: compute CH via BMP, 2: via lignin

%% load data of substrate characterization:
path_name = '../data/data_in/'; 
file_name = 'substrate_characterization_dbfz.xlsx'; 
pathNFile = [path_name,file_name]; 
substrate_characterization_dbfz = read_exceltabs_to_struct(pathNFile);

% create a table with full ADM1-R3 input for all substrates
substrates = {'maize_silage','grass_silage','sugarbeet_silage','cattle_manure'}; 
n_substrates = numel(substrates); 
n_states = 18; % ADM1-R3-frac
adm1_inputs = nan(n_states,n_substrates); 
n_samples = nan(1,n_substrates); % # samples used for computation

for k = 1:n_substrates
    % get substrate data and compute mean: 
    my_sub = substrates{k}; 
    my_sub_table = substrate_characterization_dbfz.(my_sub); 
    
    % optionally delete incomplete rows (except pH): 
    if flag_delete_rows_with_nan == 1 
        idx_nan_critical_cols = [3:10,13:20]; 
        idx_complete_rows = ~any(isnan(my_sub_table{:,idx_nan_critical_cols}),2); 
        my_sub_table = my_sub_table(idx_complete_rows,:); 
    end
    
    idx_numeric_cols = [3:11,13:20];
    my_sub_table_num = my_sub_table(:,idx_numeric_cols); % get only numeric columns
    my_sub_table_means = mean(my_sub_table_num,'omitmissing'); % get column averages (ignore NaNs)

    % compute ADM1-R3 input vector: 
    adm1_inputs(:,k) = compute_adm1_input(my_sub_table_means,my_sub, ...
            model_structure,flag_CH_method); 

    n_samples(k) = size(my_sub_table_num,1); 
end

% turn arrays into table: 
adm1_inputs_tab = array2table(adm1_inputs, 'VariableNames', substrates); 
n_samples_tab = array2table(n_samples, 'VariableNames', substrates); 

% save: 
writetable(adm1_inputs_tab,'../data/data_out/adm1_inputs.xlsx');
writetable(n_samples_tab,'../data/data_out/n_samples.xlsx');

%% helper functions (thanks youChat) 
function data_struct = read_exceltabs_to_struct(filename)
    % Read the names of the sheets in the Excel file
    sheet_names = sheetnames(filename);
    
    % Initialize an empty struct
    data_struct = struct();
    
    % Loop through each sheet and read the table
    for i = 1:length(sheet_names)
        % Read the table from the current sheet
        data = readtable(filename, 'Sheet',sheet_names{i}, 'VariableNamingRule','preserve');
        
        % simplify column headings: 
        simple_col_names = {'Probeneingangsnummer','Beschreibung','TS','oTS','NH4_N','XA','XP','XL','XLig','TKN','pH','GC-Probenbezeichnung','dilution','AC','PRO','i_BU','n_BU','i_VA','n_VA','AC_eff','Kommentar'};
        data.Properties.VariableNames = simple_col_names; 
        
        % Assign the table to the struct with the sheet name as the field name
        data_struct.(sheet_names{i}) = data;
    end
end