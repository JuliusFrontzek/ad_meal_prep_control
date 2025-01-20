import pickle

# Specify the path to your .pkl file

#file_path = "../ad_meal_prep_control/scenarios/results/Scenario_1a_quadratic_no_feedback_mismatch_1std_ch_scenario_meta_data.pkl"
file_path = "../ad_meal_prep_control/scenarios/results/"
#file_name = "Scenario_2c_dynamic_mpc_results.pkl"
file_name_results = "Scenario_1a_quadratic_no_feedback_mismatch_1std_ch_scenario_meta_data.pkl"
file_name_meta = "Scenario_2c_dynamic_scenario_meta_data.pkl"
path_name_results = file_path + file_name_results
path_name_meta = file_path + file_name_meta

# Load the .pkl file
with open(path_name_results, 'rb') as file:
    results_data = pickle.load(file)

# Now you can use the loaded data