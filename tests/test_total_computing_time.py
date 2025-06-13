def read_computation_times(file_path):
    total_micro_seconds = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Assuming each line contains a time in seconds
            try:
                total_micro_seconds += float(line.strip())
            except ValueError:
                print(f"Skipping invalid line: {line.strip()}")

    # Convert total microseconds to hh:mm:ss
    total_seconds = total_micro_seconds / 1e6
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60

    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}", round(total_seconds)

my_scenario = '2a_nominal_1sigma'  # '1b' '2a_robust', '2a_nominal', '2c', '2c_uninhibited', '2a_nominal_1sigma', '2a_robust_1sigma'
if my_scenario == '1b':
    scen_name = 'Scenario_1b_quadratic_robust_feedback_mismatch_1.5std_3tsap'
elif my_scenario == '2a_nominal':
    scen_name = 'Scenario_2a_dynamic_nominal_feedback_mismatch_2std_8tsap'
elif my_scenario == '2a_nominal_1sigma':
    scen_name = 'Scenario_2a_dynamic_nominal_feedback_mismatch_1std_12tsap'
elif my_scenario == '2a_robust':
    scen_name = 'Scenario_2a_dynamic_robust_feedback_mismatch_2std_8tsap'
elif my_scenario == '2a_robust_1sigma':
    scen_name = 'Scenario_2a_dynamic_robust_feedback_mismatch_1std_12tsap'
elif my_scenario == '2c':
    scen_name = 'Scenario_2c_dynamic'
elif my_scenario == '2c_uninhibited':
    scen_name = 'Scenario_2c_uninhibited'

my_path = '../ad_meal_prep_control/scenarios/results/'
file_name = scen_name + '_mpc_computation_times_mikro_secs.txt'
run_time_hh_mm_ss, run_time_sec = read_computation_times(my_path + file_name)
print(f"\ncomputation time was: {run_time_hh_mm_ss} or {run_time_sec} seconds.")