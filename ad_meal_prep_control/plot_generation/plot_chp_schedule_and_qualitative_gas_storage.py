import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import ad_meal_prep_control.utils as utils
from ad_meal_prep_control.params_R3 import P_el_chp

# adjust default font sizes:
my_fs = 23
plt.rcParams.update({
    'font.size': 10,           # Default font size
    'axes.labelsize': my_fs,      # Axis label font size
    'axes.titlesize': 14,      # Title font size
    'xtick.labelsize': my_fs,     # X-axis tick label font size
    'ytick.labelsize': my_fs,     # Y-axis tick label font size
})

# get the fuckin Schedule data:
n_days = 30
t_step = 0.5 / 24.0
n_steps = round(n_days / t_step)
max_power = P_el_chp
return_schedule = True
schedule = utils.typical_ch4_vol_flow_rate(max_power, n_steps, t_step, return_schedule)

# Function to convert Time objects to hours in fraction
def time_to_hours(t):
    return t.hour + t.minute / 60

# Prepare data for plotting
days = list(range(7))
timeline = []
states = []

# Start in off condition
states.append(0)  # Initial off state
timeline.append(0)  # Start from hour 0 (beginning of the week)

# Iterate through each day
for day in days:
    for time_slot in schedule[day]:
        start_hour = time_to_hours(time_slot.start_time) + day * 24
        end_hour = time_to_hours(time_slot.end_time) + day * 24

        # Append On state
        timeline.append(start_hour)
        states.append(1)  # On state

        # Append Off state immediately after
        timeline.append(end_hour)
        states.append(0)  # Off state

# Create a higher resolution timeline (1 minute intervals)
high_res_timeline = np.arange(0, 24 * 60 * 7, 1)  # 1 minute resolution
gas_storage = []
current_level = 50  # Starting gas storage level at 50%

# Generate gas storage levels based on states
for minute in high_res_timeline:
    # Determine the current hour and minute
    current_hour = minute // 60
    current_minute = minute % 60
    current_time = current_hour + current_minute / 60

    # Find the current state based on the timeline
    current_state = 0  # Default to off
    for i in range(len(timeline) - 1):
        if timeline[i] <= current_time < timeline[i + 1]:
            current_state = states[i]
            break

    if current_state == 1:  # CHP is ON
        current_level += (150 - current_level) / 1000  # Linear increase towards 80%
    else:  # CHP is OFF
        current_level -= (current_level - 20) / 200  # Linear decrease towards 30%

    # Ensure the gas storage level stays within bounds
    current_level = max(0, min(100, current_level))
    gas_storage.append(current_level)

# Create a figure for gas storage level
plt.figure(figsize=(10, 3))
plt.plot(high_res_timeline, gas_storage, color='black', label='Gas Storage Level', linewidth=2)
plt.xlabel('Time [d]')
plt.ylabel('GS level [%]')
plt.xlim(0, 24 * 60 * 7)  # Limit x-axis to one week (0 to 10080 minutes)
plt.ylim(0, 100)  # Assuming gas storage level is between 0 and 100%
plt.xticks(np.arange(0, 24 * 60 * 7 + 1, 1440), labels=[day_k for day_k in range(8)])  # Label days

# Save the gas storage plot
plt.savefig(fname=str(
    Path(
        "../../results",
        "plots",
        f"gas_storage_level.png",
    )
), dpi=600, bbox_inches='tight')

plt.show()