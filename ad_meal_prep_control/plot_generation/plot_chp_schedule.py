from matplotlib import colors
import pandas as pd
from pathlib import Path
import ad_meal_prep_control.utils as utils
from ad_meal_prep_control.utils import Time
from ad_meal_prep_control.utils import TimeSlot
from ad_meal_prep_control.params_R3 import P_el_chp

import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass

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

# Create the figure with specified size (width, height)
plt.figure(figsize=(10, 3))

# Create a stair graph showing ON and OFF states
plt.step(timeline, states, where='post', label='State', color='black')

plt.xlabel('Time [d]')
plt.ylabel('CHP status')

# Adjust x-axis to show days
hour_ticks = np.arange(0, 24 * 7 + 1, 24)
plt.xticks(hour_ticks, labels=[day_k for day_k in range(len(hour_ticks))])

plt.yticks([0, 1], ['Off', 'On'])
#plt.grid()

plt.xlim(0, 24 * 7)  # Limit x-axis to one week (0 to 168 hours)
plt.ylim(-0.05, 1.05)  # Adjust y-axis limits for clarity

# Save the plot as a high-resolution PNG file
plt.savefig(fname=str(
        Path(
            "../../results",
            "plots",
            f"chp_schedule.png",
        )
    ), dpi=600, bbox_inches='tight')

plt.show()