from matplotlib import colors
import pandas as pd
from pathlib import Path
from ad_meal_prep_control.utils import CHP
#from ad_meal_prep_control.utils import Time
#from ad_meal_prep_control.utils import TimeSlot

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

@dataclass(order=True)
class Time:
    hour: int
    minute: int = 0

@dataclass
class TimeSlot:
    start_time: Time
    end_time: Time

# Schedule data
schedule = {
    0: (TimeSlot(Time(7), Time(15)), TimeSlot(Time(16), Time(22))),
    1: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(15), Time(22))),
    2: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(16), Time(22))),
    3: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(15), Time(22))),
    4: (TimeSlot(Time(7), Time(14)), TimeSlot(Time(16), Time(23))),
    5: (TimeSlot(Time(9), Time(12)), TimeSlot(Time(17), Time(23))),
    6: (
        TimeSlot(Time(0), Time(1)),
        TimeSlot(Time(11), Time(12)),
        TimeSlot(Time(17), Time(23, 59)),
    ),
}

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