import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
import pandas as pd
from pathlib import Path
from webencodings import labels

##### get the fuckin data from excel, motherfucker!
# draw a single boxplot only for carbs of corn silage
# Read all sheets into a dictionary of DataFrames
excel_file = '../../data/data_out/v2023/macro_nutrients_concentrations.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)

df_cornsilage = all_sheets['Maissilage']
df_grassilage = all_sheets['Grassilage']
df_sugarbeet = all_sheets['Zuckerruebensilage']
df_manure = all_sheets['Rinderguelle']
dfs_substrates = [df_cornsilage, df_grassilage, df_sugarbeet, df_manure]

sub_names = ['corn silage', 'grass silage', 'sugar beet silage', 'cattle manure']
sub_colors = ['orange', 'limegreen', 'deeppink', 'saddlebrown']
sub_colors_dark = ['darkorange', 'green', 'mediumvioletred', 'maroon']

df_carbs = pd.concat(
    [df_cornsilage.iloc[:, 0], df_grassilage.iloc[:, 0], df_sugarbeet.iloc[:, 0], df_manure.iloc[:, 0]], axis=1)
df_protein = pd.concat(
    [df_cornsilage.iloc[:, 1], df_grassilage.iloc[:, 1], df_sugarbeet.iloc[:, 1], df_manure.iloc[:, 1]], axis=1)
df_lipids = pd.concat(
    [df_cornsilage.iloc[:, 2], df_grassilage.iloc[:, 2], df_sugarbeet.iloc[:, 2], df_manure.iloc[:, 2]], axis=1)
dfs = [df_carbs, df_protein, df_lipids]
for k in range(3):
    dfs[k].columns = sub_names

# Create the figure and axes
fig, axes = plt.subplots(nrows=3, figsize=(4, 5))
right_axes = []  # Store the right y-axes for later iteration

# iterate over macro nutrients:
for k, ax in enumerate(axes):
    axR = ax.twinx()  # instantiate a second Axes that shares the same x-axis
    # iterate over substrates:
    for subk in range(len(dfs_substrates)):
        axR.boxplot(dfs_substrates[subk].iloc[:, k],
                    vert=False,  # horizontal boxplots
                    patch_artist=True,
                    boxprops=dict(facecolor=sub_colors[subk], color=sub_colors_dark[subk], alpha=0.4),
                    # Box fill color with transparency
                    whiskerprops=dict(color=sub_colors_dark[subk]),  # Whisker color
                    capprops=dict(color=sub_colors[subk]),  # Cap color
                    medianprops=dict(color=sub_colors_dark[subk]),  # Median line color
                    positions=[subk],
                    showfliers=False,  # ignore outliers
                    # flierprops=dict(marker='o', color='pink', alpha=0.5)  # outlier color
                    )
    right_axes.append(axR)  # Store the right y-axis (of complete subplot) for later use

subtitles = ['carbohydrates', 'proteins', 'lipids']
xlabels = [r"$\xi_{ch}$", r"$\xi_{pr}$", r"$\xi_{li}$"]

# Iterate through all subplots and adjust labels
for k, ax, axR in zip(range(3), axes, right_axes):
    ax.set_title(subtitles[k])
    ax.set_xlabel(xlabels[k])
    # Mute the left y-axis
    ax.set_ylabel("")  # Remove the left y-axis label
    ax.tick_params(left=False, labelleft=False)  # Disable left ticks and labels
    # Add labels to the right y-axis:
    axR.set_yticks(list(range(4)), sub_names)  # Set y-ticks for groups

# Show the plot
plt.tight_layout()
fig.savefig(
    fname=str(
        Path(
            "../../results",
            "plots",
            f"xi_distributions.png",
        )
    ),
    dpi=600,
    format="png",
)
plt.show()
