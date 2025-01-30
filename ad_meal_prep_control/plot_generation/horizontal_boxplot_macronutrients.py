import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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

cornsilage_carbs = df_cornsilage['Xch']
sugarbeet_carbs = df_sugarbeet['Xch']

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
fig, axes = plt.subplots(nrows=3, figsize=(5, 4))

# iterate over macro nutrients:
for k, ax in enumerate(axes):
    # iterate over substrates: 
    for subk in range(len(dfs_substrates)):
        ax.boxplot(dfs_substrates[subk].iloc[:, k],
                   vert=False,  # horizontal boxplots
                   patch_artist=True,
                   boxprops=dict(facecolor=sub_colors[subk], color=sub_colors[subk], alpha=0.5),
                   # Box fill color with transparency
                   whiskerprops=dict(color=sub_colors[subk]),  # Whisker color
                   capprops=dict(color=sub_colors[subk]),  # Cap color
                   medianprops=dict(color=sub_colors[subk]),  # Median line color
                   showfliers=False,  # ignore outliers
                   positions=[subk]
                   # flierprops=dict(marker='o', color='pink', alpha=0.5)  # outlier color
                   )

subtitles = ['carbohydrates', 'proteins', 'lipids']
xlabels = [r"$\xi_{ch}$", r"$\xi_{pr}$", r"$\xi_{li}$"]
ylabels = ['mais', 'ruebe']
for k, ax in enumerate(axes):
    ax.set_title(subtitles[k])
    ax.set_xlabel(xlabels[k])
    ax.set_yticks(list(range(4)), sub_names)  # Set y-ticks for groups
    ax.grid(True, linestyle="--")

# Show the plot
plt.tight_layout()
plt.show()

bla = 1
# Customize boxplots
# for idx, patch in enumerate(box['boxes']):
#    #patch.set_facecolor('orange')  # Set box color: specific of each substrate
#    patch.set_edgecolor('orange')
#    patch.set_wiskers('orange')
#    patch.set_alpha(0.3)  # Set transparency (alpha)

# Add labels for each boxplot
# ax.set_yticklabels(['substrate'])  # ggf weglassen XY

# # Add title and labels
# ax.set_title(['carbohydrates', 'proteins', 'lipids'])
# ax.set_xlabel(r"$\xi_{ch} [\mathrm{g} \mathrm{L}^{-1}]$")
# ax.set_ylabel('corn silage')
