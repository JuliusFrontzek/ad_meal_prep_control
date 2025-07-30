import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from pathlib import Path
from matplotlib import colors
from ad_meal_prep_control.utils import Scenario
from ad_meal_prep_control import substrates

### preliminaries:
substrate_names = [
    "MAIZE_SILAGE",
    "GRASS_SILAGE",
    "SUGAR_BEET_SILAGE",
    "CATTLE_MANURE",
]
n_substrates = len(substrate_names)

substrate_colors = {
    "MAIZE_SILAGE": "orange",
    "GRASS_SILAGE": "limegreen",
    "SUGAR_BEET_SILAGE": "deeppink",
    "CATTLE_MANURE": "saddlebrown",
}
substrate_colors_dark = {
    "MAIZE_SILAGE": "darkorange",
    "GRASS_SILAGE": "green",
    "SUGAR_BEET_SILAGE": "mediumvioletred",
    "CATTLE_MANURE": "maroon",
}

### get data
## theoretical data:
# load input concentrations from excel:
file_path = '../../data/data_out/v2023/nominal_inlet_concentrations.xlsx'
macro_nutrients_data = pd.read_excel(file_path)

# Extract macro nutrients from excel:
substrate_cols = ['Maissilage', 'Grassilage', 'Zuckerruebensilage', 'Rinderguelle']
carbs = macro_nutrients_data.loc[5, substrate_cols].values
proteins = macro_nutrients_data.loc[7, substrate_cols].values
lipids = macro_nutrients_data.loc[8, substrate_cols].values

means_all = {  # ignore cattle manure very uncertain
    "ch": carbs,
    "pr": proteins,
    "li": lipids,
}

std_devs_all = {
    "ch": [40.085336, 36.613181, 49.423171, 5.378375],
    "pr": [1.540789, 2.473993, 0.560114, 0.778969],
    "li": [0.817209, 0.780554, 0.062160, 0.205133],
}

## measured data:
# Read all sheets into a dictionary of DataFrames
excel_file = '../../data/data_out/v2023/macro_nutrients_concentrations.xlsx'
all_sheets = pd.read_excel(excel_file, sheet_name=None)

df_cornsilage = all_sheets['Maissilage']
df_grassilage = all_sheets['Grassilage']
df_sugarbeet = all_sheets['Zuckerruebensilage']
df_manure = all_sheets['Rinderguelle']
dfs_substrates = [df_cornsilage, df_grassilage, df_sugarbeet, df_manure]

### draw plot:
fig, axes = plt.subplots(nrows=len(means_all), ncols=1)
right_axes = []  # Store the right y-axes for later iteration
offset_y = 7  # position boxplots shifted vertically

# iterate through subplot windows (ch/pr/li):
for idx_k, (means, std_devs, ax) in enumerate(zip(means_all.values(), std_devs_all.values(), axes)):
    axR = ax.twinx()  # instantiate a second Axes that shares the same x-axis
    # iterate through substrates:
    for sub_k, mean, std_dev, sub_name in zip(range(n_substrates), means, std_devs, substrate_names):
        values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, 100)  # draw +/- 3 sigmas

        distribution = norm.pdf(values, mean, std_dev)

        # add legend only to first plot
        if idx_k == 0:
            sns.lineplot(
                x=values,
                y=distribution,
                ax=ax,
                # make lowercase, replace "_", add sample size:
                label=[sub_name.lower().replace("_", " ").capitalize() + fr" ($n={dfs_substrates[sub_k].shape[0]}$)"],
                color=substrate_colors[sub_name],
            )
        else:
            sns.lineplot(
                x=values, y=distribution, ax=ax, color=substrate_colors[sub_name]
            )

        # add vertical lines at 2 sigmas for each substrate:
        # for i in [-1, 1]:
        #    ax.vlines(
        #        x=mean + 2.0 * i * std_dev,
        #        ymin=0,
        #        ymax=max(distribution),
        #        color=colors.to_rgba(substrate_colors[sub_name]),
        #    )

        # add vertical blue lines at 1.5 sigmas:
        # ax.vlines(
        #    x=mean + 1.5 * i * std_dev,
        #    ymin=0,
        #    ymax=max(distribution),
        #    color="blue",
        # )

        # add transparent color fill into normal plots:
        ax.fill_between(
            values,
            distribution,
            # where=(values > mean - 2.0 * std_dev) & (values < mean + 2.0 * std_dev),
            color=colors.to_rgba(substrate_colors[sub_name]),
            alpha=0.4,
        )

        # add boxplots:
        axR.boxplot(dfs_substrates[sub_k].iloc[:, idx_k],
                    vert=False,  # horizontal boxplots
                    patch_artist=True,
                    boxprops=dict(facecolor=substrate_colors[sub_name], color=substrate_colors_dark[sub_name],
                                  alpha=0.4),  # Box fill color with transparency
                    widths=0.5,
                    whiskerprops=dict(color=substrate_colors_dark[sub_name]),  # Whisker color
                    capprops=dict(color=substrate_colors[sub_name]),  # Cap color
                    medianprops=dict(color=substrate_colors_dark[sub_name]),  # Median line color
                    positions=[offset_y + sub_k],  # vertical positioning
                    showfliers=False,  # ignore outliers
                    # flierprops=dict(marker='o', color='pink', alpha=0.5)  # outlier color
                    )
    right_axes.append(axR)  # Store the right y-axis (of complete subplot) for later use

# fig.suptitle("Normalverteilung")
subtitles = ['Carbohydrates', 'Proteins', 'Lipids']
for k, ax, axR, nutrient_name in zip(range(len(axes)), axes, right_axes, means_all.keys()):
    ax.set_xlabel(
        r"$\xi" + "_" + "\mathrm{" + f"{nutrient_name}" + "}$ $[\mathrm{g} \; \mathrm{L}^{-1}]$"
    )
    ax.set_ylabel("Density")
    ax.set_title(subtitles[k])
    ax.grid(True, linestyle="--")
    # axR.grid(True, linestyle="--", axis="y")
    # adjust ylims of both axes:
    ylim_orig = ax.get_ylim()
    ax.set_ylim(ylim_orig[0], ylim_orig[1] * 1.5)
    axR.set_ylim([0, offset_y + 4])
    axR.set_axis_off()  # mute right ylabels
    # axR.set_yticks(list(range(7, 7+4)), substrate_names)  # Set y-ticks for groups

axes[0].legend(ncols=2, loc="lower right", bbox_to_anchor=(1, 0.2))  # make legend 2-column

fig.set_size_inches(w=8, h=8)
fig.tight_layout()
fig.savefig(
    fname=str(
        Path(
            "../../results",
            "plots",
            f"xi_pdfs_and_boxplots.png",
        )
    ),
    dpi=1000,
    format="png",
)
# __SH: additionally save as SVG:
fig.savefig(
    fname=str(
        Path(
            "../../results",
            "plots",
            f"xi_pdfs_and_boxplots.svg",
        )
    ),
    format="svg",
)
plt.show()
