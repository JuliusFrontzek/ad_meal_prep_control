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
    "CORN_SILAGE",
    "GRASS_SILAGE",
    "SUGAR_BEET_SILAGE",
    "CATTLE_MANURE",
]

substrate_colors = {
    "CORN_SILAGE": "orange",
    "GRASS_SILAGE": "limegreen",
    "SUGAR_BEET_SILAGE": "deeppink",
    "CATTLE_MANURE": "saddlebrown",
}
substrate_colors_dark = {
    "CORN_SILAGE": "darkorange",
    "GRASS_SILAGE": "green",
    "SUGAR_BEET_SILAGE": "maroon",
    "CATTLE_MANURE": "mediumvioletred",
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

means_all = {   # ignore cattle manure very uncertain
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

#sub_names = ['corn silage', 'grass silage', 'sugar beet silage', 'cattle manure']
#sub_colors = ['orange', 'limegreen', 'deeppink', 'saddlebrown']

df_carbs = pd.concat(
    [df_cornsilage.iloc[:, 0], df_grassilage.iloc[:, 0], df_sugarbeet.iloc[:, 0], df_manure.iloc[:, 0]], axis=1)
df_protein = pd.concat(
    [df_cornsilage.iloc[:, 1], df_grassilage.iloc[:, 1], df_sugarbeet.iloc[:, 1], df_manure.iloc[:, 1]], axis=1)
df_lipids = pd.concat(
    [df_cornsilage.iloc[:, 2], df_grassilage.iloc[:, 2], df_sugarbeet.iloc[:, 2], df_manure.iloc[:, 2]], axis=1)
# rename column headings:
dfs = [df_carbs, df_protein, df_lipids]
for k in range(3):
    dfs[k].columns = [substrate_names[k].lower() for k in range(4)]

### draw plot:
fig, axes = plt.subplots(nrows=len(means_all), ncols=1)

for idx, (means, std_devs, ax) in enumerate(zip(means_all.values(), std_devs_all.values(), axes)):
    ax.grid(True, linestyle="--")
    for mean, std_dev, sub_name in zip(means, std_devs, substrate_names):
        values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, 100)  # draw +/- 3 sigmas

        distribution = norm.pdf(values, mean, std_dev)

        # add legend only to first plot
        if idx == 0:
            sns.lineplot(
                x=values,
                y=distribution,
                ax=ax,
                label=sub_name.lower().replace("_", " "),  # make lowercase and replace "_"
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

        ax.fill_between(
            values,
            distribution,
            # where=(values > mean - 2.0 * std_dev) & (values < mean + 2.0 * std_dev),
            color=colors.to_rgba(substrate_colors[sub_name]),
            alpha=0.4,
        )

# fig.suptitle("Normalverteilung")
subtitles = ['carbohydrates', 'proteins', 'lipids']
for ax, nutrient_name, k in zip(axes, means_all.keys(), range(len(axes))):
    ax.set_xlabel(
        r"$\xi" + "_" + "{" + f"{nutrient_name}" + "}$ $[\mathrm{g} \mathrm{L}^{-1}]$"
    )
    ax.set_ylabel("density")
    ax.set_title(subtitles[k])

axes[0].legend(ncols=2)  # make legend 2-column

fig.set_size_inches(w=8, h=5)
fig.tight_layout()
fig.savefig(
    fname=str(
        Path(
            "../../results",
            "plots",
            f"xi_pdfs.png",
        )
    ),
    dpi=600,
    format="png",
)
plt.show()
