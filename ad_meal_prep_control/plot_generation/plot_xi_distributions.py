import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from pathlib import Path
from matplotlib import colors

means_all = {
    "ch": [239.75428353, 161.63261293, 18.46756344, 443.44779286, 18.46756344],
    "pr": [26.3338422, 42.2833745, 13.3134649, 9.57299184, 13.3134649],
    "li": [7.99171082, 7.63325168, 2.00604756, 0.60788009, 2.00604756],
}
std_devs_all = {
    "ch": [25.35351729, 14.15546532, 3.10502837, 18.33059836, 7.76257093],
    "pr": [1.48845018, 2.38995494, 0.75250809, 0.54108783, 1.88127022],
    "li": [1.04091797, 0.99422877, 0.2612871, 0.0791762, 0.65321775],
}
substrate_names = [
    "CORN_SILAGE",
    "GRASS_SILAGE",
    "CATTLE_MANURE",
    "SUGAR_BEET_SILAGE",
    "CATTLE_MANURE_\nLARGE_UNCERTAINTY",
]

substrate_colors = {
    "CORN_SILAGE": "yellow",
    "GRASS_SILAGE": "green",
    "CATTLE_MANURE": "brown",
    "SUGAR_BEET_SILAGE": "orange",
    "CATTLE_MANURE_\nLARGE_UNCERTAINTY": "black",
}

fig, axes = plt.subplots(nrows=len(means_all), ncols=1)

for idx, (means, std_devs, ax) in enumerate(
    zip(means_all.values(), std_devs_all.values(), axes)
):
    ax.grid(True, linestyle="--")
    for mean, std_dev, sub_name in zip(means, std_devs, substrate_names):
        values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, 100)

        distribution = norm.pdf(values, mean, std_dev)

        if idx == 2:
            sns.lineplot(
                x=values,
                y=distribution,
                ax=ax,
                label=sub_name.lower().replace("_", " "),
                color=substrate_colors[sub_name],
            )
        else:
            sns.lineplot(
                x=values, y=distribution, ax=ax, color=substrate_colors[sub_name]
            )

        for i in [-1, 1]:
            ax.vlines(
                x=mean + 2.0 * i * std_dev,
                ymin=0,
                ymax=max(distribution),
                color=colors.to_rgba(substrate_colors[sub_name]),
            )

            ax.vlines(
                x=mean + 1.5 * i * std_dev,
                ymin=0,
                ymax=max(distribution),
                color="blue",
            )

        ax.fill_between(
            values,
            distribution,
            # where=(values > mean - 2.0 * std_dev) & (values < mean + 2.0 * std_dev),
            color=colors.to_rgba(substrate_colors[sub_name]),
            alpha=0.4,
        )


# fig.suptitle("Normalverteilung")

for ax, nutrient_name in zip(axes, means_all.keys()):
    ax.set_xlabel(r"$\xi" + "_" + "{" + f"{nutrient_name}" + "}$ $[g/kg_{TS}]$")
    ax.set_ylabel("density")

axes[-1].legend(ncols=2)

fig.set_size_inches(w=8, h=5)
fig.tight_layout()
fig.savefig(
    fname=str(
        Path(
            "/home/julius/Masterarbeit/ad_meal_prep_control/results",
            "plots",
            f"xi_pdfs.png",
        )
    ),
    dpi=600,
    format="png",
)
plt.show()
