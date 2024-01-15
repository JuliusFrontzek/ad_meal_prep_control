from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np
from pathlib import Path
import sys

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = True

plt.rcParams["text.usetex"]
cmap = plt.colormaps.get_cmap("viridis")
result_directory = "/home/julius/Projects/ad_meal_prep_control/results"
fig, ax = plt.subplots(3, 1, gridspec_kw={"height_ratios": [3, 2, 1]})

num_u = 4
y_num = 4
num_dictated_subs = 1

control_horizons = [5, 10, 15]
robust_horizons = [0, 1]

num_scenarios = len(control_horizons) * len(robust_horizons)

colors = cmap.resampled(num_scenarios).colors

# x1, x2, y1, y2 = 4.0, 6.0, 500.0, 700.0
# axins = ax[0].inset_axes(
#     [0.5, 0.5, 0.47, 0.47], xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[]
# )

color_idx = 0
for n_horizon in control_horizons:
    for n_robust in robust_horizons:
        color = colors[color_idx]
        scenario_name = f"Scenario_1b_nc_{n_horizon}_nr_{n_robust}"

        data_simulation = do_mpc.data.load_results(
            str(
                Path(result_directory, f"{scenario_name}_mpc_results.pkl")
            )  # "./results/{self.scenario_name}_mpc_results.pkl"
        )
        data_simulator = data_simulation["simulator"]
        data_mpc = data_simulation["mpc"]
        with open(
            Path(result_directory, f"{scenario_name}_scenario_meta_data.pkl"),
            "rb",
        ) as fp:
            scenario_meta_data = pickle.load(fp)

        graphics_mpc = do_mpc.graphics.Graphics(data_mpc)

        aux_expression_idx = np.where(
            np.array(scenario_meta_data["aux_var_names"]) == "v_ch4_dot_tank_in"
        )[0][0]
        ax[0].plot(
            data_simulator._time,
            data_simulator._aux[:, aux_expression_idx],
            label=r"$N_c=$" + str(n_horizon) + r"$N_r=$" + str(n_robust),
            color=color,
        )
        # axins.plot(
        #     data_simulator._time,
        #     data_simulator._aux[:, aux_expression_idx],
        #     color=color,
        # )
        ax[1].plot(
            data_simulator._time,
            data_simulator._y[:, 4 + y_num + num_dictated_subs - 1],
            color=color,
        )

        ax[2].bar(
            color_idx,
            np.sum(
                data_simulator._u
                @ scenario_meta_data["Tu"][:-1]
                * scenario_meta_data["t_step"]
                / scenario_meta_data["n_days_mpc"]
            ),
            color=color,
            width=0.5,
        )

        color_idx += 1


graphics_mpc.add_line(
    var_type=f"_tvp",
    var_name="v_ch4_dot_tank_in_setpoint",
    axis=ax[0],
    linestyle="-.",
    color="black",
    label="Reference",
)

# graphics_mpc.add_line(
#     var_type=f"_tvp",
#     var_name="v_ch4_dot_tank_in_setpoint",
#     axis=axins,
#     linestyle="-.",
#     color="black",
#     label="Reference",
# )


# ax[0].indicate_inset_zoom(axins, edgecolor="black")

ax[0].set_ylim(300, 700)
ax[0].set_ylabel(r"$\dot V_{CH_4}$" + "\n" + r"$[m^3/d]$", rotation=0, labelpad=20)
ax[0].legend(ncol=2)
ax[1].set_ylabel(r"$pH$" + "\n" + r"$[1]$", rotation=0, labelpad=20)
ax[1].set_xlabel("Time [d]")
ax[2].set_ylabel(r"$u_{mean}$" + "\n" + r"$[m^3/d]$", rotation=0, labelpad=20)

for ax_ in ax:
    ax_.grid
    ax_.grid(True, linestyle="--")

fig.set_size_inches(w=10, h=7.5)
plt.tight_layout()
# plt.savefig(
#     fname=str(Path(result_directory, "plots", f"scenario_1b_horizon_comparison.png")),
#     dpi=dpi,
#     format="png",
# )

if show_plot:
    plt.show()
