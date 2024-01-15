from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np
from pathlib import Path

plt.rcParams["text.usetex"]


@dataclass_json
@dataclass(kw_only=True)
class MPLProperties:
    linewidth: float
    linestyle: str
    color: str = None


@dataclass_json
@dataclass(kw_only=True)
class PlotVarProperty:
    mpl_properties: MPLProperties
    label: str = None


@dataclass
class Constraint:
    value: float
    ax_idx: int


@dataclass(kw_only=True)
class PostProcessing:
    result_directory: str
    scenario_name: str
    _default_mpl_properties = MPLProperties(color=None, linewidth=2.0, linestyle="-")
    _default_plot_var_properties = PlotVarProperty(
        mpl_properties=_default_mpl_properties
    )

    _cmap = plt.colormaps.get_cmap("viridis")

    _SUBSTRATE_COLORS = {
        "CORN_SILAGE": "yellow",
        "GRASS_SILAGE": "green",
        "CATTLE_MANURE": "brown",
        "SUGAR_BEET_SILAGE": "orange",
    }

    def __post_init__(self):
        self._data_simulation = do_mpc.data.load_results(
            str(
                Path(self.result_directory, f"{self.scenario_name}_mpc_results.pkl")
            )  # "./results/{self.scenario_name}_mpc_results.pkl"
        )
        self._data_simulator = self._data_simulation["simulator"]
        self._data_mpc = self._data_simulation["mpc"]
        with open(
            Path(self.result_directory, f"{self.scenario_name}_scenario_meta_data.pkl"),
            "rb",
        ) as fp:
            self._scenario_meta_data = pickle.load(fp)

        self._num_u = len(self._scenario_meta_data["sub_names"])
        try:
            self._num_dictated_subs = len(
                self._scenario_meta_data["disturbances"]["dictated_feeding"]
            )
        except TypeError:
            self._num_dictated_subs = 0

        self._graphics_mpc = do_mpc.graphics.Graphics(self._data_mpc)

        self.MEASUREMENT_COLOR = self._cmap.resampled(1).colors[0]

    def plot(
        self,
        subplot_labels_and_vars: list[tuple[str, dict[str, PlotVarProperty]]],
        plot_inputs: bool = True,
        time_range: tuple[float] = None,
        plot_save_name: str = None,
        constraints: list[Constraint] = None,
        dpi: int = 600,
        show_plot: bool = True,
    ):
        if plot_inputs:
            subplot_labels_and_vars.insert(
                0, (r"$u_{feed}$" + "\n" + r"$[m^3/d]$", {f"u_norm": None})
            )
        fig, axes = plt.subplots(len(subplot_labels_and_vars), sharex=True)

        axes[-1].set_xlabel("Time [d]")

        for ax_idx, (ax, subplot_label_and_vars) in enumerate(
            zip(axes, subplot_labels_and_vars)
        ):
            y_label, plot_var_properties = subplot_label_and_vars

            labels = []
            constraints_drawn = False
            for plot_var_name, plot_var_property in plot_var_properties.items():
                try:
                    if plot_var_property.label is not None:
                        if not plot_var_property.label == "":
                            labels.append(plot_var_property.label)
                        label_set = True
                    else:
                        raise AttributeError
                except AttributeError:
                    label_set = False

                plt_kwargs = (
                    plot_var_property.mpl_properties.to_dict()
                    if isinstance(plot_var_property, PlotVarProperty)
                    else self._default_mpl_properties.to_dict()
                )
                if plot_var_name[0] == "x":
                    self._graphics_mpc.add_line(
                        var_type=f"_x",
                        var_name=plot_var_name,
                        axis=ax,
                        **plt_kwargs,
                    )

                    x_num = int(plot_var_name.split("_")[-1])
                    if not label_set:
                        labels.append(
                            self._scenario_meta_data["_state_names"][x_num - 1]
                        )
                elif plot_var_name[0] == "y":
                    y_num = int(plot_var_name.split("_")[-1])

                    if plt_kwargs["color"] is None:
                        plt_kwargs["color"] = self.MEASUREMENT_COLOR

                    ax.plot(
                        self._data_simulator._time,
                        self._data_simulator._y[
                            :, self._num_u + y_num + self._num_dictated_subs - 1
                        ],
                        **plt_kwargs,
                    )
                    if not label_set:
                        labels.append(
                            self._scenario_meta_data["_meas_names"][y_num - 1]
                        )
                elif plot_var_name[0] == "u":
                    colors = [
                        self._SUBSTRATE_COLORS[sub_name]
                        for sub_name in self._scenario_meta_data["sub_names"]
                    ]

                    for i in range(self._num_u):
                        plt_kwargs["color"] = colors[i]
                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._u[:, i]
                            * self._scenario_meta_data["Tu"][i],
                            **plt_kwargs,
                        )

                    # Add constraints to plot
                    ax.hlines(
                        min(self._scenario_meta_data["Tu"]),
                        xmin=0,
                        xmax=self._scenario_meta_data["n_days_mpc"],
                        color="black",
                        linestyle="--",
                    )

                    ax.hlines(
                        max(self._scenario_meta_data["Tu"]),
                        xmin=0,
                        xmax=self._scenario_meta_data["n_days_mpc"],
                        color="blue",
                        linestyle="--",
                    )

                    labels = [
                        sub.lower().replace("_", " ")
                        for sub in self._scenario_meta_data["sub_names"]
                    ] + [r"$u_{max,solid}$", r"$u_{max,liquid}$"]

                elif plot_var_name.startswith("dictated_sub_feed"):
                    feed_num = int(plot_var_name.split("_")[-1])

                    plt_kwargs["drawstyle"] = "steps-post"

                    ax.plot(
                        self._data_simulator._time,
                        self._data_simulator._tvp[:, feed_num - 1]
                        * self._scenario_meta_data["Tu"][self._num_u + feed_num - 1],
                        **plt_kwargs,
                    )
                    if not label_set:
                        labels.append(f"Dictated substrate num. {feed_num}")
                elif plot_var_name[0] != "u" and plot_var_name[0] != "x":
                    if plot_var_name in self._scenario_meta_data["aux_var_names"]:
                        aux_expression_idx = np.where(
                            np.array(self._scenario_meta_data["aux_var_names"])
                            == plot_var_name
                        )[0][0]
                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._aux[:, aux_expression_idx],
                            **plt_kwargs,
                        )
                        if not label_set:
                            labels.append(plot_var_name)
                    else:
                        self._graphics_mpc.add_line(
                            var_type=f"_tvp",
                            var_name=plot_var_name,
                            axis=ax,
                            **plt_kwargs,
                        )
                        if not label_set:
                            labels.append(plot_var_name)

            if constraints is not None:
                for constraint in constraints:
                    if constraint.ax_idx == ax_idx:
                        ax.hlines(
                            constraint.value,
                            0.0,
                            self._scenario_meta_data["n_days_mpc"],
                            color="red",
                            linestyle="--",
                        )
                        constraints_drawn = True

            if constraints_drawn:
                labels.append(r"$constraints$")
            if labels:
                ax.legend(labels=labels)
            ax.yaxis.set_label_coords(-0.1, 0)
            ax.set_ylabel(y_label, rotation=0)
            ax.grid(True, linestyle="--")

        time_start = 0.0
        if time_range is not None:
            plt.setp(axes, xlim=time_range)
            time_start = time_range[0]

        axes[0].annotate(
            text="",
            xy=(
                self._scenario_meta_data["controller_params"]["mpc_n_horizon"]
                * self._scenario_meta_data["t_step"],
                1,
            ),
            xytext=(time_start, 1),
            arrowprops=dict(arrowstyle="->"),
        )

        fig.set_size_inches(w=10, h=2.5 * len(axes))

        plt.tight_layout()
        if plot_save_name is not None:
            plt.savefig(
                fname=str(
                    Path(self.result_directory, "plots", f"{plot_save_name}.png")
                ),
                dpi=dpi,
                format="png",
            )
        if show_plot:
            plt.show()


if __name__ == "__main__":
    default_mpl_properties = MPLProperties(color="red", linewidth=3.0, linestyle="-")
    default_plot_property = PlotVarProperty(
        mpl_properties=default_mpl_properties, label=None
    )
    post_processing = PostProcessing(
        result_directory="./results", scenario_name="Scenario_2a_test"
    )
    post_processing.plot(
        [
            ("V'g", {"y_1": default_plot_property}),
            (
                "Volume flow ch4",
                {
                    "v_ch4_dot_tank_in": None,
                },
            ),
            ("pH", {"y_4": default_plot_property}),
        ]
    )
