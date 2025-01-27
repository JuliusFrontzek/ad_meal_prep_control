from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np
from pathlib import Path
import math
from ad_meal_prep_control import substrates
from collections.abc import Iterable

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
    scaling: float = 1.0


@dataclass
class Constraint:
    value: float
    ax_idx: int
    color: str = None


@dataclass(kw_only=True)
class PostProcessing:
    result_directory: str
    scenario_name: str
    _default_mpl_properties = MPLProperties(color=None, linewidth=0.8, linestyle="-")
    _default_plot_var_properties = PlotVarProperty(
        mpl_properties=_default_mpl_properties
    )

    _cmap = plt.colormaps.get_cmap("viridis")

    _SUBSTRATE_COLORS = {
        "CORN_SILAGE": "orange",
        "GRASS_SILAGE": "limegreen",
        "CATTLE_MANURE": "sienna",
        "SUGAR_BEET_SILAGE": "deeppink",
        "SUGAR_BEET_SILAGE_VERY_UNCERTAIN": "saddlebrown",
    }

    def __post_init__(self):
        self._data_simulation = do_mpc.data.load_results(
            str(Path(self.result_directory, f"{self.scenario_name}_mpc_results.pkl"))
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
        height_ratios: list[float] = None,
        input_inset_axis: dict[str, tuple] = None,
        other_inset_axes: list[dict[str, tuple]] = None,
        color_background_indices: tuple[int] = None,
    ):
        if plot_inputs:
            subplot_labels_and_vars.insert(
                0, (r"$\dot V_{feed,silage}$" + "\n" + r"$[m^3/d]$", {f"u_norm": None})
            )

        if height_ratios is None:
            height_ratios = [1 for _ in range(len(subplot_labels_and_vars))]
        fig, axes = plt.subplots(
            len(subplot_labels_and_vars),
            sharex=True,
            gridspec_kw={"height_ratios": height_ratios},
        )

        if plot_inputs:
            ax_inputs_liquid = axes[0].twinx()

            if input_inset_axis is not None:
                # unpack elements "days" and "ylimit" of input_inset_axis and save:
                x1, x2, y1, y2 = (
                    *input_inset_axis["days"],
                    *input_inset_axis["ylimit"],
                )
                axins_input_feed = axes[0].inset_axes(
                    input_inset_axis["inset_axis_specs"],
                    xlim=(x1, x2),
                    ylim=(y1, y2),
                    # xticklabels=[],
                    # yticklabels=[],
                )
                axins_input_feed_liquid = axins_input_feed.twinx()

                axins_input_feed.grid(True, linestyle="--")

            inset_axes = {}
            if other_inset_axes is not None:
                for inset_ax in other_inset_axes:
                    x1, x2, y1, y2 = (
                        *inset_ax["days"],
                        *inset_ax["ylimit"],
                    )
                    _inset_ax = axes[inset_ax["plot_idx"]].inset_axes(
                        inset_ax["inset_axis_specs"],
                        xlim=(x1, x2),
                        ylim=(y1, y2),
                        # xticklabels=[],
                        # yticklabels=[],
                    )
                    _inset_ax.grid(True, linestyle="--")

                    if inset_ax["plot_idx"] in inset_axes.keys():
                        inset_axes[inset_ax["plot_idx"]].append(_inset_ax)
                    else:
                        inset_axes[inset_ax["plot_idx"]] = [_inset_ax]

        axes[-1].set_xlabel("Time [d]")

        for ax_idx, (axis, subplot_label_and_vars) in enumerate(
            zip(axes, subplot_labels_and_vars)
        ):
            axes_stacked = [axis]

            for ax in axes_stacked:
                y_label, plot_var_properties = subplot_label_and_vars

                for plot_var_name, plot_var_property in plot_var_properties.items():
                    plt_kwargs = (
                        plot_var_property.mpl_properties.to_dict()
                        if isinstance(plot_var_property, PlotVarProperty)
                        else self._default_mpl_properties.to_dict()
                    )
                    if plot_var_name[0] == "x":
                        x_num = int(plot_var_name.split("_")[-1])

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._x[:, x_num - 1]
                            * plot_var_property.scaling,
                            label=plot_var_property.label,
                            **plt_kwargs,
                        )

                        if ax_idx in inset_axes:
                            for inset_ax in inset_axes[ax_idx]:
                                inset_ax.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._x[:, x_num - 1]
                                    * plot_var_property.scaling,
                                    label=plot_var_property.label,
                                    **plt_kwargs,
                                )

                                ax.indicate_inset_zoom(
                                    inset_ax,
                                    edgecolor="black",
                                    linewidth=1.0,
                                )

                        x_num = int(plot_var_name.split("_")[-1])

                    elif plot_var_name[0] == "y":
                        y_num = int(plot_var_name.split("_")[-1])

                        if plt_kwargs["color"] is None:
                            plt_kwargs["color"] = self.MEASUREMENT_COLOR

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._y[
                                :, self._num_u + y_num + self._num_dictated_subs - 1
                            ],
                            label=plot_var_property.label,
                            **plt_kwargs,
                        )

                        if ax_idx in inset_axes:
                            for inset_ax in inset_axes[ax_idx]:
                                inset_ax.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._y[
                                        :,
                                        self._num_u
                                        + y_num
                                        + self._num_dictated_subs
                                        - 1,
                                    ],
                                    label=plot_var_property.label,
                                    **plt_kwargs,
                                )

                                ax.indicate_inset_zoom(
                                    inset_ax,
                                    edgecolor="black",
                                    linewidth=1.0,
                                )

                    elif plot_var_name[0] == "u":
                        colors = [
                            self._SUBSTRATE_COLORS[sub_name]
                            for sub_name in self._scenario_meta_data["sub_names"]
                        ]

                        for i in range(self._num_u):
                            plt_kwargs["color"] = colors[i]

                            sub_name = self._scenario_meta_data["sub_names"][i]
                            sub = getattr(substrates, sub_name)

                            if sub.state == "solid":
                                if input_inset_axis is not None:
                                    axins_input_feed.plot(
                                        self._data_simulator._time,
                                        self._data_simulator._u[:, i]
                                        * self._scenario_meta_data["Tu"][i],
                                        # * self._scenario_meta_data["u_max"]["solid"]
                                        # / self._scenario_meta_data["u_max"][sub.state],
                                        **plt_kwargs,
                                    )
                                ax.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._u[:, i]
                                    * self._scenario_meta_data["Tu"][i],
                                    label=sub_name.lower().replace("_", " "),
                                    **plt_kwargs,
                                )

                            else:
                                ax_inputs_liquid.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._u[:, i]
                                    * self._scenario_meta_data["Tu"][i],
                                    label=sub_name.lower().replace("_", " "),
                                    **plt_kwargs,
                                )
                                if input_inset_axis is not None:
                                    axins_input_feed_liquid.plot(
                                        self._data_simulator._time,
                                        self._data_simulator._u[:, i]
                                        * self._scenario_meta_data["Tu"][i],
                                        **plt_kwargs,
                                    )

                                    axins_input_feed_liquid.set_ylim(
                                        axins_input_feed.get_ylim()
                                    )

                            if input_inset_axis is not None:
                                ax.indicate_inset_zoom(
                                    axins_input_feed, edgecolor="black", linewidth=1.0
                                )

                        # Add constraints to plot
                        # ax.hlines(
                        #    self._scenario_meta_data["u_max"]["solid"],
                        #    xmin=0,
                        #    xmax=self._scenario_meta_data["n_days_mpc"],
                        #    color="black",
                        #    linestyle="--",
                        #    linewidth= 1
                        #    # label=r"$u_{max}$",
                        # )

                        ax_inputs_liquid.hlines(
                            self._scenario_meta_data["u_max"]["liquid"],
                            xmin=0,
                            xmax=self._scenario_meta_data["n_days_mpc"],
                            color="black",
                            linestyle=(0, (5, 5)),
                            linewidth=1,
                            # label=r"$u_{max}$",
                        )

                        if "Scenario_1" in plot_save_name:
                            ax.set_ylim(
                                0.0, self._scenario_meta_data["u_max"]["solid"] * 0.1
                            )
                            ax_inputs_liquid.set_ylim(
                                0.0, self._scenario_meta_data["u_max"]["liquid"] * 0.05
                            )
                        elif "Scenario_2a" or "Scenario_2c" in plot_save_name:
                            ax.set_ylim(
                                0.0, self._scenario_meta_data["u_max"]["solid"] * 1.1
                            )
                            ax_inputs_liquid.set_ylim(
                                0.0, self._scenario_meta_data["u_max"]["liquid"] * 1.1
                            )

                    elif plot_var_name.startswith("dictated_sub_feed"):
                        feed_num = int(plot_var_name.split("_")[-1])

                        plt_kwargs["drawstyle"] = "steps-post"

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._tvp[:, feed_num - 1]
                            * self._scenario_meta_data["Tu"][
                                self._num_u + feed_num - 1
                            ],
                            label=plot_var_property.label,
                            **plt_kwargs,
                        )

                        if ax_idx in inset_axes:
                            for inset_ax in inset_axes[ax_idx]:
                                inset_ax.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._tvp[:, feed_num - 1]
                                    * self._scenario_meta_data["Tu"][
                                        self._num_u + feed_num - 1
                                    ],
                                    label=plot_var_property.label,
                                    **plt_kwargs,
                                )

                                ax.indicate_inset_zoom(
                                    inset_ax,
                                    edgecolor="black",
                                    linewidth=1.0,
                                )

                    elif plot_var_name[0] != "u" and plot_var_name[0] != "x":
                        if plot_var_name in self._scenario_meta_data["aux_var_names"]:
                            aux_expression_idx = np.where(
                                np.array(self._scenario_meta_data["aux_var_names"])
                                == plot_var_name
                            )[0][0]
                            ax.plot(
                                self._data_simulator._time,
                                self._data_simulator._aux[:, aux_expression_idx],
                                label=plot_var_property.label,
                                **plt_kwargs,
                            )

                            if ax_idx in inset_axes:
                                for inset_ax in inset_axes[ax_idx]:
                                    inset_ax.plot(
                                        self._data_simulator._time,
                                        self._data_simulator._aux[
                                            :, aux_expression_idx
                                        ],
                                        label=plot_var_property.label,
                                        **plt_kwargs,
                                    )

                                    ax.indicate_inset_zoom(
                                        inset_ax, edgecolor="black", linewidth=1.0
                                    )

                        else:
                            self._graphics_mpc.add_line(
                                var_type=f"_tvp",
                                var_name=plot_var_name,
                                axis=ax,
                                label=plot_var_property.label,
                                **plt_kwargs,
                            )

                            if ax_idx in inset_axes:
                                for inset_ax in inset_axes[ax_idx]:
                                    self._graphics_mpc.add_line(
                                        var_type=f"_tvp",
                                        var_name=plot_var_name,
                                        axis=inset_ax,
                                        label=plot_var_property.label,
                                        **plt_kwargs,
                                    )

                                    ax.indicate_inset_zoom(
                                        inset_ax,
                                        edgecolor="black",
                                        linewidth=1.0,
                                    )

                if constraints is not None:
                    for constraint in constraints:
                        if constraint.color is None:
                            color = "black"
                        else:
                            color = constraint.color

                        if constraint.ax_idx == ax_idx:
                            # make white line practically invisible:
                            if constraint.color == "white":
                                ax.hlines(
                                    constraint.value,
                                    0.0,
                                    self._scenario_meta_data["n_days_mpc"],
                                    color=color,
                                    linestyle="-",
                                    linewidth=0.1,
                                )
                            else:
                                ax.hlines(
                                    constraint.value,
                                    0.0,
                                    self._scenario_meta_data["n_days_mpc"],
                                    color=color,
                                    linestyle="--",
                                )

                            if ax_idx in inset_axes:
                                for inset_ax in inset_axes[ax_idx]:
                                    inset_ax.hlines(
                                        constraint.value,
                                        0.0,
                                        self._scenario_meta_data["n_days_mpc"],
                                        color=color,
                                        linestyle="--",
                                    )

                                    ax.indicate_inset_zoom(
                                        inset_ax, edgecolor="black", linewidth=1.0
                                    )

            for ax in axes_stacked:
                if plot_var_name[0] == "u":
                    loc = 2
                else:
                    loc = 0
                labels = [line.get_label() for line in ax.get_lines()]
                if not labels[0].startswith("_"):
                    temp_legend = ax.legend(
                        ncol=max(1, len(labels) // 3), loc="upper left"
                    )
                    temp_legend.remove()
                    ax_inputs_liquid.legend(
                        ncol=max(1, len(labels) // 3), loc="upper right"
                    )
                    if "silage" in labels[0] and "Scenario_2" in plot_save_name:
                        ax_inputs_liquid.add_artist(temp_legend)
                ax.grid(True, linestyle="--")

            # axis.yaxis.set_label_coords(-0.1, 0)
            axis.set_ylabel(
                y_label,
                rotation=0,
                labelpad=30.0,
            )

            if plot_inputs:
                # ax_inputs_liquid.yaxis.set_label_coords(0.1, 0)
                ax_inputs_liquid.set_ylabel(
                    r"$\dot V_{feed,manure}$" + "\n" + r"$[m^3/d]$",
                    rotation=0,
                    labelpad=30.0,
                )

        # Gray coloring of plot background
        if isinstance(color_background_indices, Iterable):
            x = self._data_simulator._time
            y = self._data_simulator["_tvp", "v_ch4_dot_tank_out"]

            mask = y > 0
            if mask[0]:
                start = x[0]
            for i in range(1, len(x)):
                if mask[i] and not mask[i - 1]:  # Start of a region
                    start = x[i]
                if not mask[i] and mask[i - 1]:  # End of a region
                    end = x[i]

                    for ax_idx in color_background_indices:
                        axes[ax_idx].axvspan(start[0], end[0], color="gray", alpha=0.3)

        if time_range is not None:
            plt.setp(axes, xlim=time_range)

        fig.set_size_inches(w=8, h=2 * len(axes))

        fig.tight_layout()
        if plot_save_name is not None:
            fig.savefig(
                fname=str(
                    Path(self.result_directory, "plots", f"{plot_save_name}.png")
                ),
                dpi=dpi,
                format="png",
            )
            # Save plot to a pickle so we can add the plant output later
            with open(
                Path(self.result_directory, "plots", f"{plot_save_name}.pkl"), "wb"
            ) as file:
                pickle.dump(fig, file)
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
