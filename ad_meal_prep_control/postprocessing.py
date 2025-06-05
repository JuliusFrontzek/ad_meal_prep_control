from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np
from pathlib import Path
import math
from ad_meal_prep_control import (
    substrates,
    params_R3)
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
        "MAIZE_SILAGE": "orange",
        "GRASS_SILAGE": "limegreen",
        "CATTLE_MANURE": "sienna",
        "SUGAR_BEET_SILAGE": "deeppink",
        "SUGAR_BEET_SILAGE_VERY_UNCERTAIN": "saddlebrown",
    }

    # substrates' linewidth specs:
    _LINEWIDTH_sub_min = 0.5
    _LINEWIDTH_sub_max = 2.5

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

        # compute organic loading rate (OLR):
        # get absolute volume flows
        u_controlled = self._data_simulator._u * self._scenario_meta_data["Tu"][0:self._num_u]
        if self._num_dictated_subs > 0:
            u_dictated = self._data_simulator._tvp[:, 0] * self._scenario_meta_data["Tu"][self._num_u]

        # get inlet concentrations:
        self._xis_controlled_subs = []
        for sub_name in self._scenario_meta_data["sub_names"]:
            self._xis_controlled_subs.append(getattr(substrates, sub_name).xi)

        if self._num_dictated_subs > 0:
            self._xis_dictated_subs = []
            dictated_sub_names = list(self._scenario_meta_data["disturbances"]["dictated_feeding"].keys())
            for sub_name in dictated_sub_names:
                self._xis_dictated_subs.append(getattr(substrates, sub_name).xi)

        # aggregate all volatile solids components
        self._S_vs_controlled = [sum(self._xis_controlled_subs[k][0:4] + self._xis_controlled_subs[k][5:11]) for k
                                 in range(self._num_u)]
        self._olr_controlled = np.sum((np.array(self._S_vs_controlled) * u_controlled)/params_R3.Vl,1)  # sum of all substrates

        if self._num_dictated_subs > 0:
            self._S_vs_dictated = [sum(self._xis_dictated_subs[k][0:4] + self._xis_controlled_subs[k][5:11]) for k
                                   in range(self._num_dictated_subs)]
            if self._num_dictated_subs > 1:
                self._olr_dictated = np.sum((np.array(self._S_vs_dictated) * u_dictated) / params_R3.Vl, 1)
            else:
                self._olr_dictated = (np.array(self._S_vs_dictated) * u_dictated) / params_R3.Vl

        # compute OLR:
        if self._num_dictated_subs > 0:
            self.olr = self._olr_controlled + self._olr_dictated
        else:
            self.olr = self._olr_controlled

        # get daily OLR by averaging over batches of 24h (thanks chatGPT):
        time = np.array(self._data_simulator._time[:, 0])
        day_indices = np.floor(time).astype(int)  # Convert time in days to integer day indices
        unique_days = np.unique(day_indices)
        olr_daily_average = np.array([self.olr[day_indices == day].mean() for day in unique_days])
        self.unique_days = unique_days
        self.olr_daily = olr_daily_average

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
            input_inset_axes: list[dict[str, tuple]] = None,
            other_inset_axes: list[dict[str, tuple]] = None,
            color_background_indices: tuple[int] = None,
            plot_olr: bool = False,  # __SH: if True, height_ratios must have 1 additional list entry
    ):

        # insert volume flows and/or OLR to subplots:
        if plot_inputs:
            subplot_labels_and_vars.insert(
                0, (r"$\dot V_{feed,silage}$" + "\n" + r"$[m^3/d]$", {f"u_norm": None})
            )

        if plot_olr and self._num_dictated_subs == 0:
            subplot_labels_and_vars.insert(
                1, (r"$OLR$" + "\n" + r"$[kg VS/m^3/d]$", {f"OLR": None})
            )

        if plot_olr and self._num_dictated_subs > 0:
            subplot_labels_and_vars.insert(
                2, (r"$OLR$" + "\n" + r"$[kg VS/m^3/d]$", {f"OLR": None})
            )

        if height_ratios is None:
            height_ratios = [1 for _ in range(len(subplot_labels_and_vars))]

        num_sub_plots = len(height_ratios)
        fig, axes = plt.subplots(
            num_sub_plots,
            sharex=True,
            gridspec_kw={"height_ratios": height_ratios},
        )

        # create inset axes for controlled substrates:
        if plot_inputs:
            ax_inputs_liquid = axes[0].twinx()

            inset_axes_u = []
            inset_axes_u_liq = []
            if input_inset_axes is not None:
                for input_inset_ax in input_inset_axes:
                    # unpack elements "days" and "ylimit" of input_inset_ax and save:
                    x1, x2, y1, y2 = (
                        *input_inset_ax["days"],
                        *input_inset_ax["ylimit"],
                    )

                    _inset_ax_input = axes[0].inset_axes(
                        input_inset_ax["inset_axis_specs"],
                        xlim=(x1, x2),
                        ylim=(y1, y2),
                        # xticklabels=[],
                        # yticklabels=[],
                    )
                    _inset_ax_input.grid(True, linestyle="--")
                    _inset_ax_input_liq = _inset_ax_input.twinx()
                    inset_axes_u.append(_inset_ax_input)
                    inset_axes_u_liq.append(_inset_ax_input_liq)

            # create all other inset axes:
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

        # draw the actual plots:
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
                    # plot states:
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

                    # plot outputs (includes dictated substrates):
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
                                    :, self._num_u + y_num + self._num_dictated_subs - 1,
                                    ],
                                    label=plot_var_property.label,
                                    **plt_kwargs,
                                )

                                ax.indicate_inset_zoom(
                                    inset_ax,
                                    edgecolor="black",
                                    linewidth=1.0,
                                )

                    # plot inputs (controlled substrates):
                    elif plot_var_name[0] == "u":
                        colors = [
                            self._SUBSTRATE_COLORS[sub_name]
                            for sub_name in self._scenario_meta_data["sub_names"]
                        ]
                        if "Scenario_1" in self.scenario_name:
                            # linearly interpolate linewidths of substrates (first thick, then thin):
                            linewidths = [(self._LINEWIDTH_sub_min - self._LINEWIDTH_sub_max) / (self._num_u - 1) * sub_k +
                                          self._LINEWIDTH_sub_max for sub_k in range(self._num_u)]
                        else:  # same std. linewidth for all substrates
                            linewidths = [1 for _ in range(self._num_u)]

                        # iterate over all substrates:
                        for i in range(self._num_u):
                            plt_kwargs["color"] = colors[i]
                            plt_kwargs["linewidth"] = linewidths[i]

                            sub_name = self._scenario_meta_data["sub_names"][i]
                            sub = getattr(substrates, sub_name)

                            if sub.state == "solid":
                                if input_inset_axes is not None:
                                    for inset_ax_u in inset_axes_u:
                                        inset_ax_u.plot(
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

                            else:  # liquid substrates
                                ax_inputs_liquid.plot(
                                    self._data_simulator._time,
                                    self._data_simulator._u[:, i]
                                    * self._scenario_meta_data["Tu"][i],
                                    label=sub_name.lower().replace("_", " "),
                                    **plt_kwargs,
                                )
                                if input_inset_axes is not None:
                                    for inset_ax_u_liq in inset_axes_u_liq:
                                        inset_ax_u_liq.plot(
                                            self._data_simulator._time,
                                            self._data_simulator._u[:, i]
                                            * self._scenario_meta_data["Tu"][i],
                                            **plt_kwargs,
                                        )
                                        # for input insets, use same ylims for left and right y-axes:
                                        inset_ax_u_liq.set_ylim(
                                            inset_axes_u[0].get_ylim()
                                        )
                                        inset_ax_u_liq.set_axis_off()  # mute right ylabels

                            if input_inset_axes is not None:
                                for inset_ax_u in inset_axes_u:
                                    ax.indicate_inset_zoom(inset_ax_u, edgecolor="black", linewidth=1.0)

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

                        # add constraints line for liquid substrates:
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

                    # plot dictated feeding: __SH: multiple dictated feeds are probably not plotted in the same subplot
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

                    # plot OLR:
                    elif plot_var_name.startswith("OLR"):

                        if plt_kwargs["color"] is None:
                            plt_kwargs["color"] = "black"

                        plt_kwargs["drawstyle"] = "steps-post"

                        # the last data point must be repeated for step plots to be shown until the end of plotted time:
                        unique_days_plot = np.append(self.unique_days, self._data_simulator._time[-1])
                        olr_daily_plot = np.append(self.olr_daily, self.olr_daily[-1])
                        ax.plot(
                            unique_days_plot,
                            olr_daily_plot,
                            **plt_kwargs,
                        )
                        ax.set_ylim(0,)
                        if max(self.olr_daily) < 6:
                            # Set y-ticks to show all integer values up to the maximum value
                            ax.set_yticks(np.arange(0, np.ceil(max(self.olr_daily)) + 1, 1))

                    # plot auxiliary variables:
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
                                        inset_ax,
                                        edgecolor="black",
                                        linewidth=1.0
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

                # __SH: add constraints:
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
                                        inset_ax,
                                        edgecolor="black",
                                        linewidth=1.0
                                    )

            # __SH: legend cosmetics:
            for ax in axes_stacked:
                # if plot_var_name[0] == "u":
                #     loc = 2
                # else:
                #     loc = 0
                labels = [line.get_label() for line in ax.get_lines()]
                if not labels[0].startswith("_"):
                    temp_legend = ax.legend(
                        ncol=max(1, len(labels) // 3), loc="upper left"
                    )
                    temp_legend.remove()
                    ax_inputs_liquid.legend(
                        ncol=max(1, len(labels) // 3), loc="upper right"
                    )
                    if "silage" in labels[0]: # and "Scenario_2" in plot_save_name: # __SH
                        ax_inputs_liquid.add_artist(temp_legend)
                ax.grid(True, linestyle="--")

            # axis.yaxis.set_label_coords(-0.1, 0)
            axis.set_ylabel(
                y_label,
                rotation=0,
                labelpad=self._scenario_meta_data["n_days_mpc"],
            )

            # __SH: adjust ylabel for liquid substrates:
            if plot_inputs:
                ax_inputs_liquid.set_ylabel(
                    r"$\dot V_{feed,manure}$" + "\n" + r"$[m^3/d]$",
                    rotation=0,
                    labelpad=self._scenario_meta_data["n_days_mpc"],
                )

        # Gray coloring of plot background
        if isinstance(color_background_indices, Iterable):
            x = self._data_simulator._time
            y = self._data_simulator["_tvp", "v_ch4_dot_tank_out"]

            mask = y > 0
            if mask[0]:  # __SH: edge case CHP on at beginning
                start = x[0]
            for i in range(1, len(x)):
                if mask[i] and not mask[i - 1]:  # Start of a region
                    start = x[i]
                if not mask[i] and mask[i - 1]:  # End of a region
                    end = x[i]
                    for ax_idx in color_background_indices:
                        axes[ax_idx].axvspan(start[0], end[0], color="gray", alpha=0.3, lw=0)
            if mask[-1] and start > end:  # __SH: edge case CHP on at end
                end = x[-1]
                for ax_idx in color_background_indices:
                    axes[ax_idx].axvspan(start[0], end[0], color="gray", alpha=0.3, lw=0)

        if time_range is not None:
            plt.setp(axes, xlim=time_range)

        # __SH: set xticks for time axis at integer multiples of 5:
        axes[-1].set_xlim(0, np.ceil(self._scenario_meta_data["n_days_mpc"]/5)*5)

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
