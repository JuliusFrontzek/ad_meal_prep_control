from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np

plt.rcParams["axes.grid"] = True


@dataclass_json
@dataclass(kw_only=True)
class MPLProperties:
    color: str
    linewidth: float
    linestyle: str


@dataclass_json
@dataclass(kw_only=True)
class PlotVarProperty:
    mpl_properties: MPLProperties
    label: str = None


@dataclass
class PostProcessing:
    scenario_name: str
    _default_mpl_properties = MPLProperties(color="blue", linewidth=2.0, linestyle="-")
    _default_plot_var_properties = PlotVarProperty(
        mpl_properties=_default_mpl_properties
    )

    def __post_init__(self):
        self._data_simulation = do_mpc.data.load_results(
            f"./results/{self.scenario_name}_mpc_results.pkl"
        )
        self._data_simulator = self._data_simulation["simulator"]
        self._data_mpc = self._data_simulation["mpc"]
        with open(f"./results/{self.scenario_name}_scenario_meta_data.pkl", "rb") as fp:
            self._scenario_meta_data = pickle.load(fp)

        self._num_u = len(self._scenario_meta_data["sub_names"])
        try:
            self._num_dictated_subs = len(
                self._scenario_meta_data["disturbances"]["dictated_feeding"]
            )
        except TypeError:
            self._num_dictated_subs = 0

        self._graphics_mpc = do_mpc.graphics.Graphics(self._data_mpc)

    def plot(
        self,
        subplot_labels_and_vars: list[tuple[str, dict[str, PlotVarProperty]]],
        plot_inputs: bool = True,
        time_range: tuple[float] = None,
    ):
        if plot_inputs:
            subplot_labels_and_vars.insert(0, ("u_norm", {f"u_norm": None}))
        fig, ax = plt.subplots(len(subplot_labels_and_vars), sharex=True)

        ax[-1].set_xlabel("Time [d]")

        for ax, subplot_label_and_vars in zip(ax, subplot_labels_and_vars):
            y_label, plot_var_properties = subplot_label_and_vars

            labels = []
            for plot_var_name, plot_var_property in plot_var_properties.items():
                plt_kwargs = (
                    plot_var_property.mpl_properties.to_dict()
                    if isinstance(plot_var_property, PlotVarProperty)
                    else self._default_mpl_properties.to_dict()
                )
                if plot_var_name[0] == "u":
                    plt_kwargs.pop("color")

                plot_var_type = (
                    plot_var_name[0] if plot_var_name[0] in ["x", "u", "tvp"] else None
                )
                if plot_var_type is not None:
                    self._graphics_mpc.add_line(
                        var_type=f"_{plot_var_type}",
                        var_name=plot_var_name,
                        axis=ax,
                        **plt_kwargs,
                    )
                    if plot_var_type == "u":
                        labels = [
                            sub.lower().replace("_", " ")
                            for sub in self._scenario_meta_data["sub_names"]
                        ]
                    elif plot_var_type == "x":
                        x_num = int(plot_var_name.split("_")[-1])
                        labels.append(
                            self._scenario_meta_data["_state_names"][x_num - 1]
                        )
                else:
                    if plot_var_name[0] == "y":
                        y_num = int(plot_var_name.split("_")[-1])

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._y[
                                :, self._num_u + y_num + self._num_dictated_subs - 1
                            ],
                            **plt_kwargs,
                        )
                        labels.append(
                            self._scenario_meta_data["_meas_names"][y_num - 1]
                        )

                    elif plot_var_name.startswith("dictated_sub_feed"):
                        feed_num = int(plot_var_name.split("_")[-1])

                        plt_kwargs["drawstyle"] = "steps-post"

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._tvp[:, feed_num - 1],
                            **plt_kwargs,
                        )
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
                            labels.append(plot_var_name)
                        else:
                            self._graphics_mpc.add_line(
                                var_type=f"_tvp",
                                var_name=plot_var_name,
                                axis=ax,
                                **plt_kwargs,
                            )
                            if plot_var_property.label is None:
                                labels.append(plot_var_name)
                            else:
                                labels.append(plot_var_property.label)
            ax.legend(labels=labels)
            ax.set_ylabel(y_label)

        if time_range is not None:
            plt.setp(ax, xlim=time_range)
        plt.show()


if __name__ == "__main__":
    default_mpl_properties = MPLProperties(color="red", linewidth=3.0, linestyle="-")
    default_plot_property = PlotVarProperty(
        mpl_properties=default_mpl_properties, label=None
    )
    post_processing = PostProcessing("Scenario_1a_weird_end")
    post_processing.plot(
        [
            (
                "States",
                {
                    # "x_19": PlotVarProperty(
                    #     mpl_properties=MPLProperties(
                    #         color="blue", linewidth=1.0, linestyle="-"
                    #     )
                    # ),
                    # "x_20": PlotVarProperty(
                    #     mpl_properties=MPLProperties(
                    #         color="green", linewidth=1.0, linestyle="-"
                    #     )
                    # ),
                    "x_2": None,
                    "x_3": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="green", linewidth=1.0, linestyle="-"
                        )
                    ),
                },
            ),
            ("V'g", {"y_1": default_plot_property}),
            # (
            #     "Forced substrates",
            #     {"dictated_sub_feed_1": None, "dictated_sub_feed_2": None},
            # ),
            (
                "Volume flow ch4",
                {
                    "v_ch4_dot_tank_in": None,
                    "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange", linewidth=1.0, linestyle="-."
                        )
                    ),
                },
            ),
            ("pH", {"y_4": default_plot_property}),
        ]
    )
