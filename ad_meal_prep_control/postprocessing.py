from dataclasses import dataclass
from dataclasses_json import dataclass_json
import do_mpc
import matplotlib.pyplot as plt
import pickle
import numpy as np

plt.rcParams["axes.grid"] = True


@dataclass_json
@dataclass(kw_only=True)
class PlotVarProperty:
    color: str
    linewidth: float
    linestyle: str


@dataclass
class PostProcessing:
    scenario_name: str

    _default_plot_var_property = PlotVarProperty(color="blue", linewidth=2., linestyle="-")
    
    def __post_init__(self):
        self._data_simulation = do_mpc.data.load_results(f"./results/{self.scenario_name}_mpc_results.pkl")
        self._data_simulator = self._data_simulation["simulator"]
        self._data_mpc = self._data_simulation["mpc"]
        with open(f'./results/{self.scenario_name}_scenario_meta_data.pkl', 'rb') as fp:
            self._scenario_meta_data = pickle.load(fp)
        
        self._num_u = len(self._scenario_meta_data["sub_names"])
        self._num_dictated_subs = len(self._scenario_meta_data["disturbances"]["dictated_feeding"])
        
        self._graphics_mpc = do_mpc.graphics.Graphics(self._data_mpc)

    def plot(self, subplot_labels_and_vars: list[tuple[str,dict[str,PlotVarProperty]]], plot_inputs: bool = True, time_range: tuple[float] = None):
        if plot_inputs:
            subplot_labels_and_vars.insert(0, ("u_norm",{f"u_norm": None}))
        fig, ax = plt.subplots(len(subplot_labels_and_vars), sharex=True)

        ax[-1].set_xlabel("Time [d]")
        
        for ax, subplot_label_and_vars in zip(ax, subplot_labels_and_vars):
            y_label, plot_vars = subplot_label_and_vars

            labels = []
            for plot_var_name, plot_var_property in plot_vars.items():
                plt_kwargs = plot_var_property.to_dict() if isinstance(plot_var_property, PlotVarProperty) else self._default_plot_var_property.to_dict()
                if plot_var_name[0] == "u":
                    plt_kwargs.pop("color")

                plot_var_type = plot_var_name[0] if plot_var_name[0] in ["x", "u"] else None
                if plot_var_type is not None:
                    self._graphics_mpc.add_line(var_type=f"_{plot_var_type}", var_name=plot_var_name, axis=ax, **plt_kwargs)
                    if plot_var_type == "u":
                        labels=[
                                sub.lower().replace("_", " ") for sub in self._scenario_meta_data["sub_names"]
                            ]
                    elif plot_var_type == "x":
                        x_num = int(plot_var_name.split("_")[-1])
                        labels.append(self._scenario_meta_data["_state_names"][x_num - 1]
                        )
                else:
                    if plot_var_name[0] == "y":
                        y_num = int(plot_var_name.split("_")[-1])

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._y[:,self._num_u + y_num + self._num_dictated_subs - 1],
                            **plt_kwargs
                        )
                        labels.append(self._scenario_meta_data["_meas_names"][y_num - 1]
                        )

                    elif plot_var_name.startswith("dictated_sub_feed"):
                        feed_num = int(plot_var_name.split("_")[-1])

                        plt_kwargs["drawstyle"] = "steps-post"

                        ax.plot(
                            self._data_simulator._time,
                            self._data_simulator._y[
                                :, self._num_u + feed_num - 1
                            ],**plt_kwargs
                        )
                        labels.append(f"Dictated substrate num. {feed_num}")
                    elif plot_var_name[0] != "u" and plot_var_name[0] != "x":
                        try:
                            aux_expression_idx = np.where(
                                np.array(self._scenario_meta_data["aux_var_names"]) == plot_var_name
                            )[0][0]
                            ax.plot(
                                self._data_simulator._time,
                                self._data_simulator._aux[:, aux_expression_idx],
                                **plt_kwargs
                            )
                            labels.append(plot_var_name)
                        except IndexError:
                            raise ValueError(
                                f"'{plot_var_name}' was detected as an aux expression. However, it was either wrongly identified or is not defined as an aux expression in the do-mpc model."
                            )
            ax.legend(labels=labels)
            ax.set_ylabel(y_label)
        
        if time_range is not None:
            plt.setp(ax, xlim=time_range)
        plt.show()

if __name__ == "__main__":
    default_plot_property = PlotVarProperty(color="red", linewidth=3., linestyle="-")
    post_processing = PostProcessing("Methanation_test_12_12")
    post_processing.plot([("States",{"x_2": None, "x_3": PlotVarProperty(color="green", linewidth=1., linestyle="-")}), ("Test", {"y_1": default_plot_property}), ("Test2", {"dictated_sub_feed_1": None, "dictated_sub_feed_2": None}), ("Volume flow ch4", {"v_ch4_dot_tank_in": None})], time_range=(1.5,4.))
