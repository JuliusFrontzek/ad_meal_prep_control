from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)

scenario_names = ["Scenario_1b_multi_stage", "Scenario_1b_non_robust"]
plot_save_names = ["scenario_1b_multi_stage", "scenario_1b_non_robust"]

for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing_multi_stage = PostProcessing(
        result_directory="/home/julius/Projects/ad_meal_prep_control/results",
        scenario_name=scenario_name,
    )

    post_processing_multi_stage.plot(
        [
            (
                r"$\dot d_{feed}$" + "\n" + r"$[m^3/d]$",
                {
                    "dictated_sub_feed_1": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                        label="cattle manure\nwith large uncertainty",
                    ),
                },
            ),
            (
                r"$\dot V_{CH_4}$" + "\n" + r"$[m^3/d]$",
                {
                    "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange", linewidth=1.5, linestyle="-."
                        ),
                        label="Reference",
                    ),
                    "v_ch4_dot_tank_in": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                        label="Actual",
                    ),
                },
            ),
            (
                r"$\dot V_g$" + "\n" + r"$[m^3/d]$",
                {
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                        label="",
                    )
                },
            ),
            (
                r"$pH$" + "\n" + r"$[1]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                        label="",
                    )
                },
            ),
        ],
        plot_save_name=plot_save_name,
        time_range=(0.0, 30.0),
    )