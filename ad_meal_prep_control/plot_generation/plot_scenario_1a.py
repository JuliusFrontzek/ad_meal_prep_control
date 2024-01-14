from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import matplotlib.pyplot as plt

scenario_names = ["Scenario_1a_linear", "Scenario_1a_quadratic"]
plot_save_names = [
    "scenario_1a_linear_substrate_costs",
    "Scenario_1a_quadratic_substrate_costs",
]

cmap = plt.colormaps.get_cmap("plasma_r").resampled(7).colors

for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing = PostProcessing(
        result_directory="/home/julius/Projects/ad_meal_prep_control/results",
        scenario_name=scenario_name,
    )

    post_processing.plot(
        [
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
                        mpl_properties=MPLProperties(
                            color=post_processing.MEASUREMENT_COLOR,
                            linewidth=1.5,
                            linestyle="-",
                        ),
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
