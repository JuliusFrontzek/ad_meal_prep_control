from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import sys

scenario_names = [
    "Scenario_1b_quadratic_nominal_feedback_mismatch_5std_3tsap",
]
plot_save_names = [
    "Scenario_1b_quadratic_nominal_feedback_mismatch_5std_3tsap",
]

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = False

for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing = PostProcessing(
        result_directory="../scenarios/results",
        scenario_name=scenario_name,
    )

    post_processing.plot(
        [
            (
                r"$\dot d_{feed}$" + "\n" + r"$[m^3/d]$",
                {
                    "dictated_sub_feed_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color=post_processing.MEASUREMENT_COLOR,
                            linewidth=1.5,
                            linestyle="-",
                        ),
                        label="cattle manure\nwith large uncertainty",
                    ),
                },
            ),
            (
                r"$\dot V$" + "\n" + r"$[m^3/d]$",
                {
                    "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange", linewidth=1.5, linestyle="-."
                        ),
                        label=r"Reference ($\dot V_{CH_4}$)",
                    ),
                    "v_ch4_dot_tank_in": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color=post_processing.MEASUREMENT_COLOR,
                            linewidth=1.5,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{CH_4}$",
                    ),
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="blue", linewidth=1.5, linestyle="-"
                        ),
                        label=r"$\dot V_g$",
                    ),
                },
            ),
            (
                r"$pH$" + "\n" + r"$[-]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                        label='$pH_{plant}$',
                    )
                },
            ),
            (
                r"$Inhibition$" + "\n" + r"$[-]$",
                {
                    "inhibition_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="blue", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_1$",
                    ),
                    "inhibition_2": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="green", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_2$",
                    ),
                    "inhibition_3": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="red", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_3$",
                    ),
                },
            ),
        ],
        plot_save_name=plot_save_name,
        time_range=(0.0, 30.0),
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[2, 1, 2, 1, 1],
        input_inset_axis={
            "days": (12, 18),
            "ylimit": (-1.0, 5.0),
            "inset_axis_specs": (0.4, 0.4, 0.4, 0.3),
        },
        other_inset_axes=[
            {
                "plot_idx": 2,
                "days": (12, 18),
                "ylimit": (345, 355),
                "inset_axis_specs": (0.37, 0.7, 0.3, 0.2),
            },
        ],
    )
