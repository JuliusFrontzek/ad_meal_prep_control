from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint
)
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_1b import controller_plotting_1b

scenario_name = ["Scenario_1b_quadratic_nominal_feedback_mismatch_2std_3tsap"]
plot_save_name = scenario_name

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = False

for scenario_name, plot_save_name in zip(scenario_name, plot_save_name):
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
                            color='saddlebrown',
                            linewidth=1,
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
                            color="grey", linewidth=1, linestyle='dashed'
                        ),
                        label=r"Reference ($\dot V_{CH_4}$)",
                    ),
                    "v_ch4_dot_tank_in": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color='blue',
                            linewidth=1,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{CH_4, plant}$",
                    ),
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="black", linewidth=1, linestyle="-"
                        ),
                        label=r"$\dot V_{g, plant}$",
                    ),
                },
            ),
            (
                r"$pH$" + "\n" + r"$[-]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1, linestyle="-", color = 'black'),
                        label='$pH_{plant}$',
                    )
                },
            ),
            (
                r"$Inhibition$" + "\n" + r"$[-]$",
                {
                    "inhibition_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="black", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_1$",
                    ),
                    "inhibition_2": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="mediumblue", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_2$",
                    ),
                    "inhibition_3": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="cornflowerblue", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_3$",
                    ),
                },
            ),
        ],
        plot_save_name=plot_save_name,
        constraints=[
            # adapt ylim of plots by adding invisible horizontal lines:
            Constraint(value=0, ax_idx=2, color="white"),  # gas production lower bound
            Constraint(value=7.3, ax_idx=3, color="white"),  # pH lower bound
            Constraint(value=7.5, ax_idx=3, color="white"),  # pH upper bound
            Constraint(value=0, ax_idx=4, color="white"),  # inhibtion lower bound
        ],
        time_range=(0.0, 30.0),
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[1, 1, 2, 1, 1],
        input_inset_axis={
            "days": (4.7, 6),
            "ylimit": (0.5, 1.3),
            "inset_axis_specs": (0.6, 0.45, 0.2, 0.2),
        },
        other_inset_axes=[
            {
                "plot_idx": 2,
                "days": (12, 18),
                "ylimit": (345, 355),
                "inset_axis_specs": (0.37, 0.7, 0.3, 0.2),
            },
            {
                "plot_idx": 2,
                "days": (2.6, 3.6),
                "ylimit": (320, 580),
                "inset_axis_specs": (0.15, 0.1, 0.1, 0.2),
            },
        ],
    )
controller_plotting_1b()