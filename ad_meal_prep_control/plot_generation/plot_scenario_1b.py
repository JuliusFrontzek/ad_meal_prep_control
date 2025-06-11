from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint
)
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_1b import controller_plotting_1b

# scenario_names = ["Scenario_1b_linear_robust_feedback_mismatch_1.5std_3tsap"]
scenario_names = ["Scenario_1b_quadratic_robust_feedback_mismatch_1.5std_3tsap"]
plot_save_name = scenario_names

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 1000
    show_plot = False

for scenario_name, plot_save_name in zip(scenario_names, plot_save_name):
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
                    "v_dot_ch4_AD_norm_condition": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color='blue',
                            linewidth=1,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{CH_4}$",
                    ),
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="black",
                            linewidth=1,
                            linestyle="-"
                        ),
                        label=r"$\dot V_{g}$",
                    ),
                    "v_ch4_dot_AD_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color='dimgrey',
                            linewidth=1,
                            linestyle="-.",
                        ),
                        label=r"Reference ($\dot V_{CH_4}$)",
                    ),
                },
            ),
            (
                r"$pH$" + "\n" + r"$[-]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1, linestyle="-", color='black'),
                        label='',
                    )
                },
            ),
            # (
            #     r"$Inhibition$" + "\n" + r"$[-]$",
            #     {
            #         "inhibition_1": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="black", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$I_{N-lim}$",
            #         ),
            #         "inhibition_2": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="mediumblue", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$I_{pH}$",
            #         ),
            #         "inhibition_3": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="cornflowerblue", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$I_{NH3}$",
            #         ),
            #     },
            # ),
        ],
        plot_save_name=plot_save_name,
        constraints=[
            # adapt ylim of plots by adding invisible horizontal lines:
            Constraint(value=0, ax_idx=3, color="white"),  # gas production lower bound
            Constraint(value=7.1, ax_idx=4, color="white"),  # pH lower bound
            Constraint(value=7.4, ax_idx=4, color="white"),  # pH upper bound
            # Constraint(value=0, ax_idx=5, color="white"),  # inhibition lower bound
        ],
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[1.5, 1, 1, 1.5, 1],
        input_inset_axes=[
            {
                "days": (5.5, 6.5),
                "ylimit": (0.1, 0.3),
                "inset_axis_specs": (0.3, 0.3, 0.1, 0.3),
            },
            {
                "days": (25, 27),
                "ylimit": (0.2, 0.4),
                "inset_axis_specs": (0.85, 0.3, 0.1, 0.3),
            }],
        other_inset_axes=[
            {
                "plot_idx": 3,
                "days": (2.9, 4),
                "ylimit": (330, 570),
                "inset_axis_specs": (0.18, 0.55, 0.1, 0.2),
            },
            {
                "plot_idx": 3,
                "days": (8.9, 10),
                "ylimit": (330, 480),
                "inset_axis_specs": (0.3, 0.4, 0.1, 0.2),
            },
            {
                "plot_idx": 3,
                "days": (21, 27),
                "ylimit": (345, 355),
                "inset_axis_specs": (0.65, 0.11, 0.3, 0.15),
            },
        ],
        plot_olr=True,
    )
#controller_plotting_1b(scenario_names)
