from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import sys

# scenario_names = ["Scenario_1a_quadratic_no_feedback_mismatch_1std_ch",
#                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr',
#                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_li']
# plot_save_names = ["Scenario_1a_quadratic_no_feedback_mismatch_1std_ch_substrate_costs",
#                   'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr_substrate_costs',
#                   'Scenario_1a_quadratic_no_feedback_mismatch_1std_li_substrate_costs']

# scenario_names = ["Scenario_1a_quadratic_feedback_mismatch_5std_3tsap",
#                   'Scenario_1a_quadratic_nominal_ideal_feedback_3tsap',
#                   'Scenario_1a_quadratic_robust_feedback_mismatch_5std_3tsap']
# plot_save_names = ["Scenario_1a_quadratic_feedback_mismatch_5std_3tsap_substrate_costs",
#                    'Scenario_1a_quadratic_nominal_ideal_feedback_3tsap_substrate_costs',
#                    'Scenario_1a_quadratic_robust_feedback_mismatch_5std_3tsap_substrate_costs']

scenario_names = ["Scenario_1a_quadratic"]
plot_save_names = ["Scenario_1a_quadratic_test"]

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = True


for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing = PostProcessing(
        result_directory="./results",
        scenario_name=scenario_name,
    )

    post_processing.plot(
        [
            (
                r"$pH$" + "\n" + r"$[1]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1, linestyle="-"),
                        label="$pH_{plant}$",
                    )
                },
            ),
            (
                r"$\dot V$" + "\n" + r"$[m^3/d]$",
                {
                    "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange", linewidth=1, linestyle="-."
                        ),
                        label=r"Reference" "\n" "($\dot V_{CH_4}$)",
                    ),
                    "v_ch4_dot_tank_in": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color=post_processing.MEASUREMENT_COLOR,
                            linewidth=1,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{CH_4, plant}$",
                    ),
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="blue", linewidth=1, linestyle="-"
                        ),
                        label=r"$\dot V_{g, plant}$",
                    ),
                },
            ),
            (
                r"$Inhibition$" + "\n" + r"$[1]$",
                {
                    "inhibition_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="blue", linewidth=1, linestyle="-"
                        ),
                        label=r"$Inhibition_1$",
                    ),
                    "inhibition_3": PlotVarProperty(
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
        input_inset_axis={
            "days": (15, 25),
            "ylimit": (-1.0, 20.0),
            "inset_axis_specs": (0.4, 0.4, 0.27, 0.3),
        },
        other_inset_axes=[
            {
                "plot_idx": 2,
                "days": (15, 25),
                "ylimit": (440, 460),
                "inset_axis_specs": (0.4, 0.7, 0.27, 0.2),
            },
        ],
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[2, 1, 3, 1],
    )
