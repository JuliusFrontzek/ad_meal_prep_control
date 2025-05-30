from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_1a import controller_plotting_1a

scenario_names = ["Scenario_1a_quadratic_nominal_feedback_mismatch_3std_3tsap",
                  "Scenario_1a_quadratic_nominal_ideal_feedback_3tsap",
                  "Scenario_1a_quadratic_robust_feedback_mismatch_3std_3tsap"
                  ]
plot_save_names = ["Scenario_1a_quadratic_nominal_feedback_mismatch_3std_3tsap_substrate_costs",
                   "Scenario_1a_quadratic_nominal_ideal_feedback_3tsap_substrate_costs",
                   "Scenario_1a_quadratic_robust_feedback_mismatch_3std_3tsap_substrate_costs"
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
                r"$\dot V$" + "\n" + r"$[m^3/d]$",
                {
                    "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="grey", linewidth=1, linestyle="dashed"
                        ),
                        label=r"Reference" "\n" "($\dot V_{CH_4}$)",
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
                        mpl_properties=MPLProperties(linewidth=1,
                                                     linestyle="-",
                                                     color="black"),
                        label="$pH_{plant}$",
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
        # input_inset_axis={
        #    "days": (0, 10),
        #    "ylimit": (-1.0, 30.0),
        #    "inset_axis_specs": (0.4, 0.4, 0.27, 0.3),
        # },
        # other_inset_axes=[
        #    {
        #        "plot_idx": 1,
        #        "days": (15, 25),
        #        "ylimit": (340, 360),
        #        "inset_axis_specs": (0.4, 0.7, 0.27, 0.2),
        #    },
        # ],
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[1, 2, 1, 1],
    )
controller_plotting_1a()
