from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint
)
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_1a import controller_plotting_1a

scenario_names = ['Scenario_1a_quadratic_no_feedback_mismatch_1std_ch',
                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr',
                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_li']
plot_save_names = ['Scenario_1a_quadratic_no_feedback_mismatch_1std_ch_substrate_costs',
                   'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr_substrate_costs',
                   'Scenario_1a_quadratic_no_feedback_mismatch_1std_li_substrate_costs']

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 1000
    show_plot = False


for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing = PostProcessing(
        result_directory="../scenarios/results",
        scenario_name=scenario_name,
    )

    post_processing.plot(
        [
            (
                r"$\dot V \; [m^3/d]$",
                {
                    "v_ch4_dot_AD_setpoint": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="dimgrey", linewidth=1, linestyle="-."
                        ),
                        label=r"Reference ($\dot V_{CH_4}$)",
                    ),
                    "v_dot_ch4_AD_norm_condition": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color='blue',
                            linewidth=1,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{CH_4, nom+\sigma}$",
                    ),
                    "y_1": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="black", linewidth=1, linestyle="-"
                        ),
                        label=r"$\dot V_{g, nom+\sigma}$",
                    ),
                },
            ),
            (
                r"$pH \; [-]$",
                {
                    "y_4": PlotVarProperty(
                        mpl_properties=MPLProperties(linewidth=1, linestyle="-",
                                                     color="black"),
                        label="$pH_{nom+\sigma}$",
                    )
                },
            ),
            #(
            #    r"$Inhibition \; [-]$",
            #    {
            #        "inhibition_1": PlotVarProperty(
            #            mpl_properties=MPLProperties(
            #                color="black", linewidth=1, linestyle="-"
            #            ),
            #            label=r"$Inhibition_1$",
            #        ),
            #        "inhibition_2": PlotVarProperty(
            #            mpl_properties=MPLProperties(
            #                color="mediumblue", linewidth=1, linestyle="-"
            #            ),
            #            label=r"$Inhibition_2$",
            #        ),
            #        "inhibition_3": PlotVarProperty(
            #            mpl_properties=MPLProperties(
            #                color="cornflowerblue", linewidth=1, linestyle="-"
            #            ),
            #            label=r"$Inhibition_3$",
            #        ),
            #    },
            #),
        ],
        plot_save_name=plot_save_name,
        constraints=[
            # adapt ylim of plots by adding invisible horizontal lines:
            Constraint(value=0, ax_idx=1, color="white"),    # gas production lower bound
            Constraint(value=1300, ax_idx=1, color="white"), # gas production upper bound
            Constraint(value=7.1, ax_idx=2, color="white"),  # pH lower bound
            Constraint(value=7.4, ax_idx=2, color="white"),  # pH upper bound
            # Constraint(value=0, ax_idx=5, color="white"),  # inhibition lower bound
        ],
        #input_inset_axis={
        #    "days": (0, 10),
        #    "ylimit": (-1.0, 30.0),
        #    "inset_axis_specs": (0.4, 0.4, 0.27, 0.3),
        #},
        #other_inset_axes=[
        #    {
        #        "plot_idx": 1,
        #        "days": (15, 25),
        #        "ylimit": (340, 360),
        #        "inset_axis_specs": (0.4, 0.7, 0.27, 0.2),
        #    },
        #],
        dpi=dpi,
        show_plot=show_plot,
        #height_ratios=[2, 2, 1, 1],
        height_ratios=[2, 2, 1],
    )
controller_plotting_1a(scenario_names)